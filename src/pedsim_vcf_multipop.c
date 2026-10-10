/*
 * pedsim_vcf_multipop.c - Multi-population VCF Simulation Pipeline
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include "pedsim_vcf_multipop.h"

/* Run a command and return exit status */
static int run_command(const char* cmd, int verbose)
{
    if (verbose) {
        fprintf(stderr, "Running: %s\n", cmd);
    }
    int ret = system(cmd);
    if (ret == -1) {
        fprintf(stderr, "Error: failed to execute command\n");
        return -1;
    }
    return WEXITSTATUS(ret);
}

/* Count founders in multi-population pedigree file */
static int count_founders(const char* ped_file)
{
    FILE* fp = fopen(ped_file, "r");
    if (!fp) return -1;

    int count = 0;
    char line[1024];

    while (fgets(line, sizeof(line), fp)) {
        /* Skip comments and empty lines */
        if (line[0] == '#' || line[0] == '\n') continue;

        char id[256], father[256], mother[256];
        int pop_id;

        /* Try 4-column format first (with population ID) */
        int n = sscanf(line, "%s %s %s %d", id, father, mother, &pop_id);
        if (n >= 3) {
            if (strcmp(father, "0") == 0 && strcmp(mother, "0") == 0) {
                count++;
            }
        }
    }

    fclose(fp);
    return count;
}

/* Build migration matrix string for command line */
static char* build_migration_string(double** migration, int n_pops)
{
    /* Estimate size: each entry is at most 10 chars, plus separators */
    int buf_size = n_pops * n_pops * 12 + n_pops * 2 + 10;
    char* str = malloc(buf_size);
    if (!str) return NULL;

    str[0] = '\0';
    char temp[32];

    for (int i = 0; i < n_pops; i++) {
        if (i > 0) strcat(str, ";");
        for (int j = 0; j < n_pops; j++) {
            if (j > 0) strcat(str, ",");
            snprintf(temp, sizeof(temp), "%.6f", migration[i][j]);
            strcat(str, temp);
        }
    }

    return str;
}

/* Build population sizes string */
static char* build_pop_sizes_string(multipop_vcf_pipeline_params* params)
{
    int buf_size = params->n_populations * 16 + 16;
    char* str = malloc(buf_size);
    if (!str) return NULL;

    str[0] = '\0';
    char temp[32];

    for (int i = 0; i < params->n_populations; i++) {
        if (i > 0) strcat(str, ",");
        snprintf(temp, sizeof(temp), "%d", params->populations[i].pop_size);
        strcat(str, temp);
    }

    return str;
}

/* Build sample sizes string */
static char* build_sample_sizes_string(multipop_vcf_pipeline_params* params)
{
    int buf_size = params->n_populations * 16 + 16;
    char* str = malloc(buf_size);
    if (!str) return NULL;

    str[0] = '\0';
    char temp[32];

    for (int i = 0; i < params->n_populations; i++) {
        if (i > 0) strcat(str, ",");
        snprintf(temp, sizeof(temp), "%d", params->populations[i].sample_size);
        strcat(str, temp);
    }

    return str;
}

/* ------------------------------------------------------------------------
 * Helpers for founder divergence, column mapping, thinning and metadata
 * ------------------------------------------------------------------------ */

typedef struct { char* name; int pop; } name_pop;

static int cmp_name_pop(const void* a, const void* b)
{
    return strcmp(((const name_pop*)a)->name, ((const name_pop*)b)->name);
}

/* Read "indiv father mother pop" lines of a multipop pedigree file.
 * Returns an array sorted by name (for bsearch); *n set to its length. */
static name_pop* read_ped_pops(const char* ped_file, int* n)
{
    FILE* fp = fopen(ped_file, "r");
    if (!fp) return NULL;
    int cap = 1024;
    name_pop* np = malloc(cap * sizeof(name_pop));
    *n = 0;
    char* line = NULL;
    size_t lcap = 0;
    char id[256], fa[256], mo[256];
    int pop;
    while (getline(&line, &lcap, fp) != -1) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%255s %255s %255s %d", id, fa, mo, &pop) != 4) continue;
        if (*n >= cap) { cap *= 2; np = realloc(np, cap * sizeof(name_pop)); }
        np[*n].name = strdup(id);
        np[*n].pop = pop;
        (*n)++;
    }
    free(line);
    fclose(fp);
    qsort(np, *n, sizeof(name_pop), cmp_name_pop);
    return np;
}

static int lookup_pop(name_pop* np, int n, const char* name)
{
    name_pop key = { (char*)name, 0 };
    name_pop* hit = bsearch(&key, np, n, sizeof(name_pop), cmp_name_pop);
    return hit ? hit->pop : -1;
}

static void free_ped_pops(name_pop* np, int n)
{
    for (int i = 0; i < n; i++) free(np[i].name);
    free(np);
}

/* Founder names in the order pedtrans lists them ("# Founders:" line).
 * vcfassemble maps founder k to founder-VCF haplotype columns 2k and 2k+1. */
static char** read_founder_order(const char* seg_file, int* n)
{
    FILE* fp = fopen(seg_file, "r");
    if (!fp) return NULL;
    char* line = NULL;
    size_t lcap = 0;
    int cap = 1024;
    char** names = malloc(cap * sizeof(char*));
    *n = 0;
    while (getline(&line, &lcap, fp) != -1) {
        if (strncmp(line, "# Founders:", 11) == 0) {
            char* save = NULL;
            for (char* tok = strtok_r(line + 11, " \t\n", &save); tok;
                 tok = strtok_r(NULL, " \t\n", &save)) {
                if (*n >= cap) { cap *= 2; names = realloc(names, cap * sizeof(char*)); }
                names[(*n)++] = strdup(tok);
            }
            break;
        }
    }
    free(line);
    fclose(fp);
    return names;
}

/* Contig length from a single-chromosome VCF header (0 if absent) */
static long vcf_contig_length(const char* vcf)
{
    FILE* fp = fopen(vcf, "r");
    if (!fp) return 0;
    char* line = NULL;
    size_t lcap = 0;
    long len = 0;
    while (getline(&line, &lcap, fp) != -1) {
        if (line[0] != '#') break;
        if (strncmp(line, "##contig", 8) == 0) {
            char* lp = strstr(line, "length=");
            if (lp) len = atol(lp + 7);
        }
    }
    free(line);
    fclose(fp);
    return len;
}

/* Merge per-chromosome founder VCFs (temp_dir/chrC.vcf, C = 1..n_chrom) into
 * out_vcf. Each contig keeps the length written by the simulator that made it
 * (vcfassemble converts positions to fractions of that length). If colmap is
 * non-NULL, output haplotype column h is input column colmap[h]. */
static int merge_founder_vcfs(const char* temp_dir, int n_chrom, const int* colmap,
                              int n_hap, const char* out_vcf)
{
    FILE* out = fopen(out_vcf, "w");
    if (!out) return 1;
    char path[600];
    fprintf(out, "##fileformat=VCFv4.2\n##source=pedsim_vcf_multipop\n");
    for (int c = 1; c <= n_chrom; c++) {
        snprintf(path, sizeof(path), "%s/chr%d.vcf", temp_dir, c);
        long len = vcf_contig_length(path);
        if (len <= 0) {
            fprintf(stderr, "Error: no contig length in %s\n", path);
            fclose(out);
            return 1;
        }
        fprintf(out, "##contig=<ID=chr%d,length=%ld>\n", c, len);
    }
    fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int h = 0; h < n_hap; h++) fprintf(out, "\thap%d", h);
    fprintf(out, "\n");

    char** cols = malloc((n_hap + 16) * sizeof(char*));
    char* line = NULL;
    size_t lcap = 0;
    for (int c = 1; c <= n_chrom; c++) {
        snprintf(path, sizeof(path), "%s/chr%d.vcf", temp_dir, c);
        FILE* in = fopen(path, "r");
        if (!in) { free(cols); free(line); fclose(out); return 1; }
        while (getline(&line, &lcap, in) != -1) {
            if (line[0] == '#') continue;
            int nc = 0;
            char* save = NULL;
            for (char* tok = strtok_r(line, "\t\n", &save); tok && nc < n_hap + 16;
                 tok = strtok_r(NULL, "\t\n", &save))
                cols[nc++] = tok;
            if (nc != 9 + n_hap) {
                fprintf(stderr, "Error: %s record has %d haplotypes, expected %d\n",
                        path, nc - 9, n_hap);
                fclose(in); free(cols); free(line); fclose(out);
                return 1;
            }
            fprintf(out, "chr%d", c);
            for (int k = 1; k < 9; k++) fprintf(out, "\t%s", cols[k]);
            for (int h = 0; h < n_hap; h++)
                fprintf(out, "\t%s", cols[9 + (colmap ? colmap[h] : h)]);
            fprintf(out, "\n");
        }
        fclose(in);
    }
    free(cols);
    free(line);
    fclose(out);
    return 0;
}

/* Keep at most k SNPs per chromosome, chosen uniformly at random among the
 * records that are polymorphic in the output samples (sites where every
 * allele is derived are dropped too). Order within a chromosome is kept. */
static int thin_vcf(const char* in_vcf, const char* out_vcf, int k, unsigned int seed)
{
    FILE* in = fopen(in_vcf, "r");
    FILE* out = fopen(out_vcf, "w");
    if (!in || !out) { if (in) fclose(in); if (out) fclose(out); return 1; }
    srand48(seed ? seed : 12345);

    char* line = NULL;
    size_t lcap = 0;
    char** keep = NULL;            /* reservoir for the current chromosome */
    long seen = 0;
    char cur_chrom[64] = "";
    long total_in = 0, total_out = 0;
    keep = calloc(k, sizeof(char*));

    /* flush reservoir: records were stored with their position for sorting */
    #define FLUSH_RES() do { \
        int m = (seen < k) ? (int)seen : k; \
        for (int a = 1; a < m; a++) { char* t = keep[a]; int b = a - 1; \
            long pt = atol(strchr(t, '\t') + 1); \
            while (b >= 0 && atol(strchr(keep[b], '\t') + 1) > pt) { keep[b + 1] = keep[b]; b--; } \
            keep[b + 1] = t; } \
        for (int a = 0; a < m; a++) { fputs(keep[a], out); free(keep[a]); keep[a] = NULL; } \
        total_out += m; seen = 0; } while (0)

    while (getline(&line, &lcap, in) != -1) {
        if (line[0] == '#') { fputs(line, out); continue; }
        char chrom[64];
        if (sscanf(line, "%63s", chrom) != 1) continue;
        /* polymorphic among samples? look at genotype alleles after column 9 */
        int n0 = 0, n1 = 0, col = 0;
        for (char* q = line; *q; q++) {
            if (*q == '\t') { col++; continue; }
            if (col >= 9) { if (*q == '0') n0++; else if (*q == '1') n1++; }
        }
        if (n0 == 0 || n1 == 0) continue;
        total_in++;
        if (strcmp(chrom, cur_chrom) != 0) {
            if (cur_chrom[0]) FLUSH_RES();
            snprintf(cur_chrom, sizeof(cur_chrom), "%s", chrom);
        }
        seen++;
        if (seen <= k) {
            keep[seen - 1] = strdup(line);
        } else {
            long r = (long)(drand48() * seen);
            if (r < k) { free(keep[r]); keep[r] = strdup(line); }
        }
    }
    if (cur_chrom[0]) FLUSH_RES();
    #undef FLUSH_RES
    free(keep);
    free(line);
    fclose(in);
    fclose(out);
    fprintf(stderr, "Thinning: kept %ld of %ld polymorphic SNPs (max %d per chromosome)\n",
            total_out, total_in, k);
    return 0;
}

/* Write "sample<TAB>popJ" for every sample in the VCF header (for BA3 -M) */
static int write_meta(const char* vcf, name_pop* np, int n_np, const char* meta_file)
{
    FILE* in = fopen(vcf, "r");
    FILE* out = fopen(meta_file, "w");
    if (!in || !out) { if (in) fclose(in); if (out) fclose(out); return 1; }
    char* line = NULL;
    size_t lcap = 0;
    while (getline(&line, &lcap, in) != -1) {
        if (strncmp(line, "#CHROM", 6) == 0) {
            int col = 0;
            char* save = NULL;
            for (char* tok = strtok_r(line, "\t\n", &save); tok;
                 tok = strtok_r(NULL, "\t\n", &save), col++) {
                if (col < 9) continue;
                int pop = lookup_pop(np, n_np, tok);
                fprintf(out, "%s\tpop%d\n", tok, pop);
            }
            break;
        }
    }
    free(line);
    fclose(in);
    fclose(out);
    return 0;
}

/* Output base name: output path without a trailing ".vcf" */
static void output_base(const char* out_vcf, char* base, size_t n)
{
    snprintf(base, n, "%s", out_vcf);
    size_t L = strlen(base);
    if (L > 4 && strcmp(base + L - 4, ".vcf") == 0) base[L - 4] = '\0';
}

int run_multipop_vcf_pipeline(multipop_vcf_pipeline_params* params)
{
    char cmd[8192];
    char* temp_dir = NULL;
    int ret = 0;
    name_pop* ped_pops = NULL;
    int n_ped = 0;
    char** founder_order = NULL;
    int n_founder_order = 0;
    int* colmap = NULL;

    /* Create temporary directory */
    char template[] = "/tmp/pedsim_vcf_multipop_XXXXXX";
    temp_dir = mkdtemp(template);
    if (!temp_dir) {
        fprintf(stderr, "Error: failed to create temp directory\n");
        return 1;
    }

    char ped_file[512], vcf_file[512], seg_file[512], truth_tmp[512];
    char assembled[512], base[1024], meta_file[1100], truth_file[1100];
    snprintf(ped_file, sizeof(ped_file), "%s/ped.txt", temp_dir);
    snprintf(vcf_file, sizeof(vcf_file), "%s/founder.vcf", temp_dir);
    snprintf(seg_file, sizeof(seg_file), "%s/segments.txt", temp_dir);
    snprintf(truth_tmp, sizeof(truth_tmp), "%s/truth.tsv", temp_dir);
    snprintf(assembled, sizeof(assembled), "%s/assembled.vcf", temp_dir);
    output_base(params->output_vcf, base, sizeof(base));
    snprintf(meta_file, sizeof(meta_file), "%s.meta.txt", base);
    snprintf(truth_file, sizeof(truth_file), "%s.truth.tsv", base);

    /* Build parameter strings */
    char* pop_sizes = build_pop_sizes_string(params);
    char* sample_sizes = build_sample_sizes_string(params);
    char* mig_str = NULL;

    if (!pop_sizes || !sample_sizes) {
        fprintf(stderr, "Error: failed to allocate parameter strings\n");
        ret = 1;
        goto cleanup;
    }

    /* Step 1: Run pedsim_multipop to generate pedigree */
    fprintf(stderr, "\n=== Step 1: Generating multi-population pedigree ===\n");
    char mig_opt[4096] = "";
    if (params->migration) {
        mig_str = build_migration_string(params->migration, params->n_populations);
        if (!mig_str) { ret = 1; goto cleanup; }
        snprintf(mig_opt, sizeof(mig_opt), " -m \"%s\"", mig_str);
    } else if (params->island_m >= 0) {
        snprintf(mig_opt, sizeof(mig_opt), " -M %.6f", params->island_m);
    }
    snprintf(cmd, sizeof(cmd),
             "./pedsim_multipop -J %d -N %s -n %s -k %d%s%s -s %u -t %s -o %s",
             params->n_populations, pop_sizes, sample_sizes, params->n_generations,
             mig_opt, params->migrate_before_sampling ? " -a" : "",
             params->seed, truth_tmp, ped_file);
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: pedsim_multipop failed\n");
        goto cleanup;
    }

    int n_founders = count_founders(ped_file);
    if (n_founders <= 0) {
        fprintf(stderr, "Error: failed to count founders\n");
        ret = 1;
        goto cleanup;
    }
    int n_haplotypes = n_founders * 2;
    fprintf(stderr, "Pedigree has %d founders (%d haplotypes) across %d populations\n",
            n_founders, n_haplotypes, params->n_populations);

    int n_chromosomes = params->n_chromosomes;
    if (n_chromosomes < 1) n_chromosomes = 1;

    /* Step 2: pedtrans (run before the founder simulation so the founder order,
     * which fixes the founder-VCF column layout, is known) */
    fprintf(stderr, "\n=== Step 2: Simulating chromosome transmission ===\n");
    snprintf(cmd, sizeof(cmd),
             "./pedtrans -c %d -r %.4f -s %u -o %s %s",
             n_chromosomes, params->rec_rate, params->seed + 2, seg_file, ped_file);
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: pedtrans failed\n");
        goto cleanup;
    }

    /* Step 3: founder haplotypes, one simulation per chromosome */
    fprintf(stderr, "\n=== Step 3: Simulating founder haplotypes (VCF) ===\n");
    fprintf(stderr, "Simulating %d chromosome(s)...\n", n_chromosomes);
    if (params->div_time >= 0) {
        /* Founders of population j are drawn from tip Pj of a star tree with
         * divergence time T; tip sizes = pedigree population sizes, root size = -E */
        ped_pops = read_ped_pops(ped_file, &n_ped);
        founder_order = read_founder_order(seg_file, &n_founder_order);
        if (!ped_pops || !founder_order || n_founder_order != n_founders) {
            fprintf(stderr, "Error: could not read founder populations "
                    "(%d founders listed, %d expected)\n", n_founder_order, n_founders);
            ret = 1;
            goto cleanup;
        }
        int J = params->n_populations;
        int* fpop = malloc(n_founders * sizeof(int));
        int* nf = calloc(J, sizeof(int));
        for (int k = 0; k < n_founders; k++) {
            fpop[k] = lookup_pop(ped_pops, n_ped, founder_order[k]);
            if (fpop[k] < 0 || fpop[k] >= J) {
                fprintf(stderr, "Error: no population for founder %s\n", founder_order[k]);
                free(fpop); free(nf); ret = 1; goto cleanup;
            }
            nf[fpop[k]]++;
        }
        int n_tips = 0;
        for (int j = 0; j < J; j++) if (nf[j] > 0) n_tips++;
        if (n_tips < 2) {
            fprintf(stderr, "Error: -T needs founders in at least two populations\n");
            free(fpop); free(nf); ret = 1; goto cleanup;
        }
        /* coalsim_msc writes haplotypes tip by tip in tree order, so the block
         * for population j starts at offset[j] */
        int* offset = calloc(J, sizeof(int));
        for (int j = 1; j < J; j++) offset[j] = offset[j - 1] + 2 * nf[j - 1];
        int* used = calloc(J, sizeof(int));
        colmap = malloc(n_haplotypes * sizeof(int));
        for (int k = 0; k < n_founders; k++) {
            int j = fpop[k];
            colmap[2 * k]     = offset[j] + used[j]++;
            colmap[2 * k + 1] = offset[j] + used[j]++;
        }
        fprintf(stderr, "Star tree, divergence time T=%g generations; founders per pop:",
                params->div_time);
        for (int j = 0; j < J; j++) fprintf(stderr, " %d", nf[j]);
        fprintf(stderr, "\n");

        for (int chr = 1; chr <= n_chromosomes; chr++) {
            char ctl[600], chr_vcf[600];
            snprintf(ctl, sizeof(ctl), "%s/msc_chr%d.ctl", temp_dir, chr);
            snprintf(chr_vcf, sizeof(chr_vcf), "%s/chr%d.vcf", temp_dir, chr);
            FILE* cf = fopen(ctl, "w");
            if (!cf) { ret = 1; break; }
            /* The MSC parser takes binary trees only, so the star tree is written
             * as a caterpillar whose internal branches are 1e-6 generations long:
             * numerically a single split at T. (Zero-length internal branches are
             * not used: coalsim_msc can process simultaneous divergences in an
             * order that delays coalescence.) */
            const double eps = 1e-6;
            int tips[MAX_POPULATIONS], nt = 0;
            for (int j = 0; j < J; j++) if (nf[j] > 0) tips[nt++] = j;
            char* tree = malloc(64 * (nt + 2) + 64);
            int len = sprintf(tree, "(P%d:%.10f#%d,P%d:%.10f#%d)", tips[0], params->div_time,
                              params->populations[tips[0]].pop_size, tips[1],
                              params->div_time, params->populations[tips[1]].pop_size);
            for (int k = 2; k < nt; k++) {
                char* prev = strdup(tree);
                len = sprintf(tree, "(%sA%d:%.10f#%.0f,P%d:%.10f#%d)", prev, k - 1, eps,
                              params->coal_pop_size, tips[k],
                              params->div_time + (k - 1) * eps,
                              params->populations[tips[k]].pop_size);
                free(prev);
            }
            (void)len;
            fprintf(cf, "species_tree: %sROOT#%.0f;\nsamples: ", tree, params->coal_pop_size);
            free(tree);
            int first = 1;
            for (int j = 0; j < J; j++) {
                if (nf[j] == 0) continue;
                fprintf(cf, "%sP%d:%d", first ? "" : ",", j, 2 * nf[j]);
                first = 0;
            }
            fprintf(cf, "\nrecombination_rate: %.6f\nmutation_rate: %.6f\nseed: %u\n"
                        "output_vcf: true\nvcf_file: %s\n",
                    params->rec_rate, params->mut_rate, params->seed + chr, chr_vcf);
            fclose(cf);
            snprintf(cmd, sizeof(cmd), "./coalsim_msc -q %s%s", ctl,
                     params->verbose ? "" : " > /dev/null");
            ret = run_command(cmd, params->verbose);
            if (ret != 0) {
                fprintf(stderr, "Error: coalsim_msc failed for chromosome %d\n", chr);
                break;
            }
        }
        free(fpop); free(nf); free(offset); free(used);
        if (ret != 0) goto cleanup;
    } else {
        for (int chr = 1; chr <= n_chromosomes; chr++) {
            char chr_vcf[600];
            snprintf(chr_vcf, sizeof(chr_vcf), "%s/chr%d.vcf", temp_dir, chr);
            snprintf(cmd, sizeof(cmd),
                     "./coalsim -c %d -N %.0f -r %.4f -m %.4f -s %u -u b -V %s",
                     n_haplotypes, params->coal_pop_size, params->rec_rate,
                     params->mut_rate, params->seed + chr, chr_vcf);
            ret = run_command(cmd, params->verbose);
            if (ret != 0) {
                fprintf(stderr, "Error: coalsim failed for chromosome %d\n", chr);
                goto cleanup;
            }
        }
    }
    ret = merge_founder_vcfs(temp_dir, n_chromosomes, colmap, n_haplotypes, vcf_file);
    if (ret != 0) {
        fprintf(stderr, "Error: failed to merge founder VCFs\n");
        goto cleanup;
    }

    /* Step 4: Run vcfassemble to create sample VCF */
    fprintf(stderr, "\n=== Step 4: Assembling sample VCF ===\n");
    const char* asm_out = params->snps_per_chrom > 0 ? assembled : params->output_vcf;
    snprintf(cmd, sizeof(cmd),
             "./vcfassemble -v %s -p %s -o %s%s",
             vcf_file, seg_file, asm_out, params->samples_only ? " -s" : "");
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: vcfassemble failed\n");
        goto cleanup;
    }

    /* Step 5: optional thinning */
    if (params->snps_per_chrom > 0) {
        fprintf(stderr, "\n=== Step 5: Thinning to %d SNPs per chromosome ===\n",
                params->snps_per_chrom);
        ret = thin_vcf(assembled, params->output_vcf, params->snps_per_chrom,
                       params->seed + 3);
        if (ret != 0) {
            fprintf(stderr, "Error: thinning failed\n");
            goto cleanup;
        }
    }

    /* Metadata (for BA3 -V/-M) and migrant-ancestry truth */
    if (!ped_pops) ped_pops = read_ped_pops(ped_file, &n_ped);
    if (!ped_pops || write_meta(params->output_vcf, ped_pops, n_ped, meta_file) != 0) {
        fprintf(stderr, "Error: failed to write metadata file\n");
        ret = 1;
        goto cleanup;
    }
    snprintf(cmd, sizeof(cmd), "cp %s %s", truth_tmp, truth_file);
    ret = run_command(cmd, params->verbose);
    if (ret != 0) goto cleanup;

    fprintf(stderr, "\n=== Pipeline complete ===\n");
    fprintf(stderr, "Output VCF: %s\n", params->output_vcf);
    fprintf(stderr, "Population metadata: %s\n", meta_file);
    fprintf(stderr, "Migrant-ancestry truth: %s\n", truth_file);

    int total_samples = 0;
    for (int i = 0; i < params->n_populations; i++) {
        total_samples += params->populations[i].sample_size;
    }
    fprintf(stderr, "Total samples: %d across %d populations\n",
            total_samples, params->n_populations);

cleanup:
    free(pop_sizes);
    free(sample_sizes);
    free(mig_str);
    free(colmap);
    if (ped_pops) free_ped_pops(ped_pops, n_ped);
    if (founder_order) {
        for (int k = 0; k < n_founder_order; k++) free(founder_order[k]);
        free(founder_order);
    }

    /* Clean up temporary files */
    if (!params->keep_temp && temp_dir) {
        snprintf(cmd, sizeof(cmd), "rm -rf %s", temp_dir);
        int cleanup_ret = system(cmd);
        (void)cleanup_ret;
    } else if (params->keep_temp) {
        fprintf(stderr, "Temp files kept in: %s\n", temp_dir);
    }

    return ret;
}
