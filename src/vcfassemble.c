#include <ctype.h>
#include "vcfassemble.h"

/* Parse chromosome ID from name like "chr1", "chr2", etc. */
static int parse_chrom_id(const char* chrom_name)
{
    if (strncmp(chrom_name, "chr", 3) == 0) {
        return atoi(chrom_name + 3);
    }
    /* Try parsing as plain number */
    return atoi(chrom_name);
}

/* Read founder VCF file (supports multiple chromosomes) */
founder_vcf* read_founder_vcf(const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open VCF file %s\n", filename);
        return NULL;
    }

    founder_vcf* fvcf = malloc(sizeof(founder_vcf));
    fvcf->mutations = malloc(INITIAL_CAPACITY * sizeof(vcf_mutation));
    fvcf->n_mutations = 0;
    fvcf->capacity = INITIAL_CAPACITY;
    fvcf->n_chromosomes = 0;
    fvcf->chrom_lengths = malloc(64 * sizeof(long));  /* Support up to 64 chromosomes */
    memset(fvcf->chrom_lengths, 0, 64 * sizeof(long));
    fvcf->founder_names = NULL;
    fvcf->n_founder_samples = 0;

    char line[MAX_LINE_LENGTH];

    /* Parse header */
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "##contig", 8) == 0) {
            /* Parse chromosome ID and length: ##contig=<ID=chr1,length=10000000> */
            char* id_ptr = strstr(line, "ID=");
            char* len_ptr = strstr(line, "length=");
            if (id_ptr && len_ptr) {
                char chrom_name[64];
                sscanf(id_ptr + 3, "%[^,>]", chrom_name);
                int chrom_id = parse_chrom_id(chrom_name);
                long length = atol(len_ptr + 7);
                if (chrom_id > 0 && chrom_id <= 64) {
                    fvcf->chrom_lengths[chrom_id] = length;
                    if (chrom_id > fvcf->n_chromosomes) {
                        fvcf->n_chromosomes = chrom_id;
                    }
                }
            }
        }
        else if (line[0] == '#' && line[1] != '#') {
            char* tok = strtok(line, "\t\n");
            int col = 0;
            int sample_capacity = 64;
            fvcf->founder_names = malloc(sample_capacity * sizeof(char*));

            while (tok) {
                if (col >= 9) {
                    if (fvcf->n_founder_samples >= sample_capacity) {
                        sample_capacity *= 2;
                        fvcf->founder_names = realloc(fvcf->founder_names,
                                                       sample_capacity * sizeof(char*));
                    }
                    fvcf->founder_names[fvcf->n_founder_samples] = strdup(tok);
                    fvcf->n_founder_samples++;
                }
                tok = strtok(NULL, "\t\n");
                col++;
            }
            break;
        }
    }

    /* Default to 1 chromosome if none found in header */
    if (fvcf->n_chromosomes == 0) {
        fvcf->n_chromosomes = 1;
    }

    /* Parse mutation records */
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;

        if (fvcf->n_mutations >= fvcf->capacity) {
            fvcf->capacity *= 2;
            fvcf->mutations = realloc(fvcf->mutations,
                                       fvcf->capacity * sizeof(vcf_mutation));
        }

        vcf_mutation* mut = &fvcf->mutations[fvcf->n_mutations];
        mut->n_founders = fvcf->n_founder_samples;
        mut->founder_geno = malloc(mut->n_founders * sizeof(int));
        mut->chrom_id = 1;  /* Default */

        char* tok = strtok(line, "\t");
        int col = 0;
        int sample_idx = 0;

        while (tok) {
            switch (col) {
                case 0:
                    mut->chrom_id = parse_chrom_id(tok);
                    break;
                case 1:
                    mut->position = atol(tok);
                    break;
                case 3:
                    mut->ref = tok[0];
                    break;
                case 4:
                    mut->alt = tok[0];
                    break;
                default:
                    if (col >= 9) {
                        mut->founder_geno[sample_idx] = atoi(tok);
                        sample_idx++;
                    }
                    break;
            }
            tok = strtok(NULL, "\t\n");
            col++;
        }

        fvcf->n_mutations++;
    }

    fclose(fp);
    return fvcf;
}

void free_founder_vcf(founder_vcf* fvcf)
{
    if (!fvcf) return;

    for (int i = 0; i < fvcf->n_mutations; i++) {
        free(fvcf->mutations[i].founder_geno);
    }
    free(fvcf->mutations);

    for (int i = 0; i < fvcf->n_founder_samples; i++) {
        free(fvcf->founder_names[i]);
    }
    free(fvcf->founder_names);
    free(fvcf->chrom_lengths);

    free(fvcf);
}

/* Founder name to VCF index mapping */
typedef struct {
    char** names;
    int n_founders;
} founder_map;

static founder_map* parse_founder_list(FILE* fp)
{
    char line[MAX_LINE_LENGTH];
    founder_map* fm = malloc(sizeof(founder_map));
    fm->names = malloc(1024 * sizeof(char*));
    fm->n_founders = 0;

    rewind(fp);
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "# Founders:", 11) == 0) {
            char* tok = strtok(line + 11, " \t\n");
            while (tok && fm->n_founders < 1024) {
                fm->names[fm->n_founders] = strdup(tok);
                fm->n_founders++;
                tok = strtok(NULL, " \t\n");
            }
            break;
        }
    }
    rewind(fp);
    return fm;
}

static void free_founder_map(founder_map* fm)
{
    if (!fm) return;
    for (int i = 0; i < fm->n_founders; i++) {
        free(fm->names[i]);
    }
    free(fm->names);
    free(fm);
}

static int get_vcf_haplotype_idx(founder_map* fm, const char* founder_hom)
{
    char founder[MAX_NAME_LENGTH];
    char homolog[16];

    if (sscanf(founder_hom, "%[^:]:%s", founder, homolog) != 2) {
        return -1;
    }

    int founder_idx = -1;
    for (int i = 0; i < fm->n_founders; i++) {
        if (strcmp(fm->names[i], founder) == 0) {
            founder_idx = i;
            break;
        }
    }

    if (founder_idx < 0) {
        return -1;
    }

    int hom = (strcmp(homolog, "mat") == 0) ? 1 : 0;
    return founder_idx * 2 + hom;
}

pedtrans_segments* read_pedtrans_segments(const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open pedtrans file %s\n", filename);
        return NULL;
    }

    founder_map* fm = parse_founder_list(fp);

    pedtrans_segments* pts = malloc(sizeof(pedtrans_segments));
    pts->chroms = malloc(INITIAL_CAPACITY * sizeof(sample_chrom));
    pts->n_chroms = 0;
    pts->capacity = INITIAL_CAPACITY;
    pts->n_chromosomes = 1;  /* Default to 1 */
    pts->sample_names = malloc(INITIAL_CAPACITY * sizeof(char*));
    pts->n_samples = 0;
    pts->founder_map = fm;

    char line[MAX_LINE_LENGTH];
    char current_individual[MAX_NAME_LENGTH] = "";
    int current_chrom_id = 1;
    int current_homolog = -1;
    sample_chrom* current = NULL;
    int sample_capacity = INITIAL_CAPACITY;

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') {
            /* Check for chromosome count in header */
            if (strncmp(line, "# Chromosomes:", 14) == 0) {
                int n_chrom;
                if (sscanf(line + 14, " %d", &n_chrom) == 1) {
                    pts->n_chromosomes = n_chrom;
                }
            }
            continue;
        }
        if (line[0] == '\n') continue;

        if (strncmp(line, "Individual:", 11) == 0) {
            sscanf(line + 11, " %s", current_individual);
            current_chrom_id = 1;  /* Reset to first chromosome */
            current_homolog = -1;

            int found = 0;
            for (int i = 0; i < pts->n_samples; i++) {
                if (strcmp(pts->sample_names[i], current_individual) == 0) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                if (pts->n_samples >= sample_capacity) {
                    sample_capacity *= 2;
                    pts->sample_names = realloc(pts->sample_names,
                                                sample_capacity * sizeof(char*));
                }
                pts->sample_names[pts->n_samples] = strdup(current_individual);
                pts->n_samples++;
            }
            continue;
        }

        char* trimmed = line;
        while (*trimmed == ' ') trimmed++;

        /* Parse "Chromosome: N" lines */
        if (strncmp(trimmed, "Chromosome:", 11) == 0) {
            int chrom_num;
            if (sscanf(trimmed + 11, " %d", &chrom_num) == 1) {
                current_chrom_id = chrom_num;
                if (chrom_num > pts->n_chromosomes) {
                    pts->n_chromosomes = chrom_num;
                }
            }
            current_homolog = -1;  /* Reset homolog for new chromosome */
            continue;
        }

        if (strncmp(trimmed, "Paternal:", 9) == 0) {
            current_homolog = 0;

            if (pts->n_chroms >= pts->capacity) {
                pts->capacity *= 2;
                pts->chroms = realloc(pts->chroms,
                                      pts->capacity * sizeof(sample_chrom));
            }
            current = &pts->chroms[pts->n_chroms];
            strncpy(current->name, current_individual, MAX_NAME_LENGTH - 1);
            current->name[MAX_NAME_LENGTH - 1] = '\0';
            current->chrom_id = current_chrom_id;
            current->homolog = 0;
            current->segments = malloc(64 * sizeof(segment));
            current->n_segments = 0;
            current->capacity = 64;
            pts->n_chroms++;
            continue;
        }

        if (strncmp(trimmed, "Maternal:", 9) == 0) {
            current_homolog = 1;

            if (pts->n_chroms >= pts->capacity) {
                pts->capacity *= 2;
                pts->chroms = realloc(pts->chroms,
                                      pts->capacity * sizeof(sample_chrom));
            }
            current = &pts->chroms[pts->n_chroms];
            strncpy(current->name, current_individual, MAX_NAME_LENGTH - 1);
            current->name[MAX_NAME_LENGTH - 1] = '\0';
            current->chrom_id = current_chrom_id;
            current->homolog = 1;
            current->segments = malloc(64 * sizeof(segment));
            current->n_segments = 0;
            current->capacity = 64;
            pts->n_chroms++;
            continue;
        }

        if (current && current_homolog >= 0) {
            double start, end;
            char founder_hom[MAX_NAME_LENGTH];

            if (sscanf(trimmed, "%lf %lf %s", &start, &end, founder_hom) == 3) {
                if (current->n_segments >= current->capacity) {
                    current->capacity *= 2;
                    current->segments = realloc(current->segments,
                                                current->capacity * sizeof(segment));
                }

                segment* seg = &current->segments[current->n_segments];
                seg->start = start;
                seg->end = end;
                strncpy(seg->founder, founder_hom, MAX_NAME_LENGTH - 1);
                seg->founder[MAX_NAME_LENGTH - 1] = '\0';
                seg->homolog = get_vcf_haplotype_idx(fm, founder_hom);

                current->n_segments++;
            }
        }
    }

    fclose(fp);
    return pts;
}

void free_pedtrans_segments(pedtrans_segments* pts)
{
    if (!pts) return;

    for (int i = 0; i < pts->n_chroms; i++) {
        free(pts->chroms[i].segments);
    }
    free(pts->chroms);

    for (int i = 0; i < pts->n_samples; i++) {
        free(pts->sample_names[i]);
    }
    free(pts->sample_names);

    if (pts->founder_map) {
        free_founder_map(pts->founder_map);
    }

    free(pts);
}

int get_founder_haplotype_idx(founder_vcf* fvcf, const char* founder, int homolog)
{
    (void)fvcf;
    (void)founder;
    (void)homolog;
    return -1;
}

/*
 * ============================================================================
 * OPTIMIZED ASSEMBLY - O(1) sample lookup + binary search for segments
 * ============================================================================
 */

/*
 * Fast lookup table for multi-chromosome support:
 * entries[sample_idx * n_chromosomes + chrom_idx] -> (paternal_chrom_idx, maternal_chrom_idx)
 */
typedef struct {
    int pat_idx;  /* index into pts->chroms for paternal, or -1 */
    int mat_idx;  /* index into pts->chroms for maternal, or -1 */
} sample_lookup_entry;

typedef struct {
    sample_lookup_entry* entries;  /* indexed by: sample_idx * n_chromosomes + (chrom_id - 1) */
    int n_samples;
    int n_chromosomes;
} sample_lookup_table;

/* Build fast lookup table for samples (supports multiple chromosomes) */
static sample_lookup_table* build_sample_lookup(pedtrans_segments* pts,
                                                  char** output_samples,
                                                  int n_output)
{
    int n_chromosomes = pts->n_chromosomes;

    sample_lookup_table* lut = malloc(sizeof(sample_lookup_table));
    lut->entries = malloc(n_output * n_chromosomes * sizeof(sample_lookup_entry));
    lut->n_samples = n_output;
    lut->n_chromosomes = n_chromosomes;

    /* Initialize all to -1 */
    for (int i = 0; i < n_output * n_chromosomes; i++) {
        lut->entries[i].pat_idx = -1;
        lut->entries[i].mat_idx = -1;
    }

    /* Build lookup by scanning chromosomes once */
    for (int c = 0; c < pts->n_chroms; c++) {
        sample_chrom* chr = &pts->chroms[c];
        int chrom_idx = chr->chrom_id - 1;  /* 0-based index */

        if (chrom_idx < 0 || chrom_idx >= n_chromosomes) continue;

        /* Find which output sample this belongs to */
        for (int s = 0; s < n_output; s++) {
            if (strcmp(output_samples[s], chr->name) == 0) {
                int lut_idx = s * n_chromosomes + chrom_idx;
                if (chr->homolog == 0) {
                    lut->entries[lut_idx].pat_idx = c;
                } else {
                    lut->entries[lut_idx].mat_idx = c;
                }
                break;
            }
        }
    }

    return lut;
}

static void free_sample_lookup(sample_lookup_table* lut)
{
    if (!lut) return;
    free(lut->entries);
    free(lut);
}

/* Binary search to find segment containing position (as fraction 0-1) */
static inline int find_segment_binary(segment* segments, int n_segments, double pos_frac)
{
    int lo = 0, hi = n_segments - 1;

    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        segment* seg = &segments[mid];

        if (pos_frac < seg->start) {
            hi = mid - 1;
        } else if (pos_frac >= seg->end) {
            lo = mid + 1;
        } else {
            /* pos_frac is within [start, end) */
            return seg->homolog;
        }
    }

    return -1;  /* Not found (gap in segments) */
}

/* Get VCF haplotype index using precomputed lookup table (multi-chromosome) */
static inline int get_haplotype_fast(pedtrans_segments* pts,
                                      sample_lookup_table* lut,
                                      int sample_idx, int chrom_id, int homolog,
                                      double pos_frac)
{
    int chrom_offset = chrom_id - 1;  /* Convert to 0-based */
    if (chrom_offset < 0 || chrom_offset >= lut->n_chromosomes) return -1;

    int lut_idx = sample_idx * lut->n_chromosomes + chrom_offset;
    int pts_chrom_idx;

    if (homolog == 0) {
        pts_chrom_idx = lut->entries[lut_idx].pat_idx;
    } else {
        pts_chrom_idx = lut->entries[lut_idx].mat_idx;
    }

    if (pts_chrom_idx < 0) return -1;

    sample_chrom* chr = &pts->chroms[pts_chrom_idx];
    return find_segment_binary(chr->segments, chr->n_segments, pos_frac);
}

/* Optimized assembly function (supports multiple chromosomes) */
int assemble_sample_vcf(founder_vcf* fvcf, pedtrans_segments* pts,
                        FILE* out, int samples_only)
{
    /* Determine which samples to output */
    char** output_samples;
    int n_output;

    if (samples_only) {
        output_samples = malloc(pts->n_samples * sizeof(char*));
        n_output = 0;

        founder_map* fm = pts->founder_map;

        for (int i = 0; i < pts->n_samples; i++) {
            int is_founder = 0;
            if (fm) {
                for (int j = 0; j < fm->n_founders; j++) {
                    if (strcmp(pts->sample_names[i], fm->names[j]) == 0) {
                        is_founder = 1;
                        break;
                    }
                }
            }
            if (!is_founder) {
                output_samples[n_output++] = pts->sample_names[i];
            }
        }
    } else {
        output_samples = pts->sample_names;
        n_output = pts->n_samples;
    }

    /* Build fast lookup table - O(n_chroms) once */
    sample_lookup_table* lut = build_sample_lookup(pts, output_samples, n_output);

    /* Pre-allocate genotype array - avoids malloc/free per mutation */
    int* geno = malloc(n_output * 2 * sizeof(int));

    /* Write VCF header */
    fprintf(out, "##fileformat=VCFv4.2\n");
    fprintf(out, "##source=vcfassemble\n");

    /* Write contig lines for all chromosomes */
    for (int c = 1; c <= fvcf->n_chromosomes; c++) {
        long length = fvcf->chrom_lengths[c];
        if (length == 0) length = 10000000;  /* Default length */
        fprintf(out, "##contig=<ID=chr%d,length=%ld>\n", c, length);
    }

    fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for (int i = 0; i < n_output; i++) {
        fprintf(out, "\t%s", output_samples[i]);
    }
    fprintf(out, "\n");

    int n_founders = fvcf->n_founder_samples;

    /* Process each mutation */
    int variants_written = 0;
    for (int m = 0; m < fvcf->n_mutations; m++) {
        vcf_mutation* mut = &fvcf->mutations[m];
        int chrom_id = mut->chrom_id;

        /* Get chromosome length for position fraction calculation */
        long chrom_length = fvcf->chrom_lengths[chrom_id];
        if (chrom_length == 0) chrom_length = 10000000;
        double inv_chrom_length = 1.0 / (double)chrom_length;

        /* Convert position to fraction */
        double pos_frac = (double)mut->position * inv_chrom_length;

        /* Collect genotypes using fast lookup */
        int any_variant = 0;

        for (int s = 0; s < n_output; s++) {
            /* Paternal haplotype */
            int vcf_idx = get_haplotype_fast(pts, lut, s, chrom_id, 0, pos_frac);
            if (vcf_idx >= 0 && vcf_idx < n_founders) {
                geno[s * 2] = mut->founder_geno[vcf_idx];
                if (geno[s * 2]) any_variant = 1;
            } else {
                geno[s * 2] = 0;
            }

            /* Maternal haplotype */
            vcf_idx = get_haplotype_fast(pts, lut, s, chrom_id, 1, pos_frac);
            if (vcf_idx >= 0 && vcf_idx < n_founders) {
                geno[s * 2 + 1] = mut->founder_geno[vcf_idx];
                if (geno[s * 2 + 1]) any_variant = 1;
            } else {
                geno[s * 2 + 1] = 0;
            }
        }

        /* Only output if at least one sample has the variant */
        if (any_variant) {
            fprintf(out, "chr%d\t%ld\t.\t%c\t%c\t.\tPASS\t.\tGT",
                    chrom_id, mut->position, mut->ref, mut->alt);

            for (int s = 0; s < n_output; s++) {
                fprintf(out, "\t%d|%d", geno[s * 2], geno[s * 2 + 1]);
            }
            fprintf(out, "\n");
            variants_written++;
        }
    }

    /* Cleanup */
    free(geno);
    free_sample_lookup(lut);

    if (samples_only) {
        free(output_samples);
    }

    return 0;
}
