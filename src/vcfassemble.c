#include <ctype.h>
#include "vcfassemble.h"

/* Read founder VCF file */
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
    fvcf->chrom_length = 0;
    fvcf->founder_names = NULL;
    fvcf->n_founder_samples = 0;

    char line[MAX_LINE_LENGTH];

    /* Parse header */
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "##contig", 8) == 0) {
            /* Extract chromosome length */
            char* len_ptr = strstr(line, "length=");
            if (len_ptr) {
                fvcf->chrom_length = atol(len_ptr + 7);
            }
        }
        else if (line[0] == '#' && line[1] != '#') {
            /* Header line with sample names */
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

        char* tok = strtok(line, "\t");
        int col = 0;
        int sample_idx = 0;

        while (tok) {
            switch (col) {
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

/* Free founder VCF */
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

    free(fvcf);
}

/* Founder name to VCF index mapping */
typedef struct {
    char** names;
    int n_founders;
} founder_map;

/* Parse founder list from pedtrans header: "# Founders: I35 I37 ..." */
static founder_map* parse_founder_list(FILE* fp)
{
    char line[MAX_LINE_LENGTH];
    founder_map* fm = malloc(sizeof(founder_map));
    fm->names = malloc(256 * sizeof(char*));
    fm->n_founders = 0;

    rewind(fp);
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "# Founders:", 11) == 0) {
            char* tok = strtok(line + 11, " \t\n");
            while (tok && fm->n_founders < 256) {
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

/* Get VCF haplotype index from founder:homolog (e.g., "I35:pat") */
static int get_vcf_haplotype_idx(founder_map* fm, const char* founder_hom)
{
    /* Parse "I35:pat" or "I35:mat" */
    char founder[MAX_NAME_LENGTH];
    char homolog[16];

    if (sscanf(founder_hom, "%[^:]:%s", founder, homolog) != 2) {
        return -1;
    }

    /* Find founder index */
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

    /* pat=0, mat=1 */
    int hom = (strcmp(homolog, "mat") == 0) ? 1 : 0;

    return founder_idx * 2 + hom;
}

/* Read pedtrans segment output */
pedtrans_segments* read_pedtrans_segments(const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open pedtrans file %s\n", filename);
        return NULL;
    }

    /* Parse founder list first */
    founder_map* fm = parse_founder_list(fp);

    pedtrans_segments* pts = malloc(sizeof(pedtrans_segments));
    pts->chroms = malloc(INITIAL_CAPACITY * sizeof(sample_chrom));
    pts->n_chroms = 0;
    pts->capacity = INITIAL_CAPACITY;
    pts->sample_names = malloc(INITIAL_CAPACITY * sizeof(char*));
    pts->n_samples = 0;

    /* Store founder map in pts for later use */
    pts->founder_map = fm;

    char line[MAX_LINE_LENGTH];
    char current_individual[MAX_NAME_LENGTH] = "";
    int current_homolog = -1;  /* 0=paternal, 1=maternal */
    sample_chrom* current = NULL;
    int sample_capacity = INITIAL_CAPACITY;

    while (fgets(line, sizeof(line), fp)) {
        /* Skip comments and empty lines */
        if (line[0] == '#' || line[0] == '\n') continue;

        /* Check for "Individual: NAME" */
        if (strncmp(line, "Individual:", 11) == 0) {
            sscanf(line + 11, " %s", current_individual);
            current_homolog = -1;

            /* Track unique sample names */
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

        /* Check for "  Paternal:" or "  Maternal:" */
        char* trimmed = line;
        while (*trimmed == ' ') trimmed++;

        if (strncmp(trimmed, "Paternal:", 9) == 0) {
            current_homolog = 0;

            /* Create new sample_chrom entry */
            if (pts->n_chroms >= pts->capacity) {
                pts->capacity *= 2;
                pts->chroms = realloc(pts->chroms,
                                      pts->capacity * sizeof(sample_chrom));
            }
            current = &pts->chroms[pts->n_chroms];
            strncpy(current->name, current_individual, MAX_NAME_LENGTH - 1);
            current->name[MAX_NAME_LENGTH - 1] = '\0';
            current->homolog = 0;
            current->segments = malloc(64 * sizeof(segment));
            current->n_segments = 0;
            current->capacity = 64;
            pts->n_chroms++;
            continue;
        }

        if (strncmp(trimmed, "Maternal:", 9) == 0) {
            current_homolog = 1;

            /* Create new sample_chrom entry */
            if (pts->n_chroms >= pts->capacity) {
                pts->capacity *= 2;
                pts->chroms = realloc(pts->chroms,
                                      pts->capacity * sizeof(sample_chrom));
            }
            current = &pts->chroms[pts->n_chroms];
            strncpy(current->name, current_individual, MAX_NAME_LENGTH - 1);
            current->name[MAX_NAME_LENGTH - 1] = '\0';
            current->homolog = 1;
            current->segments = malloc(64 * sizeof(segment));
            current->n_segments = 0;
            current->capacity = 64;
            pts->n_chroms++;
            continue;
        }

        /* Parse segment line: "    0.000000 1.000000 I35:pat" */
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

                /* Store the founder:homolog string for later lookup */
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

/* Free pedtrans segments */
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

/* Check if position falls within segment */
static int position_in_segment(long pos, long chrom_length, segment* seg)
{
    double pos_frac = (double)pos / chrom_length;
    return (pos_frac >= seg->start && pos_frac < seg->end);
}

/* Get VCF haplotype index for a sample at a given position */
static int get_sample_haplotype_at_pos(pedtrans_segments* pts, long chrom_length,
                                       const char* sample, int homolog, long pos)
{
    for (int i = 0; i < pts->n_chroms; i++) {
        if (strcmp(pts->chroms[i].name, sample) == 0 &&
            pts->chroms[i].homolog == homolog) {

            for (int j = 0; j < pts->chroms[i].n_segments; j++) {
                segment* seg = &pts->chroms[i].segments[j];
                if (position_in_segment(pos, chrom_length, seg)) {
                    return seg->homolog;  /* Already contains VCF index */
                }
            }
        }
    }
    return -1;
}

/* Get founder haplotype index (not used with new format, kept for API) */
int get_founder_haplotype_idx(founder_vcf* fvcf, const char* founder, int homolog)
{
    (void)fvcf;
    (void)founder;
    (void)homolog;
    return -1;  /* Use segment->homolog directly instead */
}

/* Assemble sample VCF from founder VCF and pedtrans segments */
int assemble_sample_vcf(founder_vcf* fvcf, pedtrans_segments* pts,
                        FILE* out, int samples_only)
{
    /* Determine which samples to output */
    char** output_samples;
    int n_output;

    if (samples_only) {
        /* Only output non-founder samples */
        output_samples = malloc(pts->n_samples * sizeof(char*));
        n_output = 0;

        /* Get founder set for filtering */
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

    /* Write VCF header */
    fprintf(out, "##fileformat=VCFv4.2\n");
    fprintf(out, "##source=vcfassemble\n");
    fprintf(out, "##contig=<ID=chr1,length=%ld>\n", fvcf->chrom_length);
    fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for (int i = 0; i < n_output; i++) {
        fprintf(out, "\t%s", output_samples[i]);
    }
    fprintf(out, "\n");

    /* Process each mutation */
    for (int m = 0; m < fvcf->n_mutations; m++) {
        vcf_mutation* mut = &fvcf->mutations[m];

        /* Collect genotypes for all output samples */
        int* geno = malloc(n_output * 2 * sizeof(int));
        int any_variant = 0;

        for (int s = 0; s < n_output; s++) {
            for (int h = 0; h < 2; h++) {
                int vcf_idx = get_sample_haplotype_at_pos(
                    pts, fvcf->chrom_length, output_samples[s], h, mut->position);

                if (vcf_idx >= 0 && vcf_idx < mut->n_founders) {
                    geno[s * 2 + h] = mut->founder_geno[vcf_idx];
                    if (geno[s * 2 + h]) any_variant = 1;
                } else {
                    geno[s * 2 + h] = 0;
                }
            }
        }

        /* Only output if at least one sample has the variant */
        if (any_variant) {
            fprintf(out, "chr1\t%ld\t.\t%c\t%c\t.\tPASS\t.\tGT",
                    mut->position, mut->ref, mut->alt);

            for (int s = 0; s < n_output; s++) {
                fprintf(out, "\t%d|%d", geno[s * 2], geno[s * 2 + 1]);
            }
            fprintf(out, "\n");
        }

        free(geno);
    }

    if (samples_only) {
        free(output_samples);
    }

    return 0;
}
