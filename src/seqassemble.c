/*
 * seqassemble.c - Sequence Assembly from Pedigree Segments
 *
 * Combines founder sequences (from coalsim) with segment transmission data
 * (from pedtrans) to assemble final sample sequences.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include "seqassemble.h"
#include <ctype.h>

/*
 * ============================================================================
 * Founder Sequence Functions
 * ============================================================================
 */

founder_sequences* create_founder_sequences(void) {
    founder_sequences* fs = malloc(sizeof(founder_sequences));
    if (!fs) return NULL;

    fs->capacity = 16;
    fs->founders = malloc(fs->capacity * sizeof(founder_haplotypes));
    if (!fs->founders) {
        free(fs);
        return NULL;
    }

    fs->n_founders = 0;
    fs->seq_length = 0;
    return fs;
}

void free_founder_sequences(founder_sequences* fs) {
    if (!fs) return;

    for (int i = 0; i < fs->n_founders; i++) {
        free(fs->founders[i].name);
        free(fs->founders[i].pat_seq);
        free(fs->founders[i].mat_seq);
    }
    free(fs->founders);
    free(fs);
}

static int add_founder(founder_sequences* fs, const char* name) {
    if (fs->n_founders >= fs->capacity) {
        int new_cap = fs->capacity * 2;
        founder_haplotypes* new_arr = realloc(fs->founders,
                                               new_cap * sizeof(founder_haplotypes));
        if (!new_arr) return -1;
        fs->founders = new_arr;
        fs->capacity = new_cap;
    }

    int idx = fs->n_founders;
    fs->founders[idx].name = strdup(name);
    fs->founders[idx].pat_seq = NULL;
    fs->founders[idx].mat_seq = NULL;
    fs->founders[idx].seq_length = 0;
    fs->n_founders++;
    return idx;
}

static int find_founder(founder_sequences* fs, const char* name) {
    for (int i = 0; i < fs->n_founders; i++) {
        if (strcmp(fs->founders[i].name, name) == 0) {
            return i;
        }
    }
    return -1;
}

int read_founder_fasta(founder_sequences* fs, const char* filename,
                       const char** founder_names, int n_founders) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open FASTA file '%s'\n", filename);
        return -1;
    }

    /* Pre-populate founders from the provided names */
    for (int i = 0; i < n_founders; i++) {
        add_founder(fs, founder_names[i]);
    }

    char line[MAX_LINE_LENGTH];
    int current_homolog = -1;  /* 0=pat, 1=mat */
    int current_founder_idx = -1;
    char* seq_buffer = NULL;
    long seq_len = 0;
    long seq_cap = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* Remove trailing whitespace */
        int len = strlen(line);
        while (len > 0 && isspace(line[len-1])) {
            line[--len] = '\0';
        }

        if (line[0] == '>') {
            /* Save previous sequence if any */
            if (current_founder_idx >= 0 && seq_buffer) {
                if (current_homolog == 0) {
                    fs->founders[current_founder_idx].pat_seq = seq_buffer;
                } else {
                    fs->founders[current_founder_idx].mat_seq = seq_buffer;
                }
                fs->founders[current_founder_idx].seq_length = seq_len;
                if (fs->seq_length == 0) {
                    fs->seq_length = seq_len;
                }
                seq_buffer = NULL;
            }

            /* Parse header line */
            char* header = line + 1;  /* skip '>' */

            /* Try to parse "FounderName:pat" or "FounderName:mat" format */
            char* colon = strchr(header, ':');
            if (colon) {
                *colon = '\0';
                char* homolog_str = colon + 1;

                int fidx = find_founder(fs, header);
                if (fidx < 0) {
                    fidx = add_founder(fs, header);
                }
                current_founder_idx = fidx;

                if (strcmp(homolog_str, "pat") == 0) {
                    current_homolog = 0;
                } else if (strcmp(homolog_str, "mat") == 0) {
                    current_homolog = 1;
                } else {
                    fprintf(stderr, "Warning: Unknown homolog '%s', assuming paternal\n",
                            homolog_str);
                    current_homolog = 0;
                }
            } else {
                /* Try "sample0", "sample1" format */
                /* Even indices = paternal, odd = maternal */
                int sample_idx = -1;
                if (sscanf(header, "sample%d", &sample_idx) == 1) {
                    int founder_idx = sample_idx / 2;
                    current_homolog = sample_idx % 2;

                    if (founder_idx < n_founders) {
                        current_founder_idx = founder_idx;
                    } else {
                        fprintf(stderr, "Warning: sample%d exceeds founder count\n",
                                sample_idx);
                        current_founder_idx = -1;
                    }
                } else {
                    /* Unknown format - try to use as founder name with paternal */
                    int fidx = find_founder(fs, header);
                    if (fidx < 0) {
                        fidx = add_founder(fs, header);
                    }
                    current_founder_idx = fidx;
                    current_homolog = 0;
                }
            }

            /* Start new sequence buffer */
            seq_cap = 1024;
            seq_buffer = malloc(seq_cap);
            seq_len = 0;
            if (seq_buffer) seq_buffer[0] = '\0';

        } else if (line[0] != '\0' && line[0] != '#') {
            /* Sequence line - append to buffer */
            if (seq_buffer && current_founder_idx >= 0) {
                int add_len = strlen(line);
                if (seq_len + add_len + 1 > seq_cap) {
                    seq_cap = (seq_len + add_len + 1) * 2;
                    char* new_buf = realloc(seq_buffer, seq_cap);
                    if (!new_buf) {
                        free(seq_buffer);
                        seq_buffer = NULL;
                        continue;
                    }
                    seq_buffer = new_buf;
                }
                memcpy(seq_buffer + seq_len, line, add_len);
                seq_len += add_len;
                seq_buffer[seq_len] = '\0';
            }
        }
    }

    /* Save last sequence */
    if (current_founder_idx >= 0 && seq_buffer) {
        if (current_homolog == 0) {
            fs->founders[current_founder_idx].pat_seq = seq_buffer;
        } else {
            fs->founders[current_founder_idx].mat_seq = seq_buffer;
        }
        fs->founders[current_founder_idx].seq_length = seq_len;
        if (fs->seq_length == 0) {
            fs->seq_length = seq_len;
        }
    }

    fclose(fp);
    return 0;
}

const char* get_founder_sequence(founder_sequences* fs,
                                  const char* founder_name, int homolog) {
    int idx = find_founder(fs, founder_name);
    if (idx < 0) return NULL;

    if (homolog == 0) {
        return fs->founders[idx].pat_seq;
    } else {
        return fs->founders[idx].mat_seq;
    }
}

/*
 * ============================================================================
 * Pedtrans Data Functions
 * ============================================================================
 */

pedtrans_data* create_pedtrans_data(void) {
    pedtrans_data* pd = malloc(sizeof(pedtrans_data));
    if (!pd) return NULL;

    pd->capacity = 64;
    pd->individuals = malloc(pd->capacity * sizeof(parsed_individual));
    if (!pd->individuals) {
        free(pd);
        return NULL;
    }

    pd->n_individuals = 0;
    pd->founder_names = NULL;
    pd->n_founders = 0;
    return pd;
}

void free_pedtrans_data(pedtrans_data* pd) {
    if (!pd) return;

    for (int i = 0; i < pd->n_individuals; i++) {
        free(pd->individuals[i].paternal.segments);
        free(pd->individuals[i].maternal.segments);
    }
    free(pd->individuals);

    if (pd->founder_names) {
        for (int i = 0; i < pd->n_founders; i++) {
            free(pd->founder_names[i]);
        }
        free(pd->founder_names);
    }

    free(pd);
}

static void init_chromosome(parsed_chromosome* chr) {
    chr->capacity = 8;
    chr->segments = malloc(chr->capacity * sizeof(parsed_segment));
    chr->n_segments = 0;
}

static void add_segment(parsed_chromosome* chr, double start, double end,
                        const char* founder, int homolog) {
    if (chr->n_segments >= chr->capacity) {
        chr->capacity *= 2;
        chr->segments = realloc(chr->segments,
                                 chr->capacity * sizeof(parsed_segment));
    }

    parsed_segment* seg = &chr->segments[chr->n_segments++];
    seg->start = start;
    seg->end = end;
    strncpy(seg->founder, founder, MAX_NAME_LENGTH - 1);
    seg->founder[MAX_NAME_LENGTH - 1] = '\0';
    seg->homolog = homolog;
}

int parse_pedtrans_output(pedtrans_data* pd, const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open pedtrans file '%s'\n", filename);
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    parsed_individual* current = NULL;
    int parsing_paternal = 0;
    int parsing_maternal = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* Remove trailing whitespace */
        int len = strlen(line);
        while (len > 0 && isspace(line[len-1])) {
            line[--len] = '\0';
        }

        /* Skip empty lines */
        if (len == 0) continue;

        /* Parse header comments */
        if (line[0] == '#') {
            /* Look for "# Founders:" line */
            if (strncmp(line, "# Founders:", 11) == 0) {
                char* ptr = line + 11;
                /* Count and parse founder names */
                int cap = 16;
                pd->founder_names = malloc(cap * sizeof(char*));
                pd->n_founders = 0;

                char* token = strtok(ptr, " \t");
                while (token) {
                    if (pd->n_founders >= cap) {
                        cap *= 2;
                        pd->founder_names = realloc(pd->founder_names,
                                                     cap * sizeof(char*));
                    }
                    pd->founder_names[pd->n_founders++] = strdup(token);
                    token = strtok(NULL, " \t");
                }
            }
            continue;
        }

        /* Parse "Individual: NAME" */
        if (strncmp(line, "Individual:", 11) == 0) {
            /* Add new individual */
            if (pd->n_individuals >= pd->capacity) {
                pd->capacity *= 2;
                pd->individuals = realloc(pd->individuals,
                                           pd->capacity * sizeof(parsed_individual));
            }

            current = &pd->individuals[pd->n_individuals++];
            char* name = line + 11;
            while (*name == ' ') name++;
            strncpy(current->name, name, MAX_NAME_LENGTH - 1);
            current->name[MAX_NAME_LENGTH - 1] = '\0';

            init_chromosome(&current->paternal);
            init_chromosome(&current->maternal);

            /* Check if this is a founder */
            current->is_founder = 0;
            for (int i = 0; i < pd->n_founders; i++) {
                if (strcmp(current->name, pd->founder_names[i]) == 0) {
                    current->is_founder = 1;
                    break;
                }
            }

            parsing_paternal = 0;
            parsing_maternal = 0;
            continue;
        }

        /* Parse "  Paternal:" or "  Maternal:" */
        char* trimmed = line;
        while (*trimmed == ' ') trimmed++;

        if (strncmp(trimmed, "Paternal:", 9) == 0) {
            parsing_paternal = 1;
            parsing_maternal = 0;
            continue;
        }

        if (strncmp(trimmed, "Maternal:", 9) == 0) {
            parsing_paternal = 0;
            parsing_maternal = 1;
            continue;
        }

        /* Parse segment line: "    0.000000 0.500000 FounderName:pat" */
        if (current && (parsing_paternal || parsing_maternal)) {
            double start, end;
            char origin[MAX_NAME_LENGTH];

            if (sscanf(trimmed, "%lf %lf %s", &start, &end, origin) == 3) {
                /* Parse origin as "FounderName:pat" or "FounderName:mat" */
                char* colon = strchr(origin, ':');
                int homolog = 0;
                if (colon) {
                    *colon = '\0';
                    if (strcmp(colon + 1, "mat") == 0) {
                        homolog = 1;
                    }
                }

                if (parsing_paternal) {
                    add_segment(&current->paternal, start, end, origin, homolog);
                } else {
                    add_segment(&current->maternal, start, end, origin, homolog);
                }
            }
        }
    }

    fclose(fp);
    return 0;
}

/*
 * ============================================================================
 * Sequence Assembly Functions
 * ============================================================================
 */

char* assemble_chromosome_sequence(parsed_chromosome* chr,
                                    founder_sequences* fs,
                                    long seq_length) {
    char* sequence = malloc(seq_length + 1);
    if (!sequence) return NULL;

    /* Initialize with 'N' (unknown) */
    memset(sequence, 'N', seq_length);
    sequence[seq_length] = '\0';

    for (int i = 0; i < chr->n_segments; i++) {
        parsed_segment* seg = &chr->segments[i];

        /* Get founder sequence */
        const char* founder_seq = get_founder_sequence(fs, seg->founder,
                                                        seg->homolog);
        if (!founder_seq) {
            fprintf(stderr, "Warning: No sequence for %s:%s\n",
                    seg->founder, seg->homolog == 0 ? "pat" : "mat");
            continue;
        }

        /* Calculate base positions */
        long start_bp = (long)(seg->start * seq_length);
        long end_bp = (long)(seg->end * seq_length);

        /* Clip to valid range */
        if (start_bp < 0) start_bp = 0;
        if (end_bp > seq_length) end_bp = seq_length;

        /* Copy sequence */
        for (long pos = start_bp; pos < end_bp; pos++) {
            sequence[pos] = founder_seq[pos];
        }
    }

    return sequence;
}

int assemble_all_sequences(pedtrans_data* pd, founder_sequences* fs,
                           FILE* out, int samples_only) {
    long seq_length = fs->seq_length;

    for (int i = 0; i < pd->n_individuals; i++) {
        parsed_individual* indiv = &pd->individuals[i];

        /* Skip founders if samples_only is set */
        if (samples_only && indiv->is_founder) {
            continue;
        }

        /* Assemble paternal chromosome */
        char* pat_seq = assemble_chromosome_sequence(&indiv->paternal,
                                                      fs, seq_length);
        if (pat_seq) {
            fprintf(out, ">%s:pat\n", indiv->name);
            /* Write in lines of 80 characters */
            for (long pos = 0; pos < seq_length; pos += 80) {
                long line_len = (seq_length - pos < 80) ? (seq_length - pos) : 80;
                fprintf(out, "%.*s\n", (int)line_len, pat_seq + pos);
            }
            free(pat_seq);
        }

        /* Assemble maternal chromosome */
        char* mat_seq = assemble_chromosome_sequence(&indiv->maternal,
                                                      fs, seq_length);
        if (mat_seq) {
            fprintf(out, ">%s:mat\n", indiv->name);
            for (long pos = 0; pos < seq_length; pos += 80) {
                long line_len = (seq_length - pos < 80) ? (seq_length - pos) : 80;
                fprintf(out, "%.*s\n", (int)line_len, mat_seq + pos);
            }
            free(mat_seq);
        }
    }

    return 0;
}

/*
 * ============================================================================
 * VCF Output Functions
 * ============================================================================
 */

int write_vcf_output(pedtrans_data* pd, founder_sequences* fs,
                     FILE* out, int samples_only) {
    long seq_length = fs->seq_length;

    /* Assemble all sequences first */
    int n_samples = 0;
    for (int i = 0; i < pd->n_individuals; i++) {
        if (!samples_only || !pd->individuals[i].is_founder) {
            n_samples++;
        }
    }

    char** pat_seqs = malloc(pd->n_individuals * sizeof(char*));
    char** mat_seqs = malloc(pd->n_individuals * sizeof(char*));

    for (int i = 0; i < pd->n_individuals; i++) {
        parsed_individual* indiv = &pd->individuals[i];
        pat_seqs[i] = assemble_chromosome_sequence(&indiv->paternal, fs, seq_length);
        mat_seqs[i] = assemble_chromosome_sequence(&indiv->maternal, fs, seq_length);
    }

    /* Write VCF header */
    fprintf(out, "##fileformat=VCFv4.2\n");
    fprintf(out, "##source=seqassemble\n");
    fprintf(out, "##contig=<ID=chr1,length=%ld>\n", seq_length);
    fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for (int i = 0; i < pd->n_individuals; i++) {
        if (!samples_only || !pd->individuals[i].is_founder) {
            fprintf(out, "\t%s", pd->individuals[i].name);
        }
    }
    fprintf(out, "\n");

    /* Find polymorphic sites and output VCF records */
    for (long pos = 0; pos < seq_length; pos++) {
        /* Check if this position is polymorphic */
        char ref = 'N';
        int has_alt = 0;
        char alt = 'N';

        /* Find reference allele (most common) */
        int counts[5] = {0}; /* A, C, G, T, N */
        for (int i = 0; i < pd->n_individuals; i++) {
            if (samples_only && pd->individuals[i].is_founder) continue;
            if (pat_seqs[i]) {
                char c = pat_seqs[i][pos];
                if (c == 'A') counts[0]++;
                else if (c == 'C') counts[1]++;
                else if (c == 'G') counts[2]++;
                else if (c == 'T') counts[3]++;
            }
            if (mat_seqs[i]) {
                char c = mat_seqs[i][pos];
                if (c == 'A') counts[0]++;
                else if (c == 'C') counts[1]++;
                else if (c == 'G') counts[2]++;
                else if (c == 'T') counts[3]++;
            }
        }

        /* Find most common allele as reference */
        int max_count = 0;
        int max_idx = -1;
        char bases[] = "ACGT";
        for (int j = 0; j < 4; j++) {
            if (counts[j] > max_count) {
                max_count = counts[j];
                max_idx = j;
            }
        }

        if (max_idx < 0) continue;  /* all N */
        ref = bases[max_idx];

        /* Check for alternate allele */
        for (int j = 0; j < 4; j++) {
            if (j != max_idx && counts[j] > 0) {
                alt = bases[j];
                has_alt = 1;
                break;
            }
        }

        if (!has_alt) continue;  /* not polymorphic */

        /* Output VCF record */
        fprintf(out, "chr1\t%ld\t.\t%c\t%c\t.\tPASS\t.\tGT", pos + 1, ref, alt);

        for (int i = 0; i < pd->n_individuals; i++) {
            if (samples_only && pd->individuals[i].is_founder) continue;

            char p = pat_seqs[i] ? pat_seqs[i][pos] : 'N';
            char m = mat_seqs[i] ? mat_seqs[i][pos] : 'N';

            int a1 = (p == alt) ? 1 : 0;
            int a2 = (m == alt) ? 1 : 0;

            fprintf(out, "\t%d|%d", a1, a2);
        }
        fprintf(out, "\n");
    }

    /* Cleanup */
    for (int i = 0; i < pd->n_individuals; i++) {
        free(pat_seqs[i]);
        free(mat_seqs[i]);
    }
    free(pat_seqs);
    free(mat_seqs);

    return 0;
}
