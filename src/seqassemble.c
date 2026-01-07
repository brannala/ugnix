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
 * FASTA Index Functions (Streaming Mode)
 * ============================================================================
 */

fasta_index* create_fasta_index(const char* filename,
                                 const char** founder_names, int n_founders) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open FASTA file '%s'\n", filename);
        return NULL;
    }

    fasta_index* idx = malloc(sizeof(fasta_index));
    if (!idx) {
        fclose(fp);
        return NULL;
    }

    idx->capacity = n_founders * 2 + 16;
    idx->entries = malloc(idx->capacity * sizeof(fasta_index_entry));
    if (!idx->entries) {
        free(idx);
        fclose(fp);
        return NULL;
    }

    idx->n_entries = 0;
    idx->fp = fp;
    idx->filename = strdup(filename);
    idx->seq_length = 0;

    /* Scan file to build index */
    char line[MAX_LINE_LENGTH];
    fasta_index_entry* current = NULL;
    long seq_bases = 0;
    int line_width = 0;
    int sample_counter = 0;

    while (fgets(line, sizeof(line), fp)) {
        int len = strlen(line);

        if (line[0] == '>') {
            /* Finalize previous entry */
            if (current) {
                current->seq_length = seq_bases;
                current->line_width = line_width;
                if (idx->seq_length == 0) {
                    idx->seq_length = seq_bases;
                }
            }

            /* Parse header */
            /* Remove trailing whitespace */
            while (len > 0 && isspace(line[len-1])) {
                line[--len] = '\0';
            }

            char* header = line + 1;  /* skip '>' */
            char founder[MAX_NAME_LENGTH];
            int homolog = 0;

            /* Try "FounderName:pat" or "FounderName:mat" format */
            char* colon = strchr(header, ':');
            if (colon) {
                *colon = '\0';
                strncpy(founder, header, MAX_NAME_LENGTH - 1);
                founder[MAX_NAME_LENGTH - 1] = '\0';
                if (strcmp(colon + 1, "mat") == 0) {
                    homolog = 1;
                }
            } else {
                /* Try "sample0", "sample1" format */
                int sample_idx = -1;
                if (sscanf(header, "sample%d", &sample_idx) == 1) {
                    int founder_idx = sample_idx / 2;
                    homolog = sample_idx % 2;
                    if (founder_idx < n_founders && founder_names) {
                        strncpy(founder, founder_names[founder_idx], MAX_NAME_LENGTH - 1);
                        founder[MAX_NAME_LENGTH - 1] = '\0';
                    } else {
                        snprintf(founder, MAX_NAME_LENGTH, "founder%d", founder_idx);
                    }
                } else {
                    /* Use as-is, assign to founders in order */
                    int founder_idx = sample_counter / 2;
                    homolog = sample_counter % 2;
                    if (founder_idx < n_founders && founder_names) {
                        strncpy(founder, founder_names[founder_idx], MAX_NAME_LENGTH - 1);
                        founder[MAX_NAME_LENGTH - 1] = '\0';
                    } else {
                        strncpy(founder, header, MAX_NAME_LENGTH - 1);
                        founder[MAX_NAME_LENGTH - 1] = '\0';
                    }
                    sample_counter++;
                }
            }

            /* Add new entry */
            if (idx->n_entries >= idx->capacity) {
                idx->capacity *= 2;
                idx->entries = realloc(idx->entries,
                                        idx->capacity * sizeof(fasta_index_entry));
            }

            current = &idx->entries[idx->n_entries++];
            strncpy(current->founder, founder, MAX_NAME_LENGTH - 1);
            current->founder[MAX_NAME_LENGTH - 1] = '\0';
            current->homolog = homolog;
            current->file_offset = ftell(fp);  /* position after header */
            current->seq_length = 0;
            current->line_width = 0;

            seq_bases = 0;
            line_width = 0;

        } else if (line[0] != '\0' && line[0] != '#') {
            /* Sequence line */
            /* Remove newline for counting */
            while (len > 0 && isspace(line[len-1])) {
                len--;
            }
            seq_bases += len;
            if (line_width == 0 && len > 0) {
                line_width = len;  /* assume consistent line width */
            }
        }
    }

    /* Finalize last entry */
    if (current) {
        current->seq_length = seq_bases;
        current->line_width = line_width;
        if (idx->seq_length == 0) {
            idx->seq_length = seq_bases;
        }
    }

    return idx;
}

void free_fasta_index(fasta_index* idx) {
    if (!idx) return;
    if (idx->fp) fclose(idx->fp);
    free(idx->filename);
    free(idx->entries);
    free(idx);
}

static fasta_index_entry* find_index_entry(fasta_index* idx,
                                            const char* founder, int homolog) {
    for (int i = 0; i < idx->n_entries; i++) {
        if (strcmp(idx->entries[i].founder, founder) == 0 &&
            idx->entries[i].homolog == homolog) {
            return &idx->entries[i];
        }
    }
    return NULL;
}

int read_sequence_range(fasta_index* idx, const char* founder,
                        int homolog, long start, long end, char* buffer) {
    fasta_index_entry* entry = find_index_entry(idx, founder, homolog);
    if (!entry) {
        fprintf(stderr, "Warning: No index entry for %s:%s\n",
                founder, homolog == 0 ? "pat" : "mat");
        return -1;
    }

    if (start < 0) start = 0;
    if (end > entry->seq_length) end = entry->seq_length;
    if (start >= end) return 0;

    long range_len = end - start;

    /* Calculate file position accounting for newlines */
    int line_width = entry->line_width > 0 ? entry->line_width : 80;
    long start_line = start / line_width;
    long start_col = start % line_width;

    /* Seek to approximate position (start of the line containing start) */
    /* Each line has line_width chars + 1 newline */
    long file_pos = entry->file_offset + start_line * (line_width + 1) + start_col;
    fseek(idx->fp, file_pos, SEEK_SET);

    /* Read bases, skipping newlines */
    long bases_read = 0;
    int c;
    while (bases_read < range_len && (c = fgetc(idx->fp)) != EOF) {
        if (c == '\n' || c == '\r') continue;
        if (c == '>') break;  /* hit next sequence */
        buffer[bases_read++] = (char)c;
    }

    return (int)bases_read;
}

/*
 * ============================================================================
 * Founder Sequence Functions (In-Memory Mode)
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
 * Streaming assembly - writes one chromosome at a time, reading
 * founder sequences on demand from indexed FASTA file.
 */
static void write_chromosome_streaming(parsed_chromosome* chr,
                                        fasta_index* idx,
                                        const char* indiv_name,
                                        const char* homolog_name,
                                        FILE* out) {
    long seq_length = idx->seq_length;

    /* Write FASTA header */
    fprintf(out, ">%s:%s\n", indiv_name, homolog_name);

    /* Buffer for reading segments - sized for one output line */
    #define STREAM_LINE_WIDTH 80
    #define STREAM_BUFFER_SIZE 4096

    char buffer[STREAM_BUFFER_SIZE];
    char line_buffer[STREAM_LINE_WIDTH + 1];
    int line_pos = 0;

    /* Process sequence position by position through segments */
    long pos = 0;
    int seg_idx = 0;

    while (pos < seq_length) {
        /* Find segment containing this position */
        while (seg_idx < chr->n_segments - 1 &&
               chr->segments[seg_idx].end * seq_length <= pos) {
            seg_idx++;
        }

        parsed_segment* seg = NULL;
        if (seg_idx < chr->n_segments) {
            seg = &chr->segments[seg_idx];
        }

        /* Calculate how many bases we can read from this segment */
        long seg_start_bp = seg ? (long)(seg->start * seq_length) : seq_length;
        long seg_end_bp = seg ? (long)(seg->end * seq_length) : seq_length;

        /* Handle gap before segment (fill with N) */
        while (pos < seg_start_bp && pos < seq_length) {
            line_buffer[line_pos++] = 'N';
            if (line_pos >= STREAM_LINE_WIDTH) {
                line_buffer[line_pos] = '\0';
                fprintf(out, "%s\n", line_buffer);
                line_pos = 0;
            }
            pos++;
        }

        /* Read from segment */
        if (seg && pos < seg_end_bp) {
            long chunk_start = pos;
            long chunk_end = seg_end_bp;
            if (chunk_end > seq_length) chunk_end = seq_length;

            /* Read in chunks */
            while (chunk_start < chunk_end) {
                long read_len = chunk_end - chunk_start;
                if (read_len > STREAM_BUFFER_SIZE) {
                    read_len = STREAM_BUFFER_SIZE;
                }

                int bytes_read = read_sequence_range(idx, seg->founder,
                                                      seg->homolog,
                                                      chunk_start, chunk_start + read_len,
                                                      buffer);
                if (bytes_read <= 0) {
                    /* Fill with N on error */
                    for (long j = 0; j < read_len; j++) {
                        line_buffer[line_pos++] = 'N';
                        if (line_pos >= STREAM_LINE_WIDTH) {
                            line_buffer[line_pos] = '\0';
                            fprintf(out, "%s\n", line_buffer);
                            line_pos = 0;
                        }
                    }
                } else {
                    /* Write read bases */
                    for (int j = 0; j < bytes_read; j++) {
                        line_buffer[line_pos++] = buffer[j];
                        if (line_pos >= STREAM_LINE_WIDTH) {
                            line_buffer[line_pos] = '\0';
                            fprintf(out, "%s\n", line_buffer);
                            line_pos = 0;
                        }
                    }
                }

                chunk_start += read_len;
                pos = chunk_start;
            }
        }
    }

    /* Write any remaining bases in line buffer */
    if (line_pos > 0) {
        line_buffer[line_pos] = '\0';
        fprintf(out, "%s\n", line_buffer);
    }
}

int assemble_all_sequences_streaming(pedtrans_data* pd, fasta_index* idx,
                                      FILE* out, int samples_only) {
    for (int i = 0; i < pd->n_individuals; i++) {
        parsed_individual* indiv = &pd->individuals[i];

        /* Skip founders if samples_only is set */
        if (samples_only && indiv->is_founder) {
            continue;
        }

        /* Assemble and write paternal chromosome */
        write_chromosome_streaming(&indiv->paternal, idx,
                                    indiv->name, "pat", out);

        /* Assemble and write maternal chromosome */
        write_chromosome_streaming(&indiv->maternal, idx,
                                    indiv->name, "mat", out);
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
