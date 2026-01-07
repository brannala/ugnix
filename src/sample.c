/*
 * sample.c - VCF Subsampling Tool Implementation
 *
 * Subsample individuals and/or markers from a VCF file.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "sample.h"

#define MAX_LINE_LEN 1048576  /* 1MB max line length */
#define INITIAL_CAPACITY 1024

void init_sample_params(sample_params_t* params) {
    params->input_file = NULL;
    params->output_file = NULL;
    params->indiv_file = NULL;
    params->n_markers = 0;
    params->method = METHOD_RANDOM;
    params->region.chrom = NULL;
    params->region.start = 0;
    params->region.end = 0;
    params->seed = 0;
    params->verbose = 0;
}

int parse_region(const char* region_str, vcf_region_t* region) {
    if (!region_str || !region) return -1;

    region->chrom = NULL;
    region->start = 0;
    region->end = 0;

    /* Find colon separator */
    const char* colon = strchr(region_str, ':');
    if (!colon) {
        /* Just chromosome name */
        region->chrom = strdup(region_str);
        return region->chrom ? 0 : -1;
    }

    /* Extract chromosome */
    size_t chrom_len = colon - region_str;
    region->chrom = malloc(chrom_len + 1);
    if (!region->chrom) return -1;
    strncpy(region->chrom, region_str, chrom_len);
    region->chrom[chrom_len] = '\0';

    /* Parse start-end */
    const char* dash = strchr(colon + 1, '-');
    if (!dash) {
        /* Just start position */
        region->start = atol(colon + 1);
        return 0;
    }

    region->start = atol(colon + 1);
    region->end = atol(dash + 1);

    return 0;
}

void free_region(vcf_region_t* region) {
    if (region && region->chrom) {
        free(region->chrom);
        region->chrom = NULL;
    }
}

indiv_list_t* load_indiv_list(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) return NULL;

    indiv_list_t* list = malloc(sizeof(indiv_list_t));
    if (!list) {
        fclose(fp);
        return NULL;
    }

    list->capacity = INITIAL_CAPACITY;
    list->n_ids = 0;
    list->ids = malloc(sizeof(char*) * list->capacity);
    if (!list->ids) {
        free(list);
        fclose(fp);
        return NULL;
    }

    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        /* Remove trailing whitespace */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' ||
                          line[len-1] == ' ' || line[len-1] == '\t')) {
            line[--len] = '\0';
        }
        if (len == 0) continue;  /* skip empty lines */

        /* Expand if needed */
        if (list->n_ids >= list->capacity) {
            list->capacity *= 2;
            char** new_ids = realloc(list->ids, sizeof(char*) * list->capacity);
            if (!new_ids) {
                for (int i = 0; i < list->n_ids; i++) free(list->ids[i]);
                free(list->ids);
                free(list);
                fclose(fp);
                return NULL;
            }
            list->ids = new_ids;
        }

        list->ids[list->n_ids] = strdup(line);
        if (!list->ids[list->n_ids]) {
            for (int i = 0; i < list->n_ids; i++) free(list->ids[i]);
            free(list->ids);
            free(list);
            fclose(fp);
            return NULL;
        }
        list->n_ids++;
    }

    fclose(fp);
    return list;
}

void free_indiv_list(indiv_list_t* list) {
    if (!list) return;
    if (list->ids) {
        for (int i = 0; i < list->n_ids; i++) {
            free(list->ids[i]);
        }
        free(list->ids);
    }
    free(list);
}

int indiv_in_list(const indiv_list_t* list, const char* id) {
    if (!list || !id) return 0;
    for (int i = 0; i < list->n_ids; i++) {
        if (strcmp(list->ids[i], id) == 0) return 1;
    }
    return 0;
}

/* Helper: create marker collection */
static marker_collection_t* create_marker_collection(void) {
    marker_collection_t* mc = malloc(sizeof(marker_collection_t));
    if (!mc) return NULL;
    mc->capacity = 64;
    mc->n_chroms = 0;
    mc->chroms = malloc(sizeof(chrom_markers_t) * mc->capacity);
    if (!mc->chroms) {
        free(mc);
        return NULL;
    }
    return mc;
}

/* Helper: find or create chromosome in collection */
static chrom_markers_t* get_or_create_chrom(marker_collection_t* mc, const char* chrom) {
    /* Search existing */
    for (int i = 0; i < mc->n_chroms; i++) {
        if (strcmp(mc->chroms[i].chrom, chrom) == 0) {
            return &mc->chroms[i];
        }
    }

    /* Create new */
    if (mc->n_chroms >= mc->capacity) {
        mc->capacity *= 2;
        chrom_markers_t* new_chroms = realloc(mc->chroms,
                                              sizeof(chrom_markers_t) * mc->capacity);
        if (!new_chroms) return NULL;
        mc->chroms = new_chroms;
    }

    chrom_markers_t* cm = &mc->chroms[mc->n_chroms];
    cm->chrom = strdup(chrom);
    cm->capacity = INITIAL_CAPACITY;
    cm->n_markers = 0;
    cm->markers = malloc(sizeof(marker_info_t) * cm->capacity);
    if (!cm->markers) {
        free(cm->chrom);
        return NULL;
    }
    mc->n_chroms++;
    return cm;
}

/* Helper: add marker to chromosome */
static int add_marker(chrom_markers_t* cm, long pos, long file_offset) {
    if (cm->n_markers >= cm->capacity) {
        cm->capacity *= 2;
        marker_info_t* new_markers = realloc(cm->markers,
                                             sizeof(marker_info_t) * cm->capacity);
        if (!new_markers) return -1;
        cm->markers = new_markers;
    }

    marker_info_t* m = &cm->markers[cm->n_markers];
    m->chrom = cm->chrom;  /* shared pointer */
    m->pos = pos;
    m->file_offset = file_offset;
    m->selected = 0;
    cm->n_markers++;
    return 0;
}

/* Helper: free marker collection */
static void free_marker_collection(marker_collection_t* mc) {
    if (!mc) return;
    for (int i = 0; i < mc->n_chroms; i++) {
        free(mc->chroms[i].chrom);
        free(mc->chroms[i].markers);
    }
    free(mc->chroms);
    free(mc);
}

/* Helper: check if position is in region */
static int pos_in_region(const vcf_region_t* region, const char* chrom, long pos) {
    if (!region->chrom) {
        /* No region specified, include all */
        return 1;
    }
    if (strcmp(region->chrom, chrom) != 0) {
        return 0;
    }
    if (region->start > 0 && pos < region->start) {
        return 0;
    }
    if (region->end > 0 && pos > region->end) {
        return 0;
    }
    return 1;
}

/* Select markers uniformly spaced */
static void select_uniform(chrom_markers_t* cm, int n_select, gsl_rng* rng) {
    if (n_select <= 0 || cm->n_markers == 0) return;
    if (n_select >= cm->n_markers) {
        /* Select all */
        for (int i = 0; i < cm->n_markers; i++) {
            cm->markers[i].selected = 1;
        }
        return;
    }

    /* Calculate spacing to maximize average distance */
    /* Select n_select markers from n_markers positions */
    /* Use floating-point step to distribute evenly */
    double step = (double)(cm->n_markers - 1) / (n_select - 1);

    for (int i = 0; i < n_select; i++) {
        int idx = (int)(i * step + 0.5);  /* round to nearest */
        if (idx >= cm->n_markers) idx = cm->n_markers - 1;
        cm->markers[idx].selected = 1;
    }
}

/* Select markers randomly */
static void select_random(chrom_markers_t* cm, int n_select, gsl_rng* rng) {
    if (n_select <= 0 || cm->n_markers == 0) return;
    if (n_select >= cm->n_markers) {
        /* Select all */
        for (int i = 0; i < cm->n_markers; i++) {
            cm->markers[i].selected = 1;
        }
        return;
    }

    /* Fisher-Yates shuffle to select n_select random markers */
    /* Create index array */
    int* indices = malloc(sizeof(int) * cm->n_markers);
    if (!indices) return;
    for (int i = 0; i < cm->n_markers; i++) {
        indices[i] = i;
    }

    /* Partial shuffle: only need first n_select elements */
    for (int i = 0; i < n_select; i++) {
        int j = i + gsl_rng_uniform_int(rng, cm->n_markers - i);
        int tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }

    /* Mark selected */
    for (int i = 0; i < n_select; i++) {
        cm->markers[indices[i]].selected = 1;
    }

    free(indices);
}

/* Parse VCF header line to get sample names and their column indices */
static int parse_header_line(const char* line, char*** sample_names, int* n_samples) {
    /* Header format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ... */
    /* Samples start at column 9 (0-indexed) */

    *sample_names = NULL;
    *n_samples = 0;

    char* copy = strdup(line);
    if (!copy) return -1;

    int capacity = 64;
    *sample_names = malloc(sizeof(char*) * capacity);
    if (!*sample_names) {
        free(copy);
        return -1;
    }

    int col = 0;
    char* token = strtok(copy, "\t");
    while (token) {
        if (col >= 9) {  /* Sample columns start at index 9 */
            if (*n_samples >= capacity) {
                capacity *= 2;
                char** new_names = realloc(*sample_names, sizeof(char*) * capacity);
                if (!new_names) {
                    for (int i = 0; i < *n_samples; i++) free((*sample_names)[i]);
                    free(*sample_names);
                    free(copy);
                    return -1;
                }
                *sample_names = new_names;
            }
            /* Remove trailing newline */
            size_t len = strlen(token);
            while (len > 0 && (token[len-1] == '\n' || token[len-1] == '\r')) {
                token[--len] = '\0';
            }
            (*sample_names)[*n_samples] = strdup(token);
            (*n_samples)++;
        }
        token = strtok(NULL, "\t");
        col++;
    }

    free(copy);
    return 0;
}

/* Write output line with filtered samples */
static void write_filtered_line(FILE* out, const char* line,
                                const int* keep_cols, int n_keep, int n_samples) {
    char* copy = strdup(line);
    if (!copy) return;

    int col = 0;
    char* start = copy;
    char* end;
    int sample_idx = 0;

    while (*start) {
        end = start;
        while (*end && *end != '\t' && *end != '\n' && *end != '\r') end++;
        char saved = *end;
        *end = '\0';

        if (col < 9) {
            /* Fixed columns: always output */
            if (col > 0) fprintf(out, "\t");
            fprintf(out, "%s", start);
        } else {
            /* Sample column: check if we keep it */
            if (sample_idx < n_samples && keep_cols[sample_idx]) {
                fprintf(out, "\t%s", start);
            }
            sample_idx++;
        }

        if (saved == '\0') break;
        start = end + 1;
        col++;
    }
    fprintf(out, "\n");

    free(copy);
}

int run_sample(sample_params_t* params) {
    /* Initialize RNG */
    gsl_rng* rng = NULL;
    if (params->n_markers > 0) {
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_default);
        if (params->seed == 0) {
            gsl_rng_set(rng, time(NULL));
        } else {
            gsl_rng_set(rng, params->seed);
        }
    }

    /* Load individual list if specified */
    indiv_list_t* indiv_list = NULL;
    if (params->indiv_file) {
        indiv_list = load_indiv_list(params->indiv_file);
        if (!indiv_list) {
            fprintf(stderr, "Error: Cannot load individual list from '%s'\n",
                    params->indiv_file);
            if (rng) gsl_rng_free(rng);
            return 1;
        }
        if (params->verbose) {
            fprintf(stderr, "Loaded %d individuals from list\n", indiv_list->n_ids);
        }
    }

    /* Open input file */
    FILE* in = fopen(params->input_file, "r");
    if (!in) {
        fprintf(stderr, "Error: Cannot open input file '%s'\n", params->input_file);
        if (indiv_list) free_indiv_list(indiv_list);
        if (rng) gsl_rng_free(rng);
        return 1;
    }

    /* First pass: read header and collect marker positions */
    char* line = malloc(MAX_LINE_LEN);
    if (!line) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(in);
        if (indiv_list) free_indiv_list(indiv_list);
        if (rng) gsl_rng_free(rng);
        return 1;
    }

    char** sample_names = NULL;
    int n_samples = 0;
    int* keep_cols = NULL;
    int n_keep = 0;
    marker_collection_t* markers = NULL;
    char** meta_lines = NULL;
    int n_meta = 0;
    int meta_capacity = 256;
    char* header_line = NULL;

    meta_lines = malloc(sizeof(char*) * meta_capacity);
    if (!meta_lines) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free(line);
        fclose(in);
        if (indiv_list) free_indiv_list(indiv_list);
        if (rng) gsl_rng_free(rng);
        return 1;
    }

    if (params->n_markers > 0) {
        markers = create_marker_collection();
        if (!markers) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            free(line);
            free(meta_lines);
            fclose(in);
            if (indiv_list) free_indiv_list(indiv_list);
            if (rng) gsl_rng_free(rng);
            return 1;
        }
    }

    /* Read header and collect markers */
    long file_offset;
    while ((file_offset = ftell(in)) >= 0 && fgets(line, MAX_LINE_LEN, in)) {
        if (line[0] == '#' && line[1] == '#') {
            /* Meta-information line - store for output */
            if (n_meta >= meta_capacity) {
                meta_capacity *= 2;
                char** new_meta = realloc(meta_lines, sizeof(char*) * meta_capacity);
                if (!new_meta) {
                    fprintf(stderr, "Error: Memory allocation failed\n");
                    goto cleanup_error;
                }
                meta_lines = new_meta;
            }
            meta_lines[n_meta++] = strdup(line);
        } else if (line[0] == '#') {
            /* Header line with sample names */
            header_line = strdup(line);
            if (parse_header_line(line, &sample_names, &n_samples) != 0) {
                fprintf(stderr, "Error: Failed to parse header line\n");
                goto cleanup_error;
            }

            /* Determine which columns to keep */
            keep_cols = malloc(sizeof(int) * n_samples);
            if (!keep_cols) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                goto cleanup_error;
            }

            for (int i = 0; i < n_samples; i++) {
                if (indiv_list) {
                    keep_cols[i] = indiv_in_list(indiv_list, sample_names[i]);
                } else {
                    keep_cols[i] = 1;  /* keep all */
                }
                if (keep_cols[i]) n_keep++;
            }

            if (params->verbose) {
                fprintf(stderr, "VCF contains %d samples, keeping %d\n",
                        n_samples, n_keep);
            }
        } else {
            /* Data line - collect marker info */
            if (markers) {
                /* Parse CHROM and POS */
                char chrom[256];
                long pos;
                if (sscanf(line, "%255s\t%ld", chrom, &pos) == 2) {
                    /* Check if in region (or no region specified) */
                    if (pos_in_region(&params->region, chrom, pos)) {
                        chrom_markers_t* cm = get_or_create_chrom(markers, chrom);
                        if (cm) {
                            add_marker(cm, pos, file_offset);
                        }
                    }
                }
            }
        }
    }

    if (params->verbose && markers) {
        int total = 0;
        for (int i = 0; i < markers->n_chroms; i++) {
            total += markers->chroms[i].n_markers;
        }
        fprintf(stderr, "Found %d markers in region across %d chromosomes\n",
                total, markers->n_chroms);
    }

    /* Select markers if requested */
    if (markers && params->n_markers > 0) {
        if (params->region.chrom) {
            /* Single chromosome specified - select from it */
            for (int i = 0; i < markers->n_chroms; i++) {
                if (strcmp(markers->chroms[i].chrom, params->region.chrom) == 0) {
                    if (params->method == METHOD_UNIFORM) {
                        select_uniform(&markers->chroms[i], params->n_markers, rng);
                    } else {
                        select_random(&markers->chroms[i], params->n_markers, rng);
                    }
                    break;
                }
            }
        } else {
            /* Apply to each chromosome independently */
            for (int i = 0; i < markers->n_chroms; i++) {
                if (params->method == METHOD_UNIFORM) {
                    select_uniform(&markers->chroms[i], params->n_markers, rng);
                } else {
                    select_random(&markers->chroms[i], params->n_markers, rng);
                }
            }
        }

        if (params->verbose) {
            int selected = 0;
            for (int i = 0; i < markers->n_chroms; i++) {
                for (int j = 0; j < markers->chroms[i].n_markers; j++) {
                    if (markers->chroms[i].markers[j].selected) selected++;
                }
            }
            fprintf(stderr, "Selected %d markers\n", selected);
        }
    }

    /* Open output file */
    FILE* out = stdout;
    if (params->output_file) {
        out = fopen(params->output_file, "w");
        if (!out) {
            fprintf(stderr, "Error: Cannot open output file '%s'\n",
                    params->output_file);
            goto cleanup_error;
        }
    }

    /* Write meta-information lines, filtering contig lines if subsetting markers */
    for (int i = 0; i < n_meta; i++) {
        /* Check if this is a contig line and if we're subsetting */
        if (markers && strncmp(meta_lines[i], "##contig=", 9) == 0) {
            /* Extract contig ID and check if it's in our markers */
            char* id_start = strstr(meta_lines[i], "ID=");
            if (id_start) {
                id_start += 3;
                char* id_end = strchr(id_start, ',');
                if (!id_end) id_end = strchr(id_start, '>');
                if (id_end) {
                    size_t len = id_end - id_start;
                    char contig_id[256];
                    if (len < sizeof(contig_id)) {
                        strncpy(contig_id, id_start, len);
                        contig_id[len] = '\0';

                        /* Check if this contig has selected markers */
                        int has_selected = 0;
                        for (int c = 0; c < markers->n_chroms; c++) {
                            if (strcmp(markers->chroms[c].chrom, contig_id) == 0) {
                                for (int m = 0; m < markers->chroms[c].n_markers; m++) {
                                    if (markers->chroms[c].markers[m].selected) {
                                        has_selected = 1;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        if (!has_selected) continue;  /* skip this contig */
                    }
                }
            }
        }
        fprintf(out, "%s", meta_lines[i]);
    }

    /* Write header line with filtered samples */
    if (header_line) {
        write_filtered_line(out, header_line, keep_cols, n_keep, n_samples);
    }

    /* Second pass: write data lines */
    rewind(in);
    while (fgets(line, MAX_LINE_LEN, in)) {
        if (line[0] == '#') continue;  /* skip header lines */

        /* Parse CHROM and POS */
        char chrom[256];
        long pos;
        if (sscanf(line, "%255s\t%ld", chrom, &pos) != 2) continue;

        /* Check if this marker is selected */
        int include = 1;
        if (markers) {
            include = 0;
            for (int c = 0; c < markers->n_chroms; c++) {
                if (strcmp(markers->chroms[c].chrom, chrom) == 0) {
                    for (int m = 0; m < markers->chroms[c].n_markers; m++) {
                        if (markers->chroms[c].markers[m].pos == pos &&
                            markers->chroms[c].markers[m].selected) {
                            include = 1;
                            break;
                        }
                    }
                    break;
                }
            }
        } else {
            /* No marker subsetting - check region only */
            include = pos_in_region(&params->region, chrom, pos);
        }

        if (include) {
            write_filtered_line(out, line, keep_cols, n_keep, n_samples);
        }
    }

    /* Cleanup */
    if (params->output_file && out) fclose(out);
    fclose(in);
    free(line);
    if (header_line) free(header_line);
    for (int i = 0; i < n_meta; i++) free(meta_lines[i]);
    free(meta_lines);
    if (sample_names) {
        for (int i = 0; i < n_samples; i++) free(sample_names[i]);
        free(sample_names);
    }
    if (keep_cols) free(keep_cols);
    if (markers) free_marker_collection(markers);
    if (indiv_list) free_indiv_list(indiv_list);
    if (rng) gsl_rng_free(rng);

    if (params->verbose) {
        fprintf(stderr, "Done.\n");
    }

    return 0;

cleanup_error:
    fclose(in);
    free(line);
    if (header_line) free(header_line);
    for (int i = 0; i < n_meta; i++) free(meta_lines[i]);
    free(meta_lines);
    if (sample_names) {
        for (int i = 0; i < n_samples; i++) free(sample_names[i]);
        free(sample_names);
    }
    if (keep_cols) free(keep_cols);
    if (markers) free_marker_collection(markers);
    if (indiv_list) free_indiv_list(indiv_list);
    if (rng) gsl_rng_free(rng);
    return 1;
}
