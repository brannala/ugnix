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
#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "sample.h"

#define MAX_LINE_LEN 1048576  /* 1MB max line length */
#define INITIAL_CAPACITY 1024

void init_sample_params(sample_params_t* params) {
    params->input_file = NULL;
    params->output_file = NULL;
    params->input_files = NULL;
    params->n_input_files = 0;
    params->output_suffix = NULL;
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

/*
 * Multi-file marker intersection support
 */

/* Create string key "chrom:pos" for hash table */
static char* make_marker_key(const char* chrom, long pos) {
    char* key = malloc(strlen(chrom) + 32);
    if (!key) return NULL;
    sprintf(key, "%s:%ld", chrom, pos);
    return key;
}

/* Create marker set from collection */
marker_set_t* markers_to_set(marker_collection_t* mc) {
    if (!mc) return NULL;

    marker_set_t* set = malloc(sizeof(marker_set_t));
    if (!set) return NULL;

    set->table = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    if (!set->table) {
        free(set);
        return NULL;
    }
    set->count = 0;

    for (int c = 0; c < mc->n_chroms; c++) {
        chrom_markers_t* cm = &mc->chroms[c];
        for (int m = 0; m < cm->n_markers; m++) {
            char* key = make_marker_key(cm->chrom, cm->markers[m].pos);
            if (key && !g_hash_table_contains(set->table, key)) {
                g_hash_table_insert(set->table, key, GINT_TO_POINTER(1));
                set->count++;
            } else if (key) {
                free(key);  /* duplicate */
            }
        }
    }

    return set;
}

/* Intersect set with collection, removing markers not in collection */
void intersect_with_collection(marker_set_t* set, marker_collection_t* mc) {
    if (!set || !mc) return;

    /* Build hash of markers in collection */
    GHashTable* mc_set = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);

    for (int c = 0; c < mc->n_chroms; c++) {
        chrom_markers_t* cm = &mc->chroms[c];
        for (int m = 0; m < cm->n_markers; m++) {
            char* key = make_marker_key(cm->chrom, cm->markers[m].pos);
            if (key) {
                g_hash_table_insert(mc_set, key, GINT_TO_POINTER(1));
            }
        }
    }

    /* Remove from set any markers not in mc_set */
    GHashTableIter iter;
    gpointer key, value;
    GList* to_remove = NULL;

    g_hash_table_iter_init(&iter, set->table);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        if (!g_hash_table_contains(mc_set, key)) {
            to_remove = g_list_prepend(to_remove, g_strdup((char*)key));
        }
    }

    /* Remove the markers not in intersection */
    for (GList* l = to_remove; l != NULL; l = l->next) {
        g_hash_table_remove(set->table, l->data);
        set->count--;
        g_free(l->data);
    }
    g_list_free(to_remove);

    g_hash_table_destroy(mc_set);
}

/* Check if marker is in set */
int marker_in_set(marker_set_t* set, const char* chrom, long pos) {
    if (!set || !chrom) return 0;

    char* key = make_marker_key(chrom, pos);
    if (!key) return 0;

    int result = g_hash_table_contains(set->table, key);
    free(key);
    return result;
}

/* Free marker set */
void free_marker_set(marker_set_t* set) {
    if (!set) return;
    if (set->table) {
        g_hash_table_destroy(set->table);
    }
    free(set);
}

/* Helper: comparison function for sorting markers by chrom then pos */
static int compare_marker_keys(const void* a, const void* b) {
    const char* key_a = *(const char**)a;
    const char* key_b = *(const char**)b;

    /* Parse chrom:pos from keys */
    char chrom_a[256], chrom_b[256];
    long pos_a, pos_b;

    sscanf(key_a, "%255[^:]:%ld", chrom_a, &pos_a);
    sscanf(key_b, "%255[^:]:%ld", chrom_b, &pos_b);

    int cmp = strcmp(chrom_a, chrom_b);
    if (cmp != 0) return cmp;
    return (pos_a > pos_b) - (pos_a < pos_b);
}

/* Convert set back to marker collection for selection */
marker_collection_t* set_to_collection(marker_set_t* set) {
    if (!set || set->count == 0) return NULL;

    /* Get all keys and sort them */
    guint size = g_hash_table_size(set->table);
    char** keys = malloc(sizeof(char*) * size);
    if (!keys) return NULL;

    GHashTableIter iter;
    gpointer key, value;
    int i = 0;

    g_hash_table_iter_init(&iter, set->table);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        keys[i++] = (char*)key;
    }

    qsort(keys, size, sizeof(char*), compare_marker_keys);

    /* Create marker collection */
    marker_collection_t* mc = malloc(sizeof(marker_collection_t));
    if (!mc) {
        free(keys);
        return NULL;
    }

    mc->capacity = 64;
    mc->n_chroms = 0;
    mc->chroms = malloc(sizeof(chrom_markers_t) * mc->capacity);
    if (!mc->chroms) {
        free(mc);
        free(keys);
        return NULL;
    }

    /* Add markers to collection */
    char current_chrom[256] = "";
    chrom_markers_t* current_cm = NULL;

    for (guint j = 0; j < size; j++) {
        char chrom[256];
        long pos;
        sscanf(keys[j], "%255[^:]:%ld", chrom, &pos);

        /* New chromosome? */
        if (strcmp(chrom, current_chrom) != 0) {
            /* Expand if needed */
            if (mc->n_chroms >= mc->capacity) {
                mc->capacity *= 2;
                chrom_markers_t* new_chroms = realloc(mc->chroms,
                    sizeof(chrom_markers_t) * mc->capacity);
                if (!new_chroms) {
                    /* cleanup and return partial */
                    break;
                }
                mc->chroms = new_chroms;
            }

            current_cm = &mc->chroms[mc->n_chroms];
            current_cm->chrom = strdup(chrom);
            current_cm->capacity = INITIAL_CAPACITY;
            current_cm->n_markers = 0;
            current_cm->markers = malloc(sizeof(marker_info_t) * current_cm->capacity);
            mc->n_chroms++;
            strcpy(current_chrom, chrom);
        }

        /* Add marker */
        if (current_cm) {
            if (current_cm->n_markers >= current_cm->capacity) {
                current_cm->capacity *= 2;
                marker_info_t* new_markers = realloc(current_cm->markers,
                    sizeof(marker_info_t) * current_cm->capacity);
                if (!new_markers) break;
                current_cm->markers = new_markers;
            }

            marker_info_t* m = &current_cm->markers[current_cm->n_markers];
            m->chrom = current_cm->chrom;  /* shared pointer */
            m->pos = pos;
            m->file_offset = 0;
            m->selected = 0;
            current_cm->n_markers++;
        }
    }

    free(keys);
    return mc;
}

/* Generate output filename from input file and suffix */
char* generate_output_filename(const char* input_file, const char* suffix) {
    if (!input_file || !suffix) return NULL;

    /* Find the base name (remove directory) */
    const char* basename = strrchr(input_file, '/');
    if (basename) {
        basename++;  /* skip the '/' */
    } else {
        basename = input_file;
    }

    /* Find extension (.vcf or .vcf.gz) */
    const char* ext = NULL;
    const char* vcf_gz = strstr(basename, ".vcf.gz");
    const char* vcf = strstr(basename, ".vcf");

    if (vcf_gz) {
        ext = vcf_gz;
    } else if (vcf) {
        ext = vcf;
    }

    size_t base_len;
    if (ext) {
        base_len = ext - basename;
    } else {
        base_len = strlen(basename);
    }

    /* Allocate output: base + suffix + .vcf + null */
    size_t out_len = base_len + strlen(suffix) + 5;
    char* output = malloc(out_len);
    if (!output) return NULL;

    /* Copy base name */
    strncpy(output, basename, base_len);
    output[base_len] = '\0';

    /* Append suffix and extension */
    strcat(output, suffix);
    strcat(output, ".vcf");

    return output;
}

/* Collect markers from a VCF file (first pass only) */
marker_collection_t* collect_markers_from_vcf(const char* filename,
                                              const vcf_region_t* region,
                                              int verbose) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file '%s'\n", filename);
        return NULL;
    }

    marker_collection_t* mc = malloc(sizeof(marker_collection_t));
    if (!mc) {
        fclose(fp);
        return NULL;
    }

    mc->capacity = 64;
    mc->n_chroms = 0;
    mc->chroms = malloc(sizeof(chrom_markers_t) * mc->capacity);
    if (!mc->chroms) {
        free(mc);
        fclose(fp);
        return NULL;
    }

    char* line = malloc(MAX_LINE_LEN);
    if (!line) {
        free(mc->chroms);
        free(mc);
        fclose(fp);
        return NULL;
    }

    int marker_count = 0;
    while (fgets(line, MAX_LINE_LEN, fp)) {
        if (line[0] == '#') continue;  /* skip header */

        /* Parse CHROM and POS */
        char chrom[256];
        long pos;
        if (sscanf(line, "%255s\t%ld", chrom, &pos) != 2) continue;

        /* Check region filter */
        if (region && region->chrom) {
            if (strcmp(region->chrom, chrom) != 0) continue;
            if (region->start > 0 && pos < region->start) continue;
            if (region->end > 0 && pos > region->end) continue;
        }

        /* Find or create chromosome */
        chrom_markers_t* cm = NULL;
        for (int c = 0; c < mc->n_chroms; c++) {
            if (strcmp(mc->chroms[c].chrom, chrom) == 0) {
                cm = &mc->chroms[c];
                break;
            }
        }

        if (!cm) {
            /* Create new chromosome */
            if (mc->n_chroms >= mc->capacity) {
                mc->capacity *= 2;
                chrom_markers_t* new_chroms = realloc(mc->chroms,
                    sizeof(chrom_markers_t) * mc->capacity);
                if (!new_chroms) break;
                mc->chroms = new_chroms;
            }

            cm = &mc->chroms[mc->n_chroms];
            cm->chrom = strdup(chrom);
            cm->capacity = INITIAL_CAPACITY;
            cm->n_markers = 0;
            cm->markers = malloc(sizeof(marker_info_t) * cm->capacity);
            mc->n_chroms++;
        }

        /* Add marker */
        if (cm->n_markers >= cm->capacity) {
            cm->capacity *= 2;
            marker_info_t* new_markers = realloc(cm->markers,
                sizeof(marker_info_t) * cm->capacity);
            if (!new_markers) break;
            cm->markers = new_markers;
        }

        marker_info_t* m = &cm->markers[cm->n_markers];
        m->chrom = cm->chrom;
        m->pos = pos;
        m->file_offset = 0;
        m->selected = 0;
        cm->n_markers++;
        marker_count++;
    }

    free(line);
    fclose(fp);

    if (verbose) {
        fprintf(stderr, "  %s: %d markers", filename, marker_count);
        if (mc->n_chroms > 0) {
            fprintf(stderr, " across %d chromosomes", mc->n_chroms);
        }
        fprintf(stderr, "\n");
    }

    return mc;
}

/* Helper: write a single VCF file with selected markers */
static int write_filtered_vcf(const char* input_file, const char* output_file,
                              marker_set_t* selected, indiv_list_t* indiv_list,
                              int verbose) {
    FILE* in = fopen(input_file, "r");
    if (!in) {
        fprintf(stderr, "Error: Cannot open '%s'\n", input_file);
        return 1;
    }

    FILE* out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "Error: Cannot create '%s'\n", output_file);
        fclose(in);
        return 1;
    }

    char* line = malloc(MAX_LINE_LEN);
    if (!line) {
        fclose(in);
        fclose(out);
        return 1;
    }

    char** sample_names = NULL;
    int n_samples = 0;
    int* keep_cols = NULL;
    int n_keep = 0;
    int markers_written = 0;

    while (fgets(line, MAX_LINE_LEN, in)) {
        if (line[0] == '#' && line[1] == '#') {
            /* Meta line - pass through */
            fprintf(out, "%s", line);
        } else if (line[0] == '#') {
            /* Header line - parse samples and filter if needed */
            if (indiv_list) {
                /* Parse sample names */
                char* copy = strdup(line);
                int col = 0;
                int capacity = 64;
                sample_names = malloc(sizeof(char*) * capacity);
                n_samples = 0;

                char* token = strtok(copy, "\t");
                while (token) {
                    if (col >= 9) {
                        if (n_samples >= capacity) {
                            capacity *= 2;
                            sample_names = realloc(sample_names, sizeof(char*) * capacity);
                        }
                        /* Remove newline */
                        size_t len = strlen(token);
                        while (len > 0 && (token[len-1] == '\n' || token[len-1] == '\r')) {
                            token[--len] = '\0';
                        }
                        sample_names[n_samples++] = strdup(token);
                    }
                    token = strtok(NULL, "\t");
                    col++;
                }
                free(copy);

                /* Determine which columns to keep */
                keep_cols = malloc(sizeof(int) * n_samples);
                for (int i = 0; i < n_samples; i++) {
                    keep_cols[i] = indiv_in_list(indiv_list, sample_names[i]);
                    if (keep_cols[i]) n_keep++;
                }

                /* Write filtered header */
                copy = strdup(line);
                col = 0;
                int sample_idx = 0;
                char* start = copy;
                char* end;

                while (*start) {
                    end = start;
                    while (*end && *end != '\t' && *end != '\n' && *end != '\r') end++;
                    char saved = *end;
                    *end = '\0';

                    if (col < 9) {
                        if (col > 0) fprintf(out, "\t");
                        fprintf(out, "%s", start);
                    } else {
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
            } else {
                /* No individual filtering, pass through */
                fprintf(out, "%s", line);
            }
        } else {
            /* Data line - check if marker is selected */
            char chrom[256];
            long pos;
            if (sscanf(line, "%255s\t%ld", chrom, &pos) == 2) {
                if (marker_in_set(selected, chrom, pos)) {
                    if (keep_cols) {
                        /* Filter individuals */
                        char* copy = strdup(line);
                        int col = 0;
                        int sample_idx = 0;
                        char* start = copy;
                        char* end;

                        while (*start) {
                            end = start;
                            while (*end && *end != '\t' && *end != '\n' && *end != '\r') end++;
                            char saved = *end;
                            *end = '\0';

                            if (col < 9) {
                                if (col > 0) fprintf(out, "\t");
                                fprintf(out, "%s", start);
                            } else {
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
                    } else {
                        /* No individual filtering */
                        fprintf(out, "%s", line);
                    }
                    markers_written++;
                }
            }
        }
    }

    free(line);
    if (sample_names) {
        for (int i = 0; i < n_samples; i++) free(sample_names[i]);
        free(sample_names);
    }
    if (keep_cols) free(keep_cols);
    fclose(in);
    fclose(out);

    if (verbose) {
        fprintf(stderr, "  Wrote %d markers to %s\n", markers_written, output_file);
    }

    return 0;
}

/* Multi-file sampling with marker intersection */
int run_sample_multi(sample_params_t* params) {
    if (!params || params->n_input_files < 2) {
        fprintf(stderr, "Error: Multi-file mode requires at least 2 input files\n");
        return 1;
    }

    int n_files = params->n_input_files;
    int verbose = params->verbose;

    /* Initialize RNG */
    gsl_rng* rng = NULL;
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    if (params->seed == 0) {
        gsl_rng_set(rng, time(NULL));
    } else {
        gsl_rng_set(rng, params->seed);
    }

    /* Load individual list if specified */
    indiv_list_t* indiv_list = NULL;
    if (params->indiv_file) {
        indiv_list = load_indiv_list(params->indiv_file);
        if (!indiv_list) {
            fprintf(stderr, "Error: Cannot load individual list from '%s'\n",
                    params->indiv_file);
            gsl_rng_free(rng);
            return 1;
        }
        if (verbose) {
            fprintf(stderr, "Loaded %d individuals from list\n", indiv_list->n_ids);
        }
    }

    if (verbose) {
        fprintf(stderr, "Collecting markers from %d files...\n", n_files);
    }

    /* Collect markers from all files */
    marker_collection_t** collections = malloc(sizeof(marker_collection_t*) * n_files);
    if (!collections) {
        if (indiv_list) free_indiv_list(indiv_list);
        gsl_rng_free(rng);
        return 1;
    }

    for (int i = 0; i < n_files; i++) {
        collections[i] = collect_markers_from_vcf(params->input_files[i],
                                                   &params->region, verbose);
        if (!collections[i]) {
            fprintf(stderr, "Error: Failed to read markers from '%s'\n",
                    params->input_files[i]);
            for (int j = 0; j < i; j++) {
                free_marker_collection(collections[j]);
            }
            free(collections);
            if (indiv_list) free_indiv_list(indiv_list);
            gsl_rng_free(rng);
            return 1;
        }
    }

    /* Compute intersection */
    if (verbose) {
        fprintf(stderr, "Computing marker intersection...\n");
    }

    marker_set_t* intersection = markers_to_set(collections[0]);
    if (!intersection) {
        fprintf(stderr, "Error: Failed to create marker set\n");
        for (int i = 0; i < n_files; i++) {
            free_marker_collection(collections[i]);
        }
        free(collections);
        if (indiv_list) free_indiv_list(indiv_list);
        gsl_rng_free(rng);
        return 1;
    }

    for (int i = 1; i < n_files; i++) {
        intersect_with_collection(intersection, collections[i]);
    }

    if (verbose) {
        fprintf(stderr, "Intersection: %d markers\n", intersection->count);
    }

    if (intersection->count == 0) {
        fprintf(stderr, "Error: No markers shared across all input files\n");
        free_marker_set(intersection);
        for (int i = 0; i < n_files; i++) {
            free_marker_collection(collections[i]);
        }
        free(collections);
        if (indiv_list) free_indiv_list(indiv_list);
        gsl_rng_free(rng);
        return 1;
    }

    /* Convert intersection to collection for selection */
    marker_collection_t* common = set_to_collection(intersection);
    if (!common) {
        fprintf(stderr, "Error: Failed to convert marker set to collection\n");
        free_marker_set(intersection);
        for (int i = 0; i < n_files; i++) {
            free_marker_collection(collections[i]);
        }
        free(collections);
        if (indiv_list) free_indiv_list(indiv_list);
        gsl_rng_free(rng);
        return 1;
    }

    /* Apply selection to intersection */
    if (verbose) {
        fprintf(stderr, "Selecting %d markers per chromosome (%s)...\n",
                params->n_markers,
                params->method == METHOD_UNIFORM ? "uniform" : "random");
    }

    int total_selected = 0;
    for (int c = 0; c < common->n_chroms; c++) {
        chrom_markers_t* cm = &common->chroms[c];
        int n_select = params->n_markers;

        /* Warn if fewer markers available than requested */
        if (cm->n_markers < n_select) {
            if (verbose) {
                fprintf(stderr, "  Warning: %s has only %d markers (requested %d)\n",
                        cm->chrom, cm->n_markers, n_select);
            }
            n_select = cm->n_markers;
        }

        if (params->method == METHOD_UNIFORM) {
            /* Uniform selection */
            if (n_select > 0 && cm->n_markers > 0) {
                if (n_select >= cm->n_markers) {
                    for (int m = 0; m < cm->n_markers; m++) {
                        cm->markers[m].selected = 1;
                    }
                    total_selected += cm->n_markers;
                } else {
                    double step = (double)(cm->n_markers - 1) / (n_select - 1);
                    for (int i = 0; i < n_select; i++) {
                        int idx = (int)(i * step + 0.5);
                        if (idx >= cm->n_markers) idx = cm->n_markers - 1;
                        cm->markers[idx].selected = 1;
                    }
                    total_selected += n_select;
                }
            }
        } else {
            /* Random selection */
            if (n_select > 0 && cm->n_markers > 0) {
                if (n_select >= cm->n_markers) {
                    for (int m = 0; m < cm->n_markers; m++) {
                        cm->markers[m].selected = 1;
                    }
                    total_selected += cm->n_markers;
                } else {
                    int* indices = malloc(sizeof(int) * cm->n_markers);
                    for (int m = 0; m < cm->n_markers; m++) {
                        indices[m] = m;
                    }
                    /* Partial Fisher-Yates shuffle */
                    for (int i = 0; i < n_select; i++) {
                        int j = i + gsl_rng_uniform_int(rng, cm->n_markers - i);
                        int tmp = indices[i];
                        indices[i] = indices[j];
                        indices[j] = tmp;
                    }
                    for (int i = 0; i < n_select; i++) {
                        cm->markers[indices[i]].selected = 1;
                    }
                    free(indices);
                    total_selected += n_select;
                }
            }
        }
    }

    if (verbose) {
        fprintf(stderr, "Total selected: %d markers\n", total_selected);
    }

    /* Create marker set of selected markers for fast lookup */
    marker_set_t* selected = malloc(sizeof(marker_set_t));
    selected->table = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    selected->count = 0;

    for (int c = 0; c < common->n_chroms; c++) {
        chrom_markers_t* cm = &common->chroms[c];
        for (int m = 0; m < cm->n_markers; m++) {
            if (cm->markers[m].selected) {
                char* key = make_marker_key(cm->chrom, cm->markers[m].pos);
                if (key) {
                    g_hash_table_insert(selected->table, key, GINT_TO_POINTER(1));
                    selected->count++;
                }
            }
        }
    }

    /* Write filtered output files */
    if (verbose) {
        fprintf(stderr, "Writing output files...\n");
    }

    int error = 0;
    for (int i = 0; i < n_files && !error; i++) {
        char* output_file = generate_output_filename(params->input_files[i],
                                                      params->output_suffix);
        if (!output_file) {
            fprintf(stderr, "Error: Failed to generate output filename\n");
            error = 1;
            break;
        }

        error = write_filtered_vcf(params->input_files[i], output_file,
                                   selected, indiv_list, verbose);
        free(output_file);
    }

    /* Cleanup */
    free_marker_set(selected);
    free_marker_collection(common);
    free_marker_set(intersection);
    for (int i = 0; i < n_files; i++) {
        free_marker_collection(collections[i]);
    }
    free(collections);
    if (indiv_list) free_indiv_list(indiv_list);
    gsl_rng_free(rng);

    if (verbose && !error) {
        fprintf(stderr, "Done.\n");
    }

    return error;
}
