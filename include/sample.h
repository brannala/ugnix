/*
 * sample.h - VCF Subsampling Tool
 *
 * Subsample individuals and/or markers from a VCF file.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#ifndef SAMPLE_H
#define SAMPLE_H

#include <stdio.h>
#include <glib.h>

/* Marker selection methods */
typedef enum {
    METHOD_UNIFORM,   /* evenly spaced, maximizing inter-marker distance */
    METHOD_RANDOM     /* random selection with equal probability */
} sample_method_t;

/*
 * Region specification for marker subsetting.
 * If chrom is NULL, apply to all chromosomes independently.
 */
typedef struct {
    char* chrom;      /* chromosome/contig name (NULL = all) */
    long start;       /* start position (1-based, 0 = beginning) */
    long end;         /* end position (0 = end of chromosome) */
} vcf_region_t;

/*
 * Parameters for VCF sampling.
 */
typedef struct {
    /* Input/output - single file mode */
    const char* input_file;       /* input VCF file (single-file mode) */
    const char* output_file;      /* output VCF file (NULL = stdout) */

    /* Input/output - multi-file mode */
    const char** input_files;     /* array of input VCF files */
    int n_input_files;            /* number of input files (0 = single-file mode) */
    const char* output_suffix;    /* output suffix for multi-file mode (e.g., "_filtered") */

    /* Individual subsetting */
    const char* indiv_file;       /* file with individual IDs (NULL = keep all) */

    /* Marker subsetting */
    int n_markers;                /* number of markers to select (0 = keep all) */
    sample_method_t method;       /* selection method */
    vcf_region_t region;          /* region specification */

    /* Options */
    unsigned long seed;           /* random seed */
    int verbose;                  /* verbose output */
} sample_params_t;

/*
 * Individual ID list.
 */
typedef struct {
    char** ids;       /* array of individual IDs */
    int n_ids;        /* number of IDs */
    int capacity;     /* allocated capacity */
} indiv_list_t;

/*
 * Marker information for selection.
 */
typedef struct {
    char* chrom;      /* chromosome name */
    long pos;         /* position */
    long file_offset; /* byte offset in file for this record */
    int selected;     /* flag: 1 if selected, 0 otherwise */
} marker_info_t;

/*
 * Marker list for a single chromosome.
 */
typedef struct {
    char* chrom;          /* chromosome name */
    marker_info_t* markers;
    int n_markers;
    int capacity;
} chrom_markers_t;

/*
 * Collection of markers across all chromosomes.
 */
typedef struct {
    chrom_markers_t* chroms;
    int n_chroms;
    int capacity;
} marker_collection_t;

/*
 * Marker set for intersection operations (uses GLib hash table).
 * Key format: "chrom:pos" string
 */
typedef struct {
    GHashTable* table;
    int count;
} marker_set_t;

/*
 * Initialize default parameters.
 */
void init_sample_params(sample_params_t* params);

/*
 * Parse region string "chrom:start-end" or "chrom".
 * Returns 0 on success, non-zero on error.
 */
int parse_region(const char* region_str, vcf_region_t* region);

/*
 * Free region resources.
 */
void free_region(vcf_region_t* region);

/*
 * Load individual IDs from file.
 * Returns NULL on error.
 */
indiv_list_t* load_indiv_list(const char* filename);

/*
 * Free individual list.
 */
void free_indiv_list(indiv_list_t* list);

/*
 * Check if an individual ID is in the list.
 * Returns 1 if found, 0 otherwise.
 */
int indiv_in_list(const indiv_list_t* list, const char* id);

/*
 * Run the VCF sampling (single-file mode).
 * Returns 0 on success, non-zero on error.
 */
int run_sample(sample_params_t* params);

/*
 * Run the VCF sampling (multi-file mode with marker intersection).
 * Returns 0 on success, non-zero on error.
 */
int run_sample_multi(sample_params_t* params);

/*
 * Collect markers from a VCF file without individual filtering.
 * Returns marker collection, or NULL on error.
 */
marker_collection_t* collect_markers_from_vcf(const char* filename,
                                              const vcf_region_t* region,
                                              int verbose);

/*
 * Create marker set from collection.
 */
marker_set_t* markers_to_set(marker_collection_t* mc);

/*
 * Intersect set with collection, removing markers not in collection.
 */
void intersect_with_collection(marker_set_t* set, marker_collection_t* mc);

/*
 * Check if marker is in set.
 */
int marker_in_set(marker_set_t* set, const char* chrom, long pos);

/*
 * Convert set back to marker collection for selection.
 */
marker_collection_t* set_to_collection(marker_set_t* set);

/*
 * Free marker set.
 */
void free_marker_set(marker_set_t* set);

/*
 * Generate output filename from input file and suffix.
 * Caller must free the returned string.
 */
char* generate_output_filename(const char* input_file, const char* suffix);

#endif /* SAMPLE_H */
