/*
 * pedtrans.h - Pedigree Chromosome Transmission Simulator
 *
 * Data structures and function declarations for simulating chromosome
 * transmission through a pedigree with recombination, tracking founder
 * chromosome segments.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#ifndef PEDTRANS_H
#define PEDTRANS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Maximum length for individual names */
#define PED_MAX_NAME 64

/* Initial capacity for dynamic arrays */
#define PED_INITIAL_CAPACITY 1024

/* ============================================================================
 * PEDIGREE DATA STRUCTURES
 * ============================================================================ */

/*
 * Pedigree individual - stored in array for O(1) access
 */
typedef struct {
    int id;              /* Index in array (0-based) */
    int father_id;       /* Father's index (-1 if founder) */
    int mother_id;       /* Mother's index (-1 if founder) */
    int is_founder;      /* 1 if founder, 0 otherwise */
    int in_degree;       /* For topological sort: # unprocessed parents */
} ped_indiv;

/*
 * Efficient pedigree structure with O(1) lookup
 */
typedef struct {
    ped_indiv* individuals;    /* Array of individuals (contiguous memory) */
    int* topo_order;           /* Topologically sorted indices */
    int* children;             /* Flat array of child indices for topo sort */
    int* child_offsets;        /* Offset into children array for each individual */
    int* n_children;           /* Number of children for each individual */
    GHashTable* name_to_id;    /* O(1) name lookup: char* -> int* */
    char** id_to_name;         /* O(1) index to name lookup */
    int n_individuals;         /* Current number of individuals */
    int n_founders;            /* Number of founders */
    int capacity;              /* Allocated capacity */
} pedigree;

/* ============================================================================
 * CHROMOSOME SEGMENT DATA STRUCTURES
 * ============================================================================ */

/*
 * Founder chromosome identity
 * Identifies which founder and which of their two chromosomes
 */
typedef struct {
    int founder_id;    /* Index of founder (0-based among founders) */
    int homolog;       /* 0 = paternal, 1 = maternal chromosome of founder */
} founder_origin;

/*
 * Chromosome segment with single founder origin
 * Segments are stored in arrays for cache efficiency
 */
typedef struct {
    double start;           /* Start position [0, 1) */
    double end;             /* End position (start, 1] */
    founder_origin origin;  /* Which founder chromosome this came from */
} ped_segment;

/*
 * Single chromosome as dynamic array of segments
 */
typedef struct {
    ped_segment* segments;  /* Array of segments (sorted by position) */
    int n_segments;         /* Current number of segments */
    int capacity;           /* Allocated capacity */
} ped_chromosome;

/*
 * Diploid individual's pair of homologous chromosomes
 */
typedef struct {
    ped_chromosome* paternal;  /* Chromosome inherited from father */
    ped_chromosome* maternal;  /* Chromosome inherited from mother */
} diploid_chromosomes;

/*
 * Complete simulation state
 */
typedef struct {
    pedigree* ped;                     /* Pedigree structure */
    diploid_chromosomes* chromosomes;  /* Array indexed by individual ID */
    double rec_rate;                   /* Recombination rate (expected crossovers) */
    gsl_rng* rng;                      /* Random number generator */
    unsigned long seed;                /* RNG seed for reproducibility */
} ped_simulation;

/* ============================================================================
 * PEDIGREE FUNCTIONS
 * ============================================================================ */

/*
 * Create a new pedigree with given initial capacity
 * Returns NULL on allocation failure
 */
pedigree* create_pedigree(int initial_capacity);

/*
 * Free pedigree and all associated memory
 */
void free_pedigree(pedigree* ped);

/*
 * Parse pedigree from file
 * File format: "indiv father mother" per line (founders have father=mother="0")
 * Returns NULL on error
 */
pedigree* parse_pedigree(FILE* fp);

/*
 * Parse pedigree from filename
 * Convenience wrapper around parse_pedigree
 */
pedigree* parse_pedigree_file(const char* filename);

/*
 * Add individual to pedigree
 * Returns individual's ID on success, -1 on error
 */
int pedigree_add_individual(pedigree* ped, const char* name,
                            const char* father, const char* mother);

/*
 * Get individual ID by name
 * Returns -1 if not found
 */
int pedigree_get_id(pedigree* ped, const char* name);

/*
 * Get individual name by ID
 * Returns NULL if ID is invalid
 */
const char* pedigree_get_name(pedigree* ped, int id);

/*
 * Perform topological sort of pedigree (founders first, then descendants)
 * Uses Kahn's algorithm - O(n)
 * Populates ped->topo_order
 * Returns 0 on success, -1 on cycle detected (invalid pedigree)
 */
int topological_sort(pedigree* ped);

/*
 * Build child list for each individual (needed for topological sort)
 * Called automatically by topological_sort
 */
void build_child_lists(pedigree* ped);

/* ============================================================================
 * CHROMOSOME FUNCTIONS
 * ============================================================================ */

/*
 * Create a new chromosome with given initial segment capacity
 */
ped_chromosome* create_chromosome(int initial_capacity);

/*
 * Create a founder chromosome (single segment spanning [0,1])
 */
ped_chromosome* create_founder_chromosome(int founder_id, int homolog);

/*
 * Free chromosome and all segments
 */
void free_chromosome(ped_chromosome* chr);

/*
 * Add segment to chromosome (at end - caller must ensure sorted order)
 */
void chromosome_add_segment(ped_chromosome* chr, double start, double end,
                            founder_origin origin);

/*
 * Copy segments from source that overlap [range_start, range_end)
 * Clips segments at boundaries
 */
void copy_segments_in_range(ped_chromosome* dest, ped_chromosome* source,
                            double range_start, double range_end);

/*
 * Merge adjacent segments with identical founder origin
 */
void merge_adjacent_segments(ped_chromosome* chr);

/*
 * Deep copy a chromosome
 */
ped_chromosome* copy_chromosome(ped_chromosome* src);

/* ============================================================================
 * MEIOSIS AND SIMULATION FUNCTIONS
 * ============================================================================ */

/*
 * Perform meiosis: generate gamete from diploid parent
 * Uses Poisson model for crossovers
 * Returns new chromosome (caller owns memory)
 */
ped_chromosome* meiosis(ped_chromosome* parent_pat, ped_chromosome* parent_mat,
                        double rec_rate, gsl_rng* rng);

/*
 * Create simulation state
 */
ped_simulation* create_simulation(void);

/*
 * Free simulation and all associated memory
 */
void free_simulation(ped_simulation* sim);

/*
 * Run simulation on pedigree file
 * rec_rate: expected number of crossovers per meiosis
 * seed: RNG seed (0 for time-based)
 */
ped_simulation* simulate_pedigree(const char* ped_file, double rec_rate,
                                   unsigned long seed);

/* ============================================================================
 * OUTPUT FUNCTIONS
 * ============================================================================ */

/*
 * Print segments for one individual
 */
void print_individual_segments(ped_simulation* sim, int indiv_id, FILE* out);

/*
 * Print segments for all individuals
 */
void print_all_segments(ped_simulation* sim, FILE* out);

/*
 * Print simulation summary/header
 */
void print_simulation_header(ped_simulation* sim, FILE* out);

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

/*
 * Compare two founder origins for equality
 */
static inline int origin_equal(founder_origin a, founder_origin b) {
    return (a.founder_id == b.founder_id) && (a.homolog == b.homolog);
}

/*
 * Compare doubles for sorting
 */
int compare_double(const void* a, const void* b);

#endif /* PEDTRANS_H */
