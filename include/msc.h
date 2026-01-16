#ifndef MSC_H
#define MSC_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "species_tree.h"
#include "coalescent.h"

/*
 * msc.h - Multispecies Coalescent Simulation
 *
 * Simulates gene genealogies under the multispecies coalescent model.
 * Lineages can only coalesce within the same population, and populations
 * merge at divergence times (going backwards) according to the species tree.
 */

/* ============== Migration Structures ============== */

/*
 * A single migration band between two populations.
 * Migration is specified going backwards in time: a lineage in to_pop
 * migrates to from_pop at rate M/2 per lineage.
 * Going forward in time, this represents individuals moving from from_pop to to_pop.
 */
typedef struct {
    int from_pop;           /* Source population ID (backwards: destination) */
    int to_pop;             /* Dest population ID (backwards: source of lineage) */
    double M;               /* Migration rate M = 4Nm */
} migration_band;

/*
 * List of all migration bands.
 */
typedef struct {
    migration_band* bands;
    int n_bands;
} migration_band_list;

/* ============== MSC Parameters ============== */

typedef struct {
    species_tree* tree;              /* Parsed species tree */
    double recombination_rate;       /* Global recombination rate (cM) */
    double mutation_rate;            /* Mutation rate (optional, for seq output) */
    unsigned long seed;              /* RNG seed (0 for time-based) */

    /* Output options */
    int output_gene_trees;           /* Print gene trees in Newick format */
    int output_vcf;                  /* Generate VCF output */
    char* vcf_filename;              /* VCF output filename */
    int verbose;                     /* Verbose output */

    /* Derived parameters (computed after parsing) */
    int total_samples;               /* Sum of sample sizes across species */

    /* Migration parameters */
    migration_band_list* migrations; /* Migration bands (NULL if no migration) */
} msc_params;

/* ============== Control File Parsing ============== */

/*
 * Parse MSC control file.
 *
 * Control file format:
 *   # Comment lines start with #
 *   species_tree: ((A:1000#0.001,B:1000#0.002):500#0.003,C:1500#0.001)#0.002;
 *   samples: A:10,B:10,C:10
 *   recombination_rate: 1.0
 *   mutation_rate: 0.5
 *   seed: 12345
 *   output_gene_trees: true
 *   output_vcf: true
 *   vcf_file: output.vcf
 *
 * Returns msc_params on success, NULL on error.
 */
msc_params* msc_parse_control_file(const char* filename);

/*
 * Validate parsed parameters (check all required fields are set).
 * Returns 0 on success, -1 on error (with message to stderr).
 */
int msc_validate_params(msc_params* params);

/*
 * Free MSC parameters.
 */
void msc_free_params(msc_params* params);

/* ============== MSC Simulation ============== */

/*
 * Per-population rate information during simulation.
 */
typedef struct {
    int pop_id;                      /* Population ID */
    int n_lineages;                  /* Number of active lineages in population */
    double coal_rate;                /* Coalescence rate for this population */
} pop_rate_info;

/*
 * Per-migration-band rate information during simulation.
 */
typedef struct {
    int from_pop;                    /* Source population (backwards: destination) */
    int to_pop;                      /* Destination pop (backwards: source of lineage) */
    double rate;                     /* n_to * M / 2 */
} mig_rate_info;

/*
 * Rate calculation result.
 */
typedef struct {
    pop_rate_info* pop_rates;        /* Per-population rates */
    int n_active_pops;               /* Number of populations with >= 2 lineages */
    double total_coal_rate;          /* Sum of coalescence rates */

    /* Migration rates */
    mig_rate_info* mig_rates;        /* Per-band migration rates */
    int n_mig_pairs;                 /* Number of active migration pairs */
    double total_mig_rate;           /* Sum of migration rates */

    double rec_rate;                 /* Global recombination rate * active length */
    double mut_rate;                 /* Global mutation rate * active length */
    double total_rate;               /* Sum of all rates */
} msc_rates;

/*
 * Create initial sample with population assignments.
 * Each chromosome is tagged with its species' population_id.
 */
chrsample* msc_create_sample(species_tree* tree, int total_samples);

/*
 * Calculate event rates for MSC simulation.
 * Per-population coalescence: n_p(n_p-1)/(2*theta_p)
 * Per-band migration: n_to * M / 2
 */
msc_rates* msc_calculate_rates(chrsample* sample, species_tree* tree,
                                msc_params* params, double current_time,
                                const bitarray* mrca);

/*
 * Free rates structure.
 */
void msc_free_rates(msc_rates* rates);

/*
 * Select which population has the coalescence event.
 * Returns population ID, weighted by each population's coalescence rate.
 */
int msc_select_coal_population(gsl_rng* r, msc_rates* rates);

/*
 * Select a coalescence pair from a specific population.
 * Only chromosomes with matching population_id can be selected.
 * Returns 0 on success, -1 if < 2 lineages in population.
 */
int msc_get_coal_pair_in_pop(gsl_rng* r, chrsample* sample,
                              int pop_id, coalescent_pair* pair);

/*
 * Process a divergence event: merge two daughter populations into parent.
 * All chromosomes in child1 and child2 get their population_id changed
 * to parent->id.
 */
void msc_process_divergence(chrsample* sample, divergence_event* event);

/*
 * Run MSC simulation.
 * Returns 0 on success, -1 on error.
 */
int msc_simulate(msc_params* params);

/* ============== Utility Functions ============== */

/*
 * Count active lineages per population.
 * Updates active_lineages field in each species_node.
 */
void msc_count_lineages(chrsample* sample, species_tree* tree);

/*
 * Print simulation summary.
 */
void msc_print_summary(msc_params* params, FILE* out);

/* ============== Migration Functions ============== */

/*
 * Select which migration pair has the migration event.
 * Returns index into rates->mig_rates, weighted by each pair's migration rate.
 */
int msc_select_migration_pair(gsl_rng* r, msc_rates* rates);

/*
 * Select a random lineage from a specific population.
 * Returns the chromosome index, or -1 if no lineages in population.
 */
int msc_select_lineage_in_pop(gsl_rng* r, chrsample* sample, int pop_id);

/*
 * Check if a population exists at the given time.
 * A population exists if: (1) it's a tip, or (2) current_time >= its divergence_time.
 * Additionally, it must not have merged yet (current_time < parent's divergence_time).
 */
int msc_population_exists_at_time(species_node* node, double current_time);

/*
 * Free migration band list.
 */
void msc_free_migrations(migration_band_list* migrations);

#endif /* MSC_H */
