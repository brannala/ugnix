#ifndef SPECIES_TREE_H
#define SPECIES_TREE_H

#include <stdio.h>

/* Maximum species name length */
#define MSC_MAX_SPECIES_NAME 64

/* Maximum number of species/populations */
#define MSC_MAX_SPECIES 100

/*
 * Species tree node - represents a species/population
 * Internal nodes represent ancestral populations, tips represent sampled species
 */
typedef struct species_node {
    int id;                              /* Unique node ID (0..n_nodes-1) */
    char name[MSC_MAX_SPECIES_NAME];     /* Species name (tips only, empty for internal) */
    double branch_length;                /* Branch length to parent (in generations) */
    double theta;                        /* 4Nu for this population */
    double divergence_time;              /* Absolute time from present (computed) */

    struct species_node* parent;         /* Parent node (NULL for root) */
    struct species_node* left;           /* Left child (NULL for tips) */
    struct species_node* right;          /* Right child (NULL for tips) */

    int is_tip;                          /* 1 if tip (sampled species), 0 if internal */
    int sample_size;                     /* Number of sampled chromosomes (tips only) */
    int active_lineages;                 /* Current lineage count during simulation */
} species_node;

/*
 * Complete species tree structure
 */
typedef struct {
    species_node* root;                  /* Root of tree */
    species_node** nodes;                /* Array of all nodes for O(1) access */
    int n_nodes;                         /* Total nodes (tips + internal) */
    int n_tips;                          /* Number of tip species */
    species_node** tips;                 /* Array of tip nodes for quick iteration */
} species_tree;

/*
 * Divergence event for event-driven simulation
 * Events are processed going backwards in time
 */
typedef struct divergence_event {
    double time;                         /* Time of divergence (going backwards) */
    species_node* child1;                /* First daughter population */
    species_node* child2;                /* Second daughter population */
    species_node* parent;                /* Ancestral population to merge into */
    struct divergence_event* next;       /* Next event in time-sorted list */
} divergence_event;

/*
 * List of pending divergence events (sorted by time ascending)
 */
typedef struct {
    divergence_event* head;
    int n_events;
} divergence_event_list;

/* ============== Species Tree Parsing ============== */

/*
 * Parse Newick string with theta annotations.
 *
 * Extended Newick format: name:branch_length#theta
 * Example: ((A:1000#0.001,B:1000#0.002):500#0.003,C:1500#0.001)#0.002;
 *
 * Returns species_tree on success, NULL on error.
 */
species_tree* species_tree_parse_newick(const char* newick);

/*
 * Compute absolute divergence times from branch lengths.
 * Must be called after parsing before simulation.
 *
 * For each internal node, divergence_time is when (going backwards)
 * the daughter populations merge into the ancestral population.
 */
void species_tree_compute_divergence_times(species_tree* tree);

/*
 * Create sorted list of divergence events for simulation.
 * Events are sorted by time (ascending, as we go backwards).
 */
divergence_event_list* species_tree_create_divergence_events(species_tree* tree);

/* ============== Species Tree Operations ============== */

/*
 * Find species node by name, returns NULL if not found.
 * Only searches tip nodes.
 */
species_node* species_tree_find_species(species_tree* tree, const char* name);

/*
 * Find any node (tip or internal) by name, returns NULL if not found.
 * Searches all nodes in the tree.
 */
species_node* species_tree_find_node_by_name(species_tree* tree, const char* name);

/*
 * Get total sample size across all species.
 */
int species_tree_total_samples(species_tree* tree);

/*
 * Print species tree for debugging.
 */
void species_tree_print(species_tree* tree, FILE* out);

/*
 * Print tree in Newick format (without theta annotations).
 */
void species_tree_print_newick(species_tree* tree, FILE* out);

/* ============== Memory Management ============== */

/*
 * Free species tree and all nodes.
 */
void species_tree_free(species_tree* tree);

/*
 * Free divergence event list.
 */
void divergence_event_list_free(divergence_event_list* list);

#endif /* SPECIES_TREE_H */
