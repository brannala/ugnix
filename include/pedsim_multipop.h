/*
 * pedsim_multipop.h - Multi-population Backwards Pedigree Simulator
 *
 * Simulates pedigrees backwards in time from samples in J populations.
 * Each population j has its own size N_j and sample size n_j.
 * Migration: probability m_{ij} that individual in pop j has parents from pop i.
 *
 * Based on the Wright-Fisher model with N/2 males and N/2 females per population.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#ifndef PEDSIM_MULTIPOP_H
#define PEDSIM_MULTIPOP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>

/* Maximum length for individual names */
#define MULTIPOP_MAX_NAME 64

/* Maximum number of populations */
#define MULTIPOP_MAX_POPS 100

/* Initial capacity for dynamic arrays */
#define MULTIPOP_INITIAL_CAPACITY 1024

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================ */

/*
 * Parameters for a single population
 */
typedef struct {
    int pop_id;             /* Population identifier (0, 1, ..., J-1) */
    int pop_size;           /* N_j: population size (must be even) */
    int sample_size;        /* n_j: number of sampled individuals */
} multipop_pop_params;

/*
 * Individual in the simulated multi-population pedigree
 */
typedef struct {
    int id;                 /* Unique identifier (global across all populations) */
    int generation;         /* Generation (0 = sample, positive = past) */
    int sex;                /* 0 = male, 1 = female */
    int father_id;          /* Father's ID (-1 if founder) */
    int mother_id;          /* Mother's ID (-1 if founder) */
    int population_id;      /* Which population this individual belongs to */
    int is_founder;         /* 1 if founder (no parents simulated), 0 otherwise */
    int is_migrant;         /* 1 if parents drawn from different population */
} multipop_indiv;

/*
 * Container for individuals at a single generation in one population
 */
typedef struct {
    int* male_ids;          /* IDs of males at this generation */
    int* female_ids;        /* IDs of females at this generation */
    int n_males;            /* Number of males */
    int n_females;          /* Number of females */
    int male_capacity;      /* Allocated capacity for males */
    int female_capacity;    /* Allocated capacity for females */
} multipop_generation;

/*
 * Per-population generation tracking
 */
typedef struct {
    multipop_generation** generations;  /* generations[pop_id] = generation struct */
    int n_pops;
} multipop_gen_by_pop;

/*
 * Complete simulated multi-population pedigree
 */
typedef struct {
    /* All individuals across all populations */
    multipop_indiv* individuals;
    int n_individuals;
    int capacity;

    /* Per-population generation tracking */
    multipop_gen_by_pop** gens_by_pop;  /* gens_by_pop[gen][pop] */
    int n_generations;

    /* Population parameters */
    multipop_pop_params* populations;
    int n_populations;

    /* Migration matrix: migration[i][j] = prob individual in j has parents from i */
    double** migration;

    /* Shared parameters */
    int k_generations;          /* k: number of generations to simulate */

    /* Random number generator */
    gsl_rng* rng;
    unsigned long seed;
} multipop_pedigree;

/*
 * Statistics about the simulation
 */
typedef struct {
    int** n_males_by_gen;       /* n_males_by_gen[pop][gen] */
    int** n_females_by_gen;     /* n_females_by_gen[pop][gen] */
    int** n_total_by_gen;       /* n_total_by_gen[pop][gen] */
    int* n_migrants_by_gen;     /* Total migrants per generation */
    int n_populations;
    int n_generations;
} multipop_stats;

/* ============================================================================
 * PEDIGREE CREATION AND DESTRUCTION
 * ============================================================================ */

/*
 * Create a new multi-population pedigree simulation
 * Parameters:
 *   pop_params: array of population parameters (n_pops elements)
 *   n_pops: number of populations
 *   migration: J x J migration matrix (migration[i][j] = prob parents from i for ind in j)
 *              Diagonal entries are ignored (computed as 1 - sum of off-diagonal)
 *   k: number of generations to simulate backwards
 *   seed: RNG seed (0 for time-based)
 * Returns NULL on error
 */
multipop_pedigree* multipop_create(multipop_pop_params* pop_params, int n_pops,
                                    double** migration, int k, unsigned long seed);

/*
 * Create migration matrix for symmetric island model
 * All off-diagonal entries set to m/(J-1)
 * Returns allocated J x J matrix (caller must free with multipop_free_migration)
 */
double** multipop_island_migration(int n_pops, double m);

/*
 * Free migration matrix
 */
void multipop_free_migration(double** migration, int n_pops);

/*
 * Free pedigree and all associated memory
 */
void multipop_free(multipop_pedigree* ped);

/*
 * Free statistics structure
 */
void multipop_free_stats(multipop_stats* stats);

/* ============================================================================
 * SIMULATION
 * ============================================================================ */

/*
 * Run the backwards pedigree simulation
 * Populates the pedigree structure with individuals and parent relationships
 * Returns 0 on success, -1 on error
 */
int multipop_simulate(multipop_pedigree* ped);

/*
 * Get simulation statistics
 * Returns allocated structure (caller must free with multipop_free_stats)
 */
multipop_stats* multipop_get_stats(multipop_pedigree* ped);

/* ============================================================================
 * OUTPUT FUNCTIONS
 * ============================================================================ */

/*
 * Write pedigree to file in extended pedtrans-compatible format
 * Format: "indiv_name father_name mother_name pop_id" per line
 * Founders have father = mother = "0"
 * Individual names are prefixed with population: P0_I0, P1_I0, etc.
 */
void multipop_write_pedigree(multipop_pedigree* ped, FILE* out);

/*
 * Write statistics summary to file
 */
void multipop_write_stats(multipop_stats* stats, FILE* out);

/*
 * Write pedigree in DOT format for visualization
 */
void multipop_write_dot(multipop_pedigree* ped, FILE* out);

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

/*
 * Parse comma-separated list of integers
 * Returns allocated array (caller must free)
 * Sets *count to number of elements
 */
int* multipop_parse_int_list(const char* str, int* count);

/*
 * Parse migration matrix string
 * Format: "m00,m01,m02;m10,m11,m12;m20,m21,m22"
 * Returns allocated matrix (caller must free with multipop_free_migration)
 */
double** multipop_parse_migration(const char* str, int n_pops);

#endif /* PEDSIM_MULTIPOP_H */
