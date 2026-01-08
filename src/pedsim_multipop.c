/*
 * pedsim_multipop.c - Multi-population Backwards Pedigree Simulator
 *
 * Simulates pedigrees backwards in time from samples in J populations.
 * Each population j has its own size N_j and sample size n_j.
 * Migration: probability m_{ij} that individual in pop j has parents from pop i.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "pedsim_multipop.h"

/* ============================================================================
 * INTERNAL HELPER FUNCTIONS
 * ============================================================================ */

/*
 * Create a new generation structure
 */
static multipop_generation* create_generation(int initial_capacity) {
    multipop_generation* gen = malloc(sizeof(multipop_generation));
    if (!gen) return NULL;

    gen->male_ids = malloc(initial_capacity * sizeof(int));
    gen->female_ids = malloc(initial_capacity * sizeof(int));
    if (!gen->male_ids || !gen->female_ids) {
        free(gen->male_ids);
        free(gen->female_ids);
        free(gen);
        return NULL;
    }

    gen->n_males = 0;
    gen->n_females = 0;
    gen->male_capacity = initial_capacity;
    gen->female_capacity = initial_capacity;

    return gen;
}

/*
 * Free a generation structure
 */
static void free_generation(multipop_generation* gen) {
    if (!gen) return;
    free(gen->male_ids);
    free(gen->female_ids);
    free(gen);
}

/*
 * Create gen_by_pop structure for one generation
 */
static multipop_gen_by_pop* create_gen_by_pop(int n_pops, int initial_capacity) {
    multipop_gen_by_pop* gbp = malloc(sizeof(multipop_gen_by_pop));
    if (!gbp) return NULL;

    gbp->generations = malloc(n_pops * sizeof(multipop_generation*));
    if (!gbp->generations) {
        free(gbp);
        return NULL;
    }

    gbp->n_pops = n_pops;
    for (int i = 0; i < n_pops; i++) {
        gbp->generations[i] = create_generation(initial_capacity);
        if (!gbp->generations[i]) {
            for (int j = 0; j < i; j++) {
                free_generation(gbp->generations[j]);
            }
            free(gbp->generations);
            free(gbp);
            return NULL;
        }
    }

    return gbp;
}

/*
 * Free gen_by_pop structure
 */
static void free_gen_by_pop(multipop_gen_by_pop* gbp) {
    if (!gbp) return;
    for (int i = 0; i < gbp->n_pops; i++) {
        free_generation(gbp->generations[i]);
    }
    free(gbp->generations);
    free(gbp);
}

/*
 * Add individual to the pedigree array, expanding if needed
 */
static int add_individual(multipop_pedigree* ped, int generation, int sex,
                          int father_id, int mother_id, int population_id,
                          int is_founder, int is_migrant) {
    /* Expand array if needed */
    if (ped->n_individuals >= ped->capacity) {
        int new_capacity = ped->capacity * 2;
        multipop_indiv* new_arr = realloc(ped->individuals,
                                          new_capacity * sizeof(multipop_indiv));
        if (!new_arr) return -1;
        ped->individuals = new_arr;
        ped->capacity = new_capacity;
    }

    int id = ped->n_individuals;
    multipop_indiv* indiv = &ped->individuals[id];
    indiv->id = id;
    indiv->generation = generation;
    indiv->sex = sex;
    indiv->father_id = father_id;
    indiv->mother_id = mother_id;
    indiv->population_id = population_id;
    indiv->is_founder = is_founder;
    indiv->is_migrant = is_migrant;

    ped->n_individuals++;
    return id;
}

/*
 * Add ID to generation's male or female list, expanding if needed
 */
static int add_to_generation(multipop_generation* gen, int id, int sex) {
    if (sex == 0) {
        /* Male */
        if (gen->n_males >= gen->male_capacity) {
            int new_cap = gen->male_capacity * 2;
            int* new_arr = realloc(gen->male_ids, new_cap * sizeof(int));
            if (!new_arr) return -1;
            gen->male_ids = new_arr;
            gen->male_capacity = new_cap;
        }
        gen->male_ids[gen->n_males++] = id;
    } else {
        /* Female */
        if (gen->n_females >= gen->female_capacity) {
            int new_cap = gen->female_capacity * 2;
            int* new_arr = realloc(gen->female_ids, new_cap * sizeof(int));
            if (!new_arr) return -1;
            gen->female_ids = new_arr;
            gen->female_capacity = new_cap;
        }
        gen->female_ids[gen->n_females++] = id;
    }
    return 0;
}

/*
 * Sample source population for parents given migration matrix
 * Returns population ID where parents should be drawn from
 */
static int sample_parent_population(gsl_rng* rng, double** migration,
                                     int n_pops, int child_pop) {
    double u = gsl_rng_uniform(rng);
    double cumsum = 0.0;

    /* Try each other population */
    for (int i = 0; i < n_pops; i++) {
        if (i == child_pop) continue;
        cumsum += migration[i][child_pop];
        if (u < cumsum) {
            return i;
        }
    }

    /* No migration - stay in own population */
    return child_pop;
}

/* ============================================================================
 * PEDIGREE CREATION AND DESTRUCTION
 * ============================================================================ */

double** multipop_island_migration(int n_pops, double m) {
    if (n_pops < 1 || m < 0.0 || m > 1.0) return NULL;

    double** mig = malloc(n_pops * sizeof(double*));
    if (!mig) return NULL;

    for (int i = 0; i < n_pops; i++) {
        mig[i] = malloc(n_pops * sizeof(double));
        if (!mig[i]) {
            for (int j = 0; j < i; j++) free(mig[j]);
            free(mig);
            return NULL;
        }
    }

    /* Set off-diagonal to m/(J-1), diagonal to 0 */
    double m_per_pop = (n_pops > 1) ? m / (n_pops - 1) : 0.0;
    for (int i = 0; i < n_pops; i++) {
        for (int j = 0; j < n_pops; j++) {
            mig[i][j] = (i == j) ? 0.0 : m_per_pop;
        }
    }

    return mig;
}

void multipop_free_migration(double** migration, int n_pops) {
    if (!migration) return;
    for (int i = 0; i < n_pops; i++) {
        free(migration[i]);
    }
    free(migration);
}

multipop_pedigree* multipop_create(multipop_pop_params* pop_params, int n_pops,
                                    double** migration, int k, unsigned long seed) {
    /* Validate parameters */
    if (!pop_params || n_pops <= 0 || n_pops > MULTIPOP_MAX_POPS || k <= 0) {
        fprintf(stderr, "Error: invalid parameters\n");
        return NULL;
    }

    for (int i = 0; i < n_pops; i++) {
        if (pop_params[i].pop_size <= 0 || pop_params[i].sample_size <= 0) {
            fprintf(stderr, "Error: population sizes must be positive\n");
            return NULL;
        }
        if (pop_params[i].pop_size % 2 != 0) {
            fprintf(stderr, "Error: population size must be even for dioecious model\n");
            return NULL;
        }
        if (pop_params[i].sample_size > pop_params[i].pop_size) {
            fprintf(stderr, "Error: sample size cannot exceed population size\n");
            return NULL;
        }
    }

    multipop_pedigree* ped = malloc(sizeof(multipop_pedigree));
    if (!ped) return NULL;

    /* Initialize individuals array */
    ped->individuals = malloc(MULTIPOP_INITIAL_CAPACITY * sizeof(multipop_indiv));
    if (!ped->individuals) {
        free(ped);
        return NULL;
    }
    ped->n_individuals = 0;
    ped->capacity = MULTIPOP_INITIAL_CAPACITY;

    /* Copy population parameters */
    ped->populations = malloc(n_pops * sizeof(multipop_pop_params));
    if (!ped->populations) {
        free(ped->individuals);
        free(ped);
        return NULL;
    }
    memcpy(ped->populations, pop_params, n_pops * sizeof(multipop_pop_params));
    ped->n_populations = n_pops;

    /* Copy migration matrix (or create identity if NULL) */
    ped->migration = malloc(n_pops * sizeof(double*));
    if (!ped->migration) {
        free(ped->populations);
        free(ped->individuals);
        free(ped);
        return NULL;
    }
    for (int i = 0; i < n_pops; i++) {
        ped->migration[i] = malloc(n_pops * sizeof(double));
        if (!ped->migration[i]) {
            for (int j = 0; j < i; j++) free(ped->migration[j]);
            free(ped->migration);
            free(ped->populations);
            free(ped->individuals);
            free(ped);
            return NULL;
        }
        if (migration) {
            memcpy(ped->migration[i], migration[i], n_pops * sizeof(double));
        } else {
            /* No migration */
            for (int j = 0; j < n_pops; j++) {
                ped->migration[i][j] = 0.0;
            }
        }
    }

    /* Initialize per-population generation tracking */
    ped->gens_by_pop = malloc((k + 1) * sizeof(multipop_gen_by_pop*));
    if (!ped->gens_by_pop) {
        multipop_free_migration(ped->migration, n_pops);
        free(ped->populations);
        free(ped->individuals);
        free(ped);
        return NULL;
    }
    for (int t = 0; t <= k; t++) {
        ped->gens_by_pop[t] = create_gen_by_pop(n_pops, MULTIPOP_INITIAL_CAPACITY);
        if (!ped->gens_by_pop[t]) {
            for (int j = 0; j < t; j++) {
                free_gen_by_pop(ped->gens_by_pop[j]);
            }
            free(ped->gens_by_pop);
            multipop_free_migration(ped->migration, n_pops);
            free(ped->populations);
            free(ped->individuals);
            free(ped);
            return NULL;
        }
    }
    ped->n_generations = k + 1;
    ped->k_generations = k;

    /* Initialize RNG */
    ped->rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!ped->rng) {
        for (int t = 0; t <= k; t++) {
            free_gen_by_pop(ped->gens_by_pop[t]);
        }
        free(ped->gens_by_pop);
        multipop_free_migration(ped->migration, n_pops);
        free(ped->populations);
        free(ped->individuals);
        free(ped);
        return NULL;
    }

    if (seed == 0) {
        seed = (unsigned long)time(NULL);
    }
    ped->seed = seed;
    gsl_rng_set(ped->rng, seed);

    return ped;
}

void multipop_free(multipop_pedigree* ped) {
    if (!ped) return;

    free(ped->individuals);

    for (int t = 0; t < ped->n_generations; t++) {
        free_gen_by_pop(ped->gens_by_pop[t]);
    }
    free(ped->gens_by_pop);

    multipop_free_migration(ped->migration, ped->n_populations);
    free(ped->populations);

    if (ped->rng) {
        gsl_rng_free(ped->rng);
    }

    free(ped);
}

void multipop_free_stats(multipop_stats* stats) {
    if (!stats) return;
    for (int p = 0; p < stats->n_populations; p++) {
        free(stats->n_males_by_gen[p]);
        free(stats->n_females_by_gen[p]);
        free(stats->n_total_by_gen[p]);
    }
    free(stats->n_males_by_gen);
    free(stats->n_females_by_gen);
    free(stats->n_total_by_gen);
    free(stats->n_migrants_by_gen);
    free(stats);
}

/* ============================================================================
 * SIMULATION
 * ============================================================================ */

int multipop_simulate(multipop_pedigree* ped) {
    if (!ped) return -1;

    int n_pops = ped->n_populations;

    /*
     * Algorithm:
     * 1. Create sampled individuals at generation 0 for each population
     * 2. For each generation t = 1 to k:
     *    For each population p:
     *      For each individual in population p at generation t-1:
     *        - Sample source population for parents (based on migration)
     *        - Draw father from N/2 males in source population
     *        - Draw mother from N/2 females in source population
     *    Collect unique ancestors per population
     * 3. Generation k ancestors become founders
     */

    /* Step 1: Create sampled individuals at generation 0 for each population */
    for (int p = 0; p < n_pops; p++) {
        for (int i = 0; i < ped->populations[p].sample_size; i++) {
            int sex = (gsl_rng_uniform_int(ped->rng, 2) == 0) ? 0 : 1;
            int id = add_individual(ped, 0, sex, -1, -1, p, 0, 0);
            if (id < 0) return -1;
            add_to_generation(ped->gens_by_pop[0]->generations[p], id, sex);
        }
    }

    /* Allocate parent maps for each population */
    int max_pop_size = 0;
    for (int p = 0; p < n_pops; p++) {
        if (ped->populations[p].pop_size > max_pop_size) {
            max_pop_size = ped->populations[p].pop_size;
        }
    }

    /* parent_map[pop][sex][pop_id] = pedigree individual ID */
    int*** parent_map = malloc(n_pops * sizeof(int**));
    if (!parent_map) return -1;
    for (int p = 0; p < n_pops; p++) {
        parent_map[p] = malloc(2 * sizeof(int*));
        if (!parent_map[p]) {
            for (int q = 0; q < p; q++) {
                free(parent_map[q][0]);
                free(parent_map[q][1]);
                free(parent_map[q]);
            }
            free(parent_map);
            return -1;
        }
        int n_each = ped->populations[p].pop_size / 2;
        parent_map[p][0] = malloc(n_each * sizeof(int));
        parent_map[p][1] = malloc(n_each * sizeof(int));
        if (!parent_map[p][0] || !parent_map[p][1]) {
            /* cleanup and return */
            free(parent_map[p][0]);
            free(parent_map[p][1]);
            free(parent_map[p]);
            for (int q = 0; q < p; q++) {
                free(parent_map[q][0]);
                free(parent_map[q][1]);
                free(parent_map[q]);
            }
            free(parent_map);
            return -1;
        }
    }

    /* Step 2: Iterate backwards through generations */
    for (int t = 1; t <= ped->k_generations; t++) {
        int is_founder = (t == ped->k_generations) ? 1 : 0;

        /* Reset parent maps for all populations */
        for (int p = 0; p < n_pops; p++) {
            int n_each = ped->populations[p].pop_size / 2;
            for (int i = 0; i < n_each; i++) {
                parent_map[p][0][i] = -1;
                parent_map[p][1][i] = -1;
            }
        }

        /* Process each population */
        for (int p = 0; p < n_pops; p++) {
            multipop_generation* prev_gen = ped->gens_by_pop[t - 1]->generations[p];
            int n_prev = prev_gen->n_males + prev_gen->n_females;
            if (n_prev == 0) continue;

            /* Process each individual at previous generation in this population */
            for (int i = 0; i < prev_gen->n_males; i++) {
                int child_id = prev_gen->male_ids[i];

                /* Sample source population for parents */
                int src_pop = sample_parent_population(ped->rng, ped->migration,
                                                        n_pops, p);
                int is_migrant = (src_pop != p) ? 1 : 0;
                ped->individuals[child_id].is_migrant = is_migrant;

                int n_males_src = ped->populations[src_pop].pop_size / 2;
                int n_females_src = ped->populations[src_pop].pop_size / 2;

                /* Draw father and mother from source population */
                int f_pop = gsl_rng_uniform_int(ped->rng, n_males_src);
                int m_pop = gsl_rng_uniform_int(ped->rng, n_females_src);

                /* Get or create father */
                if (parent_map[src_pop][0][f_pop] < 0) {
                    int father_id = add_individual(ped, t, 0, -1, -1, src_pop,
                                                    is_founder, 0);
                    if (father_id < 0) goto cleanup_error;
                    parent_map[src_pop][0][f_pop] = father_id;
                    add_to_generation(ped->gens_by_pop[t]->generations[src_pop],
                                      father_id, 0);
                }
                ped->individuals[child_id].father_id = parent_map[src_pop][0][f_pop];

                /* Get or create mother */
                if (parent_map[src_pop][1][m_pop] < 0) {
                    int mother_id = add_individual(ped, t, 1, -1, -1, src_pop,
                                                    is_founder, 0);
                    if (mother_id < 0) goto cleanup_error;
                    parent_map[src_pop][1][m_pop] = mother_id;
                    add_to_generation(ped->gens_by_pop[t]->generations[src_pop],
                                      mother_id, 1);
                }
                ped->individuals[child_id].mother_id = parent_map[src_pop][1][m_pop];
            }

            for (int i = 0; i < prev_gen->n_females; i++) {
                int child_id = prev_gen->female_ids[i];

                /* Sample source population for parents */
                int src_pop = sample_parent_population(ped->rng, ped->migration,
                                                        n_pops, p);
                int is_migrant = (src_pop != p) ? 1 : 0;
                ped->individuals[child_id].is_migrant = is_migrant;

                int n_males_src = ped->populations[src_pop].pop_size / 2;
                int n_females_src = ped->populations[src_pop].pop_size / 2;

                /* Draw father and mother from source population */
                int f_pop = gsl_rng_uniform_int(ped->rng, n_males_src);
                int m_pop = gsl_rng_uniform_int(ped->rng, n_females_src);

                /* Get or create father */
                if (parent_map[src_pop][0][f_pop] < 0) {
                    int father_id = add_individual(ped, t, 0, -1, -1, src_pop,
                                                    is_founder, 0);
                    if (father_id < 0) goto cleanup_error;
                    parent_map[src_pop][0][f_pop] = father_id;
                    add_to_generation(ped->gens_by_pop[t]->generations[src_pop],
                                      father_id, 0);
                }
                ped->individuals[child_id].father_id = parent_map[src_pop][0][f_pop];

                /* Get or create mother */
                if (parent_map[src_pop][1][m_pop] < 0) {
                    int mother_id = add_individual(ped, t, 1, -1, -1, src_pop,
                                                    is_founder, 0);
                    if (mother_id < 0) goto cleanup_error;
                    parent_map[src_pop][1][m_pop] = mother_id;
                    add_to_generation(ped->gens_by_pop[t]->generations[src_pop],
                                      mother_id, 1);
                }
                ped->individuals[child_id].mother_id = parent_map[src_pop][1][m_pop];
            }
        }
    }

    /* Cleanup */
    for (int p = 0; p < n_pops; p++) {
        free(parent_map[p][0]);
        free(parent_map[p][1]);
        free(parent_map[p]);
    }
    free(parent_map);
    return 0;

cleanup_error:
    for (int p = 0; p < n_pops; p++) {
        free(parent_map[p][0]);
        free(parent_map[p][1]);
        free(parent_map[p]);
    }
    free(parent_map);
    return -1;
}

multipop_stats* multipop_get_stats(multipop_pedigree* ped) {
    if (!ped) return NULL;

    multipop_stats* stats = malloc(sizeof(multipop_stats));
    if (!stats) return NULL;

    stats->n_populations = ped->n_populations;
    stats->n_generations = ped->n_generations;

    stats->n_males_by_gen = malloc(ped->n_populations * sizeof(int*));
    stats->n_females_by_gen = malloc(ped->n_populations * sizeof(int*));
    stats->n_total_by_gen = malloc(ped->n_populations * sizeof(int*));
    stats->n_migrants_by_gen = calloc(ped->n_generations, sizeof(int));

    if (!stats->n_males_by_gen || !stats->n_females_by_gen ||
        !stats->n_total_by_gen || !stats->n_migrants_by_gen) {
        multipop_free_stats(stats);
        return NULL;
    }

    for (int p = 0; p < ped->n_populations; p++) {
        stats->n_males_by_gen[p] = malloc(ped->n_generations * sizeof(int));
        stats->n_females_by_gen[p] = malloc(ped->n_generations * sizeof(int));
        stats->n_total_by_gen[p] = malloc(ped->n_generations * sizeof(int));
        if (!stats->n_males_by_gen[p] || !stats->n_females_by_gen[p] ||
            !stats->n_total_by_gen[p]) {
            multipop_free_stats(stats);
            return NULL;
        }
    }

    for (int t = 0; t < ped->n_generations; t++) {
        for (int p = 0; p < ped->n_populations; p++) {
            multipop_generation* gen = ped->gens_by_pop[t]->generations[p];
            stats->n_males_by_gen[p][t] = gen->n_males;
            stats->n_females_by_gen[p][t] = gen->n_females;
            stats->n_total_by_gen[p][t] = gen->n_males + gen->n_females;
        }
    }

    /* Count migrants */
    for (int i = 0; i < ped->n_individuals; i++) {
        if (ped->individuals[i].is_migrant) {
            stats->n_migrants_by_gen[ped->individuals[i].generation]++;
        }
    }

    return stats;
}

/* ============================================================================
 * OUTPUT FUNCTIONS
 * ============================================================================ */

void multipop_write_pedigree(multipop_pedigree* ped, FILE* out) {
    if (!ped || !out) return;

    fprintf(out, "# Multipopulation pedigree generated by pedsim_multipop\n");
    fprintf(out, "# Parameters: J=%d populations, k=%d generations, seed=%lu\n",
            ped->n_populations, ped->k_generations, ped->seed);

    /* Write population info */
    for (int p = 0; p < ped->n_populations; p++) {
        fprintf(out, "# Pop%d: N=%d, n=%d\n",
                p, ped->populations[p].pop_size, ped->populations[p].sample_size);
    }

    /* Write migration matrix */
    fprintf(out, "# Migration matrix:\n");
    for (int i = 0; i < ped->n_populations; i++) {
        fprintf(out, "#  ");
        for (int j = 0; j < ped->n_populations; j++) {
            fprintf(out, " %.4f", ped->migration[i][j]);
        }
        fprintf(out, "\n");
    }

    fprintf(out, "# Total individuals: %d\n", ped->n_individuals);
    fprintf(out, "# Format: individual father mother population\n");
    fprintf(out, "#\n");

    /*
     * Output in topological order (founders first, then descendants)
     * Go from generation k down to generation 0, all populations at each gen
     */
    for (int t = ped->k_generations; t >= 0; t--) {
        for (int p = 0; p < ped->n_populations; p++) {
            multipop_generation* gen = ped->gens_by_pop[t]->generations[p];

            /* Output males at this generation in this population */
            for (int i = 0; i < gen->n_males; i++) {
                int id = gen->male_ids[i];
                multipop_indiv* indiv = &ped->individuals[id];

                if (indiv->is_founder) {
                    fprintf(out, "P%d_I%d 0 0 %d\n", p, id, p);
                } else {
                    /* Find father's and mother's populations for naming */
                    int f_pop = ped->individuals[indiv->father_id].population_id;
                    int m_pop = ped->individuals[indiv->mother_id].population_id;
                    fprintf(out, "P%d_I%d P%d_I%d P%d_I%d %d\n",
                            p, id, f_pop, indiv->father_id,
                            m_pop, indiv->mother_id, p);
                }
            }

            /* Output females at this generation in this population */
            for (int i = 0; i < gen->n_females; i++) {
                int id = gen->female_ids[i];
                multipop_indiv* indiv = &ped->individuals[id];

                if (indiv->is_founder) {
                    fprintf(out, "P%d_I%d 0 0 %d\n", p, id, p);
                } else {
                    int f_pop = ped->individuals[indiv->father_id].population_id;
                    int m_pop = ped->individuals[indiv->mother_id].population_id;
                    fprintf(out, "P%d_I%d P%d_I%d P%d_I%d %d\n",
                            p, id, f_pop, indiv->father_id,
                            m_pop, indiv->mother_id, p);
                }
            }
        }
    }
}

void multipop_write_stats(multipop_stats* stats, FILE* out) {
    if (!stats || !out) return;

    fprintf(out, "# Generation statistics by population\n");
    fprintf(out, "# Gen  ");
    for (int p = 0; p < stats->n_populations; p++) {
        fprintf(out, "Pop%d(M/F/T)      ", p);
    }
    fprintf(out, "Migrants\n");

    for (int t = 0; t < stats->n_generations; t++) {
        fprintf(out, "  %3d  ", t);
        for (int p = 0; p < stats->n_populations; p++) {
            fprintf(out, "%4d/%4d/%4d  ",
                    stats->n_males_by_gen[p][t],
                    stats->n_females_by_gen[p][t],
                    stats->n_total_by_gen[p][t]);
        }
        fprintf(out, "%4d\n", stats->n_migrants_by_gen[t]);
    }
}

void multipop_write_dot(multipop_pedigree* ped, FILE* out) {
    if (!ped || !out) return;

    /* Colors for different populations */
    const char* pop_colors[] = {
        "lightblue", "lightgreen", "lightyellow", "lightpink",
        "lavender", "peachpuff", "lightcyan", "mistyrose"
    };
    int n_colors = 8;

    fprintf(out, "digraph pedigree {\n");
    fprintf(out, "  rankdir=TB;\n");
    fprintf(out, "  node [shape=box];\n");
    fprintf(out, "\n");

    /* Define nodes with population-based colors */
    for (int i = 0; i < ped->n_individuals; i++) {
        multipop_indiv* indiv = &ped->individuals[i];
        const char* color = pop_colors[indiv->population_id % n_colors];
        const char* shape = indiv->is_founder ? "ellipse" : "box";
        fprintf(out, "  P%d_I%d [label=\"P%d_I%d\\n(gen %d)\", fillcolor=%s, "
                "style=filled, shape=%s];\n",
                indiv->population_id, i, indiv->population_id, i,
                indiv->generation, color, shape);
    }
    fprintf(out, "\n");

    /* Define edges (parent -> child) */
    for (int i = 0; i < ped->n_individuals; i++) {
        multipop_indiv* indiv = &ped->individuals[i];
        if (indiv->father_id >= 0) {
            int f_pop = ped->individuals[indiv->father_id].population_id;
            fprintf(out, "  P%d_I%d -> P%d_I%d;\n",
                    f_pop, indiv->father_id, indiv->population_id, i);
        }
        if (indiv->mother_id >= 0) {
            int m_pop = ped->individuals[indiv->mother_id].population_id;
            fprintf(out, "  P%d_I%d -> P%d_I%d;\n",
                    m_pop, indiv->mother_id, indiv->population_id, i);
        }
    }

    fprintf(out, "}\n");
}

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

int* multipop_parse_int_list(const char* str, int* count) {
    if (!str || !count) return NULL;

    /* Count commas to determine array size */
    int n = 1;
    for (const char* p = str; *p; p++) {
        if (*p == ',') n++;
    }

    int* arr = malloc(n * sizeof(int));
    if (!arr) return NULL;

    char* copy = strdup(str);
    if (!copy) {
        free(arr);
        return NULL;
    }

    char* token = strtok(copy, ",");
    int i = 0;
    while (token && i < n) {
        arr[i++] = atoi(token);
        token = strtok(NULL, ",");
    }

    free(copy);
    *count = i;
    return arr;
}

double** multipop_parse_migration(const char* str, int n_pops) {
    if (!str || n_pops <= 0) return NULL;

    double** mig = malloc(n_pops * sizeof(double*));
    if (!mig) return NULL;

    for (int i = 0; i < n_pops; i++) {
        mig[i] = malloc(n_pops * sizeof(double));
        if (!mig[i]) {
            for (int j = 0; j < i; j++) free(mig[j]);
            free(mig);
            return NULL;
        }
        /* Initialize to 0 */
        for (int j = 0; j < n_pops; j++) {
            mig[i][j] = 0.0;
        }
    }

    char* copy = strdup(str);
    if (!copy) {
        multipop_free_migration(mig, n_pops);
        return NULL;
    }

    /* Parse rows separated by semicolons */
    char* row_str = strtok(copy, ";");
    int row = 0;
    while (row_str && row < n_pops) {
        /* Parse columns separated by commas */
        char* row_copy = strdup(row_str);
        char* col_str = strtok(row_copy, ",");
        int col = 0;
        while (col_str && col < n_pops) {
            mig[row][col] = atof(col_str);
            col_str = strtok(NULL, ",");
            col++;
        }
        free(row_copy);
        row_str = strtok(NULL, ";");
        row++;
    }

    free(copy);
    return mig;
}
