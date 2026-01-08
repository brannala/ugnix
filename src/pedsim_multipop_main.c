/*
 * pedsim_multipop_main.c - CLI for multi-population pedigree simulator
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "pedsim_multipop.h"

static void print_usage(const char* progname) {
    fprintf(stderr, "Usage: %s [options] -o <output.ped>\n", progname);
    fprintf(stderr, "\nMulti-population backwards pedigree simulation\n\n");
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -o FILE   Output pedigree file\n");
    fprintf(stderr, "\nPopulation options:\n");
    fprintf(stderr, "  -J NUM    Number of populations (default: 1)\n");
    fprintf(stderr, "  -N LIST   Population sizes, comma-separated (e.g., 1000,500,2000)\n");
    fprintf(stderr, "            If single value, used for all populations\n");
    fprintf(stderr, "  -n LIST   Sample sizes, comma-separated (e.g., 10,5,20)\n");
    fprintf(stderr, "            If single value, used for all populations\n");
    fprintf(stderr, "  -k NUM    Number of generations to trace back (default: 5)\n");
    fprintf(stderr, "\nMigration options:\n");
    fprintf(stderr, "  -m STR    Full migration matrix, semicolon-separated rows\n");
    fprintf(stderr, "            e.g., \"0,0.01;0.01,0\" for 2 populations\n");
    fprintf(stderr, "            m[i][j] = prob individual in j has parents from i\n");
    fprintf(stderr, "  -M NUM    Symmetric island model migration rate (alternative to -m)\n");
    fprintf(stderr, "            Total migration rate m, divided equally among other pops\n");
    fprintf(stderr, "\nOther options:\n");
    fprintf(stderr, "  -s NUM    Random seed (default: 0 = time-based)\n");
    fprintf(stderr, "  -v        Verbose output (print statistics)\n");
    fprintf(stderr, "  -d        Output DOT format for visualization\n");
    fprintf(stderr, "  -h        Show this help message\n");
    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  # Single population (like original pedsim)\n");
    fprintf(stderr, "  %s -J 1 -N 1000 -n 10 -k 5 -o ped.txt\n", progname);
    fprintf(stderr, "\n  # Two populations with island model migration\n");
    fprintf(stderr, "  %s -J 2 -N 1000,500 -n 10,5 -M 0.01 -o ped.txt\n", progname);
    fprintf(stderr, "\n  # Three populations with custom migration matrix\n");
    fprintf(stderr, "  %s -J 3 -N 1000,500,2000 -n 10,5,20 \\\n", progname);
    fprintf(stderr, "     -m \"0,0.01,0.02;0.01,0,0.01;0.02,0.01,0\" -o ped.txt\n");
}

int main(int argc, char** argv) {
    int n_pops = 1;
    char* pop_sizes_str = NULL;
    char* sample_sizes_str = NULL;
    int k_generations = 5;
    char* migration_str = NULL;
    double island_m = -1.0;  /* -1 means not set */
    unsigned long seed = 0;
    char* output_file = NULL;
    int verbose = 0;
    int dot_output = 0;

    int c;
    while ((c = getopt(argc, argv, "J:N:n:k:m:M:s:o:vdh")) != -1) {
        switch (c) {
            case 'J':
                n_pops = atoi(optarg);
                break;
            case 'N':
                pop_sizes_str = optarg;
                break;
            case 'n':
                sample_sizes_str = optarg;
                break;
            case 'k':
                k_generations = atoi(optarg);
                break;
            case 'm':
                migration_str = optarg;
                break;
            case 'M':
                island_m = atof(optarg);
                break;
            case 's':
                seed = atol(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'd':
                dot_output = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    /* Validate required options */
    if (!output_file) {
        fprintf(stderr, "Error: output file required (-o)\n");
        print_usage(argv[0]);
        return 1;
    }

    if (n_pops < 1 || n_pops > MULTIPOP_MAX_POPS) {
        fprintf(stderr, "Error: number of populations must be 1-%d\n",
                MULTIPOP_MAX_POPS);
        return 1;
    }

    if (migration_str && island_m >= 0) {
        fprintf(stderr, "Error: cannot specify both -m and -M\n");
        return 1;
    }

    /* Parse population sizes */
    int* pop_sizes = NULL;
    int n_pop_sizes = 0;
    if (pop_sizes_str) {
        pop_sizes = multipop_parse_int_list(pop_sizes_str, &n_pop_sizes);
        if (!pop_sizes) {
            fprintf(stderr, "Error: failed to parse population sizes\n");
            return 1;
        }
    }

    /* Parse sample sizes */
    int* sample_sizes = NULL;
    int n_sample_sizes = 0;
    if (sample_sizes_str) {
        sample_sizes = multipop_parse_int_list(sample_sizes_str, &n_sample_sizes);
        if (!sample_sizes) {
            fprintf(stderr, "Error: failed to parse sample sizes\n");
            free(pop_sizes);
            return 1;
        }
    }

    /* Build population parameters */
    multipop_pop_params* params = malloc(n_pops * sizeof(multipop_pop_params));
    if (!params) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(pop_sizes);
        free(sample_sizes);
        return 1;
    }

    for (int i = 0; i < n_pops; i++) {
        params[i].pop_id = i;

        /* Population size */
        if (pop_sizes) {
            if (n_pop_sizes == 1) {
                params[i].pop_size = pop_sizes[0];
            } else if (i < n_pop_sizes) {
                params[i].pop_size = pop_sizes[i];
            } else {
                fprintf(stderr, "Error: not enough population sizes specified\n");
                free(params);
                free(pop_sizes);
                free(sample_sizes);
                return 1;
            }
        } else {
            params[i].pop_size = 1000;  /* Default */
        }

        /* Sample size */
        if (sample_sizes) {
            if (n_sample_sizes == 1) {
                params[i].sample_size = sample_sizes[0];
            } else if (i < n_sample_sizes) {
                params[i].sample_size = sample_sizes[i];
            } else {
                fprintf(stderr, "Error: not enough sample sizes specified\n");
                free(params);
                free(pop_sizes);
                free(sample_sizes);
                return 1;
            }
        } else {
            params[i].sample_size = 10;  /* Default */
        }
    }

    free(pop_sizes);
    free(sample_sizes);

    /* Build migration matrix */
    double** migration = NULL;
    if (migration_str) {
        migration = multipop_parse_migration(migration_str, n_pops);
        if (!migration) {
            fprintf(stderr, "Error: failed to parse migration matrix\n");
            free(params);
            return 1;
        }
    } else if (island_m >= 0) {
        migration = multipop_island_migration(n_pops, island_m);
        if (!migration) {
            fprintf(stderr, "Error: failed to create island migration matrix\n");
            free(params);
            return 1;
        }
    }
    /* If no migration specified, multipop_create will handle NULL */

    /* Print parameters */
    fprintf(stderr, "Multi-population Pedigree Simulation\n");
    fprintf(stderr, "====================================\n");
    fprintf(stderr, "Populations: %d\n", n_pops);
    for (int i = 0; i < n_pops; i++) {
        fprintf(stderr, "  Pop%d: N=%d, n=%d\n", i,
                params[i].pop_size, params[i].sample_size);
    }
    fprintf(stderr, "Generations: %d\n", k_generations);
    if (migration) {
        fprintf(stderr, "Migration matrix:\n");
        for (int i = 0; i < n_pops; i++) {
            fprintf(stderr, "  ");
            for (int j = 0; j < n_pops; j++) {
                fprintf(stderr, "%.4f ", migration[i][j]);
            }
            fprintf(stderr, "\n");
        }
    } else {
        fprintf(stderr, "Migration: none\n");
    }

    /* Create and run simulation */
    multipop_pedigree* ped = multipop_create(params, n_pops, migration, k_generations, seed);
    if (!ped) {
        fprintf(stderr, "Error: failed to create pedigree\n");
        if (migration) multipop_free_migration(migration, n_pops);
        free(params);
        return 1;
    }

    if (migration) multipop_free_migration(migration, n_pops);
    free(params);

    fprintf(stderr, "\nSimulating pedigree...\n");
    if (multipop_simulate(ped) < 0) {
        fprintf(stderr, "Error: simulation failed\n");
        multipop_free(ped);
        return 1;
    }

    fprintf(stderr, "Total individuals: %d\n", ped->n_individuals);

    /* Print statistics if verbose */
    if (verbose) {
        multipop_stats* stats = multipop_get_stats(ped);
        if (stats) {
            fprintf(stderr, "\n");
            multipop_write_stats(stats, stderr);
            multipop_free_stats(stats);
        }
    }

    /* Write output */
    FILE* out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "Error: cannot open output file: %s\n", output_file);
        multipop_free(ped);
        return 1;
    }

    if (dot_output) {
        multipop_write_dot(ped, out);
    } else {
        multipop_write_pedigree(ped, out);
    }
    fclose(out);

    fprintf(stderr, "Output written to: %s\n", output_file);

    multipop_free(ped);
    return 0;
}
