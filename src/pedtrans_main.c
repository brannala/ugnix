/*
 * pedtrans_main.c - CLI for Pedigree Chromosome Transmission Simulator
 *
 * Simulates chromosome transmission through a pedigree with recombination,
 * tracking founder chromosome segments.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include "pedtrans.h"

/* Default values */
#define DEFAULT_REC_RATE 1.0
#define DEFAULT_N_CHROMOSOMES 1

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [OPTIONS] PEDIGREE_FILE\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Simulate chromosome transmission through a pedigree with recombination.\n");
    fprintf(stderr, "Tracks founder chromosome origin for each genomic segment.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -c NUM     Number of chromosomes to simulate (default: %d)\n", DEFAULT_N_CHROMOSOMES);
    fprintf(stderr, "  -r RATE    Recombination rate (expected crossovers per meiosis)\n");
    fprintf(stderr, "             Default: %.1f\n", DEFAULT_REC_RATE);
    fprintf(stderr, "  -R LIST    Per-chromosome recombination rates, comma-separated\n");
    fprintf(stderr, "             (map lengths in Morgans); sets the number of chromosomes\n");
    fprintf(stderr, "             and overrides -c and -r\n");
    fprintf(stderr, "  -s SEED    Random seed for reproducibility\n");
    fprintf(stderr, "             Default: time-based\n");
    fprintf(stderr, "  -o FILE    Output file (default: stdout)\n");
    fprintf(stderr, "  -A FILE    Write the number of genealogical and genetic ancestors\n");
    fprintf(stderr, "             (founders) of each individual without children (\"-\" = stdout);\n");
    fprintf(stderr, "             segments are then written only if -o is given\n");
    fprintf(stderr, "  -a         With -A, summarize all individuals, not only those\n");
    fprintf(stderr, "             without children\n");
    fprintf(stderr, "  -v         Verbose output (print progress)\n");
    fprintf(stderr, "  -h         Show this help message\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input format (one line per individual):\n");
    fprintf(stderr, "  indivID fatherID motherID\n");
    fprintf(stderr, "  Founders have fatherID = motherID = 0\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output format:\n");
    fprintf(stderr, "  For each individual, lists chromosome segments with founder origin.\n");
    fprintf(stderr, "  Positions are normalized to [0, 1].\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "  %s -c 5 -r 1.5 -s 12345 pedigree.txt\n", program_name);
    fprintf(stderr, "  %s -r 2.0 -o output.txt pedigree.txt\n", program_name);
    fprintf(stderr, "  %s -R 2.86,2.68,2.23 -A ancestors.tsv pedigree.txt\n", program_name);
}

/* Parse comma-separated list of non-negative rates; returns count or -1 */
static int parse_rates(const char* str, double** rates) {
    int n = 1;
    for (const char* p = str; *p; p++)
        if (*p == ',') n++;
    *rates = malloc(n * sizeof(double));
    if (!*rates) return -1;
    const char* p = str;
    for (int i = 0; i < n; i++) {
        char* end;
        (*rates)[i] = strtod(p, &end);
        if (end == p || (*end != ',' && *end != '\0') || (*rates)[i] < 0) {
            free(*rates);
            *rates = NULL;
            return -1;
        }
        p = end + 1;
    }
    return n;
}

int main(int argc, char* argv[]) {
    /* Command line options */
    int n_chromosomes = DEFAULT_N_CHROMOSOMES;
    double rec_rate = DEFAULT_REC_RATE;
    unsigned long seed = 0;  /* 0 means use time-based seed */
    const char* output_file = NULL;
    const char* ancestor_file = NULL;
    int all_individuals = 0;
    double* rates = NULL;
    int n_rates = 0;
    int verbose = 0;

    /* Parse command line options */
    int opt;
    while ((opt = getopt(argc, argv, "c:r:R:s:o:A:avh")) != -1) {
        switch (opt) {
            case 'c':
                n_chromosomes = atoi(optarg);
                if (n_chromosomes < 1) {
                    fprintf(stderr, "Error: Number of chromosomes must be at least 1\n");
                    return 1;
                }
                break;
            case 'r':
                rec_rate = atof(optarg);
                if (rec_rate < 0) {
                    fprintf(stderr, "Error: Recombination rate must be non-negative\n");
                    return 1;
                }
                break;
            case 'R':
                free(rates);
                n_rates = parse_rates(optarg, &rates);
                if (n_rates < 1) {
                    fprintf(stderr, "Error: Invalid recombination rate list '%s'\n", optarg);
                    return 1;
                }
                break;
            case 'A':
                ancestor_file = optarg;
                break;
            case 'a':
                all_individuals = 1;
                break;
            case 's':
                seed = strtoul(optarg, NULL, 10);
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    /* Check for pedigree file argument */
    if (optind >= argc) {
        fprintf(stderr, "Error: No pedigree file specified\n\n");
        print_usage(argv[0]);
        return 1;
    }

    const char* ped_file = argv[optind];

    if (rates) {
        n_chromosomes = n_rates;
    }

    /* Use time-based seed if not specified */
    if (seed == 0) {
        seed = (unsigned long)time(NULL);
    }

    if (verbose) {
        fprintf(stderr, "pedtrans - Pedigree Chromosome Transmission Simulator\n");
        fprintf(stderr, "Pedigree file: %s\n", ped_file);
        fprintf(stderr, "Chromosomes: %d\n", n_chromosomes);
        fprintf(stderr, "Recombination rate: %.4f\n", rec_rate);
        fprintf(stderr, "Random seed: %lu\n", seed);
        if (output_file) {
            fprintf(stderr, "Output file: %s\n", output_file);
        }
        fprintf(stderr, "\n");
    }

    /* Run simulation */
    if (verbose) {
        fprintf(stderr, "Running simulation...\n");
    }

    ped_simulation* sim = rates
        ? simulate_pedigree_rates(ped_file, n_chromosomes, rates, seed)
        : simulate_pedigree(ped_file, n_chromosomes, rec_rate, seed);
    free(rates);
    if (!sim) {
        fprintf(stderr, "Error: Simulation failed\n");
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Simulation complete.\n");
        fprintf(stderr, "  Individuals: %d\n", sim->ped->n_individuals);
        fprintf(stderr, "  Founders: %d\n", sim->ped->n_founders);
        fprintf(stderr, "  Chromosomes: %d\n", sim->n_chromosomes);
        fprintf(stderr, "\n");
    }

    if (ancestor_file) {
        int to_stdout = strcmp(ancestor_file, "-") == 0;
        FILE* af = to_stdout ? stdout : fopen(ancestor_file, "w");
        if (!af) {
            fprintf(stderr, "Error: Cannot open ancestor file '%s'\n", ancestor_file);
            free_simulation(sim);
            return 1;
        }
        print_ancestor_summary(sim, af, all_individuals);
        if (!to_stdout) fclose(af);
        if (!output_file) {
            free_simulation(sim);
            return 0;
        }
    }

    /* Open output file */
    FILE* out = stdout;
    if (output_file) {
        out = fopen(output_file, "w");
        if (!out) {
            fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
            free_simulation(sim);
            return 1;
        }
    }

    /* Print results */
    print_all_segments(sim, out);

    /* Cleanup */
    if (output_file) {
        fclose(out);
        if (verbose) {
            fprintf(stderr, "Output written to %s\n", output_file);
        }
    }

    free_simulation(sim);

    return 0;
}
