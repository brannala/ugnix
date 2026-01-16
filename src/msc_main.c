/*
 * msc_main.c - Main entry point for multispecies coalescent simulation
 *
 * Usage: coalsim_msc <control_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "msc.h"
#include "species_tree.h"

/* Version string for uGnix banner */
char version[] = "coalsim_msc";

static void print_usage(const char* prog) {
    fprintf(stderr, "Usage: %s [options] <control_file>\n", prog);
    fprintf(stderr, "\nMultispecies Coalescent Simulation with Migration (MSC-M)\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "  -r, --reps N     Number of replicates (default: 1)\n");
    fprintf(stderr, "  -L, --tmrca      Output TMRCA and tree length (like ms -L)\n");
    fprintf(stderr, "  -q, --quiet      Suppress banner and per-replicate details\n");
    fprintf(stderr, "  -v, --verbose    Verbose output\n");
    fprintf(stderr, "  -h, --help       Print this help message\n");
    fprintf(stderr, "\nControl file format:\n");
    fprintf(stderr, "  # Comment lines start with #\n");
    fprintf(stderr, "  species_tree: ((A:1000#100,B:1000#200)AB:500#150,C:1500#100)ROOT#200;\n");
    fprintf(stderr, "  samples: A:10,B:10,C:10\n");
    fprintf(stderr, "  recombination_rate: 1.0\n");
    fprintf(stderr, "  mutation_rate: 0.5\n");
    fprintf(stderr, "  migration: A->B:0.001,B->A:0.001,AB->C:0.0005\n");
    fprintf(stderr, "  seed: 12345\n");
    fprintf(stderr, "  output_gene_trees: true\n");
    fprintf(stderr, "  output_vcf: true\n");
    fprintf(stderr, "  vcf_file: output.vcf\n");
    fprintf(stderr, "\nNewick format: (children)NAME:branch_length#N\n");
    fprintf(stderr, "  N = diploid effective population size\n");
    fprintf(stderr, "  branch_length = divergence time in generations\n");
    fprintf(stderr, "  Internal nodes can have names (e.g., AB, ROOT) for migration\n");
    fprintf(stderr, "\nCoalescence rate: n(n-1)/(4N) for n lineages in population of size N\n");
    fprintf(stderr, "\nMutation rate: per chromosome per generation\n");
    fprintf(stderr, "\nMigration format: source->dest:m\n");
    fprintf(stderr, "  m = per-lineage migration rate per generation\n");
    fprintf(stderr, "  Migration only between contemporary populations\n");
    fprintf(stderr, "  Can reference tip names or internal node names\n");
}

int main(int argc, char* argv[]) {
    int verbose = 0;
    int quiet = 0;
    int output_tmrca = 0;
    int nreps = 0;  /* 0 means use control file value or default */

    static struct option long_options[] = {
        {"reps", required_argument, 0, 'r'},
        {"tmrca", no_argument, 0, 'L'},
        {"quiet", no_argument, 0, 'q'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "r:Lqvh", long_options, NULL)) != -1) {
        switch (c) {
            case 'r':
                nreps = atoi(optarg);
                if (nreps <= 0) {
                    fprintf(stderr, "Error: Number of replicates must be positive\n");
                    return 1;
                }
                break;
            case 'L':
                output_tmrca = 1;
                break;
            case 'q':
                quiet = 1;
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

    if (optind >= argc) {
        fprintf(stderr, "Error: Control file required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    const char* control_file = argv[optind];

    /* Print banner (unless quiet) */
    if (!quiet) {
        printf("uGnix Multispecies Coalescent Simulator\n");
        printf("https://github.com/ugnix\n\n");
    }

    /* Parse control file */
    msc_params* params = msc_parse_control_file(control_file);
    if (!params) {
        return 1;
    }

    /* Override parameters from command line */
    if (verbose) {
        params->verbose = 1;
    }
    if (quiet) {
        params->quiet = 1;
    }
    if (output_tmrca) {
        params->output_tmrca = 1;
    }
    if (nreps > 0) {
        params->nreps = nreps;
    }

    /* Validate parameters */
    if (msc_validate_params(params) < 0) {
        msc_free_params(params);
        return 1;
    }

    /* Print summary */
    if (params->verbose) {
        msc_print_summary(params, stderr);
        fprintf(stderr, "\n");
    }

    /* Run simulation */
    int result = msc_simulate(params);

    /* Cleanup */
    msc_free_params(params);

    return result;
}
