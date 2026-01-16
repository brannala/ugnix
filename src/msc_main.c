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
    fprintf(stderr, "  -v, --verbose    Verbose output\n");
    fprintf(stderr, "  -h, --help       Print this help message\n");
    fprintf(stderr, "\nControl file format:\n");
    fprintf(stderr, "  # Comment lines start with #\n");
    fprintf(stderr, "  species_tree: ((A:1000#0.001,B:1000#0.002)AB:500#0.003,C:1500#0.001)ROOT#0.002;\n");
    fprintf(stderr, "  samples: A:10,B:10,C:10\n");
    fprintf(stderr, "  recombination_rate: 1.0\n");
    fprintf(stderr, "  mutation_rate: 0.5\n");
    fprintf(stderr, "  migration: A->B:0.5,B->A:0.3,AB->C:0.1\n");
    fprintf(stderr, "  seed: 12345\n");
    fprintf(stderr, "  output_gene_trees: true\n");
    fprintf(stderr, "  output_vcf: true\n");
    fprintf(stderr, "  vcf_file: output.vcf\n");
    fprintf(stderr, "\nNewick format: (children)NAME:branch_length#theta\n");
    fprintf(stderr, "  Internal nodes can have names (e.g., AB, ROOT) for migration\n");
    fprintf(stderr, "  branch_length = divergence time in generations\n");
    fprintf(stderr, "  theta = 4Nu (population-scaled mutation rate)\n");
    fprintf(stderr, "\nMigration format: source->dest:M\n");
    fprintf(stderr, "  M = 4Nm (population-scaled migration rate)\n");
    fprintf(stderr, "  Migration only between contemporary populations\n");
    fprintf(stderr, "  Can reference tip names or internal node names\n");
}

int main(int argc, char* argv[]) {
    int verbose = 0;

    static struct option long_options[] = {
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "vh", long_options, NULL)) != -1) {
        switch (c) {
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

    /* Print banner */
    printf("uGnix Multispecies Coalescent Simulator\n");
    printf("https://github.com/ugnix\n\n");

    /* Parse control file */
    msc_params* params = msc_parse_control_file(control_file);
    if (!params) {
        return 1;
    }

    /* Override verbose from command line */
    if (verbose) {
        params->verbose = 1;
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
