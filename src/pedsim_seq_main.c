/*
 * pedsim_seq_main.c - CLI for Integrated Pedigree-Coalescent Simulation
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "pedsim_seq.h"

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [OPTIONS]\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Integrated pedigree-coalescent sequence simulation.\n");
    fprintf(stderr, "Combines pedsim, coalsim, pedtrans, and seqassemble.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Pedigree options:\n");
    fprintf(stderr, "  -n NUM     Number of sample individuals (default: 10)\n");
    fprintf(stderr, "  -N NUM     Effective population size (default: 1000)\n");
    fprintf(stderr, "  -k NUM     Number of generations (default: 5)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Coalescent options:\n");
    fprintf(stderr, "  -r RATE    Recombination rate per chromosome (default: 0.1)\n");
    fprintf(stderr, "  -m RATE    Mutation rate per chromosome (default: 0.5)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Transmission options:\n");
    fprintf(stderr, "  -R RATE    Recombination rate for pedtrans (default: 1.0)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "  -o FILE    Output file (default: stdout)\n");
    fprintf(stderr, "  -f FORMAT  Output format: fasta (default) or vcf\n");
    fprintf(stderr, "  -S         Samples only (exclude founders from output)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, "  -s SEED    Random seed (default: time-based)\n");
    fprintf(stderr, "  -v         Verbose output\n");
    fprintf(stderr, "  -h         Show this help message\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "  %s -n 20 -N 1000 -k 5 -r 0.1 -m 0.5 -R 1.0 -o output.fa\n",
            program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "This tool performs a complete simulation pipeline:\n");
    fprintf(stderr, "  1. Generates a random pedigree (pedsim)\n");
    fprintf(stderr, "  2. Simulates founder sequences via coalescent (coalsim)\n");
    fprintf(stderr, "  3. Transmits chromosomes through the pedigree (pedtrans)\n");
    fprintf(stderr, "  4. Assembles sample sequences from founder segments\n");
}

int main(int argc, char* argv[]) {
    pedsim_seq_params params;
    init_pedsim_seq_params(&params);

    int opt;
    while ((opt = getopt(argc, argv, "n:N:k:r:m:R:o:f:Ss:vh")) != -1) {
        switch (opt) {
            case 'n':
                params.n_samples = atoi(optarg);
                if (params.n_samples < 1) {
                    fprintf(stderr, "Error: n must be positive\n");
                    return 1;
                }
                break;
            case 'N':
                params.N = atoi(optarg);
                if (params.N < 1) {
                    fprintf(stderr, "Error: N must be positive\n");
                    return 1;
                }
                break;
            case 'k':
                params.k = atoi(optarg);
                if (params.k < 1) {
                    fprintf(stderr, "Error: k must be positive\n");
                    return 1;
                }
                break;
            case 'r':
                params.rec_rate = atof(optarg);
                if (params.rec_rate < 0) {
                    fprintf(stderr, "Error: r must be non-negative\n");
                    return 1;
                }
                break;
            case 'm':
                params.mut_rate = atof(optarg);
                if (params.mut_rate < 0) {
                    fprintf(stderr, "Error: m must be non-negative\n");
                    return 1;
                }
                break;
            case 'R':
                params.pedtrans_rec = atof(optarg);
                if (params.pedtrans_rec < 0) {
                    fprintf(stderr, "Error: R must be non-negative\n");
                    return 1;
                }
                break;
            case 'o':
                params.output_file = optarg;
                break;
            case 'f':
                params.format = optarg;
                if (strcmp(params.format, "fasta") != 0 &&
                    strcmp(params.format, "vcf") != 0) {
                    fprintf(stderr, "Error: Unknown format '%s'\n", params.format);
                    return 1;
                }
                break;
            case 'S':
                params.samples_only = 1;
                break;
            case 's':
                params.seed = strtoul(optarg, NULL, 10);
                break;
            case 'v':
                params.verbose = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    return run_pedsim_seq(&params);
}
