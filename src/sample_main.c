/*
 * sample_main.c - CLI for VCF Subsampling Tool
 *
 * Subsample individuals and/or markers from a VCF file.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>

#include "sample.h"

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [OPTIONS] INPUT.vcf\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Subsample individuals and/or markers from a VCF file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o FILE    Output VCF file (default: stdout)\n");
    fprintf(stderr, "  -I FILE    File with individual IDs to keep (one per line)\n");
    fprintf(stderr, "  -n NUM     Number of markers to select per chromosome\n");
    fprintf(stderr, "  -m METHOD  Marker selection method:\n");
    fprintf(stderr, "               u = uniform spacing (maximize inter-marker distance)\n");
    fprintf(stderr, "               r = random selection (default)\n");
    fprintf(stderr, "  -R REGION  Region to subset markers from:\n");
    fprintf(stderr, "               CHROM or CHROM:START-END\n");
    fprintf(stderr, "             If not specified, sampling applies to each chromosome\n");
    fprintf(stderr, "  -s SEED    Random seed for reproducibility\n");
    fprintf(stderr, "  -v         Verbose output\n");
    fprintf(stderr, "  -h         Show this help message\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s -I samples.txt input.vcf -o output.vcf\n", program_name);
    fprintf(stderr, "    Keep only individuals listed in samples.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s -n 100 -m u input.vcf -o output.vcf\n", program_name);
    fprintf(stderr, "    Select 100 uniformly-spaced markers per chromosome\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s -n 50 -m r -R chr1:1000000-5000000 input.vcf\n", program_name);
    fprintf(stderr, "    Randomly select 50 markers from chr1:1000000-5000000\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s -I samples.txt -n 100 -m u input.vcf -o output.vcf\n", program_name);
    fprintf(stderr, "    Subset both individuals and markers simultaneously\n");
}

int main(int argc, char* argv[]) {
    sample_params_t params;
    init_sample_params(&params);

    /* Parse command line options */
    int opt;
    while ((opt = getopt(argc, argv, "o:I:n:m:R:s:vh")) != -1) {
        switch (opt) {
            case 'o':
                params.output_file = optarg;
                break;
            case 'I':
                params.indiv_file = optarg;
                break;
            case 'n':
                params.n_markers = atoi(optarg);
                if (params.n_markers < 0) {
                    fprintf(stderr, "Error: Number of markers must be non-negative\n");
                    return 1;
                }
                break;
            case 'm':
                if (optarg[0] == 'u' || optarg[0] == 'U') {
                    params.method = METHOD_UNIFORM;
                } else if (optarg[0] == 'r' || optarg[0] == 'R') {
                    params.method = METHOD_RANDOM;
                } else {
                    fprintf(stderr, "Error: Unknown method '%s'. Use 'u' or 'r'.\n",
                            optarg);
                    return 1;
                }
                break;
            case 'R':
                if (parse_region(optarg, &params.region) != 0) {
                    fprintf(stderr, "Error: Invalid region '%s'\n", optarg);
                    return 1;
                }
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
            case '?':
                if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option '-%c'.\n", optopt);
                } else {
                    fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
                }
                print_usage(argv[0]);
                return 1;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    /* Check for input file argument */
    if (optind >= argc) {
        fprintf(stderr, "Error: No input VCF file specified\n\n");
        print_usage(argv[0]);
        return 1;
    }

    params.input_file = argv[optind];

    /* Validate options */
    if (!params.indiv_file && params.n_markers == 0 && !params.region.chrom) {
        fprintf(stderr, "Error: No subsetting options specified.\n");
        fprintf(stderr, "Use -I for individuals, -n for markers, or -R for region.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    if (params.verbose) {
        fprintf(stderr, "sample - VCF Subsampling Tool\n");
        fprintf(stderr, "Input file: %s\n", params.input_file);
        if (params.output_file) {
            fprintf(stderr, "Output file: %s\n", params.output_file);
        } else {
            fprintf(stderr, "Output: stdout\n");
        }
        if (params.indiv_file) {
            fprintf(stderr, "Individual list: %s\n", params.indiv_file);
        }
        if (params.n_markers > 0) {
            fprintf(stderr, "Markers per chromosome: %d\n", params.n_markers);
            fprintf(stderr, "Selection method: %s\n",
                    params.method == METHOD_UNIFORM ? "uniform" : "random");
        }
        if (params.region.chrom) {
            fprintf(stderr, "Region: %s", params.region.chrom);
            if (params.region.start > 0 || params.region.end > 0) {
                fprintf(stderr, ":%ld-%ld", params.region.start, params.region.end);
            }
            fprintf(stderr, "\n");
        }
        if (params.seed > 0) {
            fprintf(stderr, "Random seed: %lu\n", params.seed);
        }
        fprintf(stderr, "\n");
    }

    /* Run sampling */
    int result = run_sample(&params);

    /* Cleanup */
    free_region(&params.region);

    return result;
}
