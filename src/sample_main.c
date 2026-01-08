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
    fprintf(stderr, "Usage: %s [OPTIONS] INPUT.vcf [INPUT2.vcf [INPUT3.vcf]]\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Subsample individuals and/or markers from VCF file(s).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o FILE    Output VCF file (single-file mode, default: stdout)\n");
    fprintf(stderr, "  -O SUFFIX  Output suffix for multi-file mode (e.g., \"_filtered\")\n");
    fprintf(stderr, "             Creates INPUT_SUFFIX.vcf for each input file\n");
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
    fprintf(stderr, "Single-file examples:\n");
    fprintf(stderr, "  %s -I samples.txt input.vcf -o output.vcf\n", program_name);
    fprintf(stderr, "    Keep only individuals listed in samples.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s -n 100 -m u input.vcf -o output.vcf\n", program_name);
    fprintf(stderr, "    Select 100 uniformly-spaced markers per chromosome\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Multi-file examples (marker intersection):\n");
    fprintf(stderr, "  %s -n 100 -m u -O \"_common\" pop1.vcf pop2.vcf\n", program_name);
    fprintf(stderr, "    Select 100 markers shared between pop1 and pop2\n");
    fprintf(stderr, "    Output: pop1_common.vcf, pop2_common.vcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s -n 50 -O \"_filt\" pop1.vcf pop2.vcf pop3.vcf\n", program_name);
    fprintf(stderr, "    Select 50 markers shared across all 3 files\n");
}

int main(int argc, char* argv[]) {
    sample_params_t params;
    init_sample_params(&params);

    /* Parse command line options */
    int opt;
    while ((opt = getopt(argc, argv, "o:O:I:n:m:R:s:vh")) != -1) {
        switch (opt) {
            case 'o':
                params.output_file = optarg;
                break;
            case 'O':
                params.output_suffix = optarg;
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

    /* Check for input file argument(s) */
    if (optind >= argc) {
        fprintf(stderr, "Error: No input VCF file specified\n\n");
        print_usage(argv[0]);
        return 1;
    }

    /* Count input files */
    int n_files = argc - optind;
    if (n_files > 3) {
        fprintf(stderr, "Error: Maximum 3 input files supported\n");
        return 1;
    }

    if (n_files == 1) {
        /* Single-file mode */
        params.input_file = argv[optind];
        params.n_input_files = 0;
    } else {
        /* Multi-file mode */
        params.input_files = (const char**)(argv + optind);
        params.n_input_files = n_files;
    }

    /* Validate options */
    if (n_files > 1) {
        /* Multi-file mode validations */
        if (!params.output_suffix) {
            fprintf(stderr, "Error: -O (output suffix) required for multi-file mode\n");
            fprintf(stderr, "Example: -O \"_filtered\"\n\n");
            print_usage(argv[0]);
            return 1;
        }
        if (params.output_file) {
            fprintf(stderr, "Error: -o and -O are mutually exclusive\n");
            fprintf(stderr, "Use -o for single file, -O for multiple files\n");
            return 1;
        }
        if (params.n_markers == 0) {
            fprintf(stderr, "Error: -n (number of markers) required for multi-file mode\n");
            return 1;
        }
    } else {
        /* Single-file mode validations */
        if (params.output_suffix) {
            fprintf(stderr, "Error: -O requires multiple input files\n");
            fprintf(stderr, "Use -o for single file output\n");
            return 1;
        }
        if (!params.indiv_file && params.n_markers == 0 && !params.region.chrom) {
            fprintf(stderr, "Error: No subsetting options specified.\n");
            fprintf(stderr, "Use -I for individuals, -n for markers, or -R for region.\n\n");
            print_usage(argv[0]);
            return 1;
        }
    }

    if (params.verbose) {
        fprintf(stderr, "sample - VCF Subsampling Tool\n");
        if (params.n_input_files > 0) {
            fprintf(stderr, "Mode: multi-file (%d files)\n", params.n_input_files);
            for (int i = 0; i < params.n_input_files; i++) {
                fprintf(stderr, "  Input %d: %s\n", i + 1, params.input_files[i]);
            }
            fprintf(stderr, "Output suffix: %s\n", params.output_suffix);
        } else {
            fprintf(stderr, "Input file: %s\n", params.input_file);
            if (params.output_file) {
                fprintf(stderr, "Output file: %s\n", params.output_file);
            } else {
                fprintf(stderr, "Output: stdout\n");
            }
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
    int result;
    if (params.n_input_files > 0) {
        result = run_sample_multi(&params);
    } else {
        result = run_sample(&params);
    }

    /* Cleanup */
    free_region(&params.region);

    return result;
}
