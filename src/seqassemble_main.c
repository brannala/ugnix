/*
 * seqassemble_main.c - CLI for Sequence Assembly
 *
 * Combines founder sequences (from coalsim) with segment transmission data
 * (from pedtrans) to assemble final sample sequences.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "seqassemble.h"

static void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [OPTIONS] SEGMENTS_FILE FOUNDERS_FASTA\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "Assemble sample sequences from pedigree segment data and founder sequences.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  SEGMENTS_FILE   Output from pedtrans (segment transmission data)\n");
    fprintf(stderr, "  FOUNDERS_FASTA  Output from coalsim (founder sequences in FASTA format)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o FILE    Output file (default: stdout)\n");
    fprintf(stderr, "  -f FORMAT  Output format: fasta (default) or vcf\n");
    fprintf(stderr, "  -s         Samples only (exclude founders from output)\n");
    fprintf(stderr, "  -m         Memory-efficient streaming mode (for large datasets)\n");
    fprintf(stderr, "  -v         Verbose output\n");
    fprintf(stderr, "  -h         Show this help message\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The tool combines:\n");
    fprintf(stderr, "  1. Segment data from pedtrans showing founder origin of each region\n");
    fprintf(stderr, "  2. Founder sequences from coalsim\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "For each sample individual, it assembles the full sequence by copying\n");
    fprintf(stderr, "the appropriate founder sequence for each inherited segment.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example (complete pipeline):\n");
    fprintf(stderr, "  # Step 1: Generate pedigree\n");
    fprintf(stderr, "  ./pedsim -n 10 -N 1000 -k 5 -s 12345 -o pedigree.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  # Step 2: Count founders (individuals with \"0 0\" parents)\n");
    fprintf(stderr, "  F=$(grep \" 0 0$\" pedigree.txt | wc -l)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  # Step 3: Generate founder sequences (2 haplotypes per founder)\n");
    fprintf(stderr, "  ./coalsim -c $((2*F)) -N 1000 -r 0.1 -m 0.5 -s 12345 -o founders.fa\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  # Step 4: Simulate chromosome transmission\n");
    fprintf(stderr, "  ./pedtrans -r 1.0 -s 12345 pedigree.txt -o segments.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  # Step 5: Assemble sample sequences\n");
    fprintf(stderr, "  ./%s segments.txt founders.fa -o samples.fa\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "  # Or output as VCF:\n");
    fprintf(stderr, "  ./%s -f vcf segments.txt founders.fa -o samples.vcf\n", program_name);
}

int main(int argc, char* argv[]) {
    /* Command line options */
    const char* output_file = NULL;
    const char* format = "fasta";
    int samples_only = 0;
    int streaming_mode = 0;
    int verbose = 0;

    /* Parse command line options */
    int opt;
    while ((opt = getopt(argc, argv, "o:f:smvh")) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            case 'f':
                format = optarg;
                if (strcmp(format, "fasta") != 0 && strcmp(format, "vcf") != 0) {
                    fprintf(stderr, "Error: Unknown format '%s'. Use 'fasta' or 'vcf'.\n",
                            format);
                    return 1;
                }
                break;
            case 's':
                samples_only = 1;
                break;
            case 'm':
                streaming_mode = 1;
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

    /* Check for required arguments */
    if (optind + 2 > argc) {
        fprintf(stderr, "Error: Missing required arguments\n\n");
        print_usage(argv[0]);
        return 1;
    }

    const char* segments_file = argv[optind];
    const char* founders_file = argv[optind + 1];

    if (verbose) {
        fprintf(stderr, "seqassemble - Sequence Assembly from Pedigree Segments\n");
        fprintf(stderr, "Segments file:  %s\n", segments_file);
        fprintf(stderr, "Founders file:  %s\n", founders_file);
        fprintf(stderr, "Output format:  %s\n", format);
        if (output_file) {
            fprintf(stderr, "Output file:    %s\n", output_file);
        }
        fprintf(stderr, "\n");
    }

    /* Parse pedtrans output */
    if (verbose) {
        fprintf(stderr, "Parsing segment data...\n");
    }

    pedtrans_data* pd = create_pedtrans_data();
    if (!pd) {
        fprintf(stderr, "Error: Failed to allocate pedtrans data\n");
        return 1;
    }

    if (parse_pedtrans_output(pd, segments_file) != 0) {
        fprintf(stderr, "Error: Failed to parse segments file\n");
        free_pedtrans_data(pd);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "  Individuals: %d\n", pd->n_individuals);
        fprintf(stderr, "  Founders: %d\n", pd->n_founders);
    }

    /* Open output file */
    FILE* out = stdout;
    if (output_file) {
        out = fopen(output_file, "w");
        if (!out) {
            fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
            free_pedtrans_data(pd);
            return 1;
        }
    }

    int result = 0;

    if (streaming_mode && strcmp(format, "fasta") == 0) {
        /* Streaming mode - memory efficient for large datasets */
        if (verbose) {
            fprintf(stderr, "Indexing founder sequences (streaming mode)...\n");
        }

        fasta_index* idx = create_fasta_index(founders_file,
                                               (const char**)pd->founder_names,
                                               pd->n_founders);
        if (!idx) {
            fprintf(stderr, "Error: Failed to index founder sequences\n");
            if (output_file) fclose(out);
            free_pedtrans_data(pd);
            return 1;
        }

        if (verbose) {
            fprintf(stderr, "  Sequence length: %ld bp\n", idx->seq_length);
            fprintf(stderr, "  Indexed entries: %d\n", idx->n_entries);
            fprintf(stderr, "Assembling sequences (streaming)...\n");
        }

        result = assemble_all_sequences_streaming(pd, idx, out, samples_only);

        if (result != 0) {
            fprintf(stderr, "Error: Failed to assemble sequences\n");
        } else if (verbose) {
            fprintf(stderr, "Done.\n");
            if (output_file) {
                fprintf(stderr, "Output written to %s\n", output_file);
            }
        }

        free_fasta_index(idx);

    } else {
        /* In-memory mode (default) */
        if (verbose) {
            fprintf(stderr, "Reading founder sequences...\n");
        }

        founder_sequences* fs = create_founder_sequences();
        if (!fs) {
            fprintf(stderr, "Error: Failed to allocate founder sequences\n");
            if (output_file) fclose(out);
            free_pedtrans_data(pd);
            return 1;
        }

        if (read_founder_fasta(fs, founders_file,
                               (const char**)pd->founder_names, pd->n_founders) != 0) {
            fprintf(stderr, "Error: Failed to read founder sequences\n");
            free_founder_sequences(fs);
            if (output_file) fclose(out);
            free_pedtrans_data(pd);
            return 1;
        }

        if (verbose) {
            fprintf(stderr, "  Sequence length: %ld bp\n", fs->seq_length);
            fprintf(stderr, "  Founders with sequences: %d\n", fs->n_founders);
            fprintf(stderr, "Assembling sequences...\n");
        }

        if (strcmp(format, "vcf") == 0) {
            result = write_vcf_output(pd, fs, out, samples_only);
        } else {
            result = assemble_all_sequences(pd, fs, out, samples_only);
        }

        if (result != 0) {
            fprintf(stderr, "Error: Failed to assemble sequences\n");
        } else if (verbose) {
            fprintf(stderr, "Done.\n");
            if (output_file) {
                fprintf(stderr, "Output written to %s\n", output_file);
            }
        }

        free_founder_sequences(fs);
    }

    /* Cleanup */
    if (output_file) {
        fclose(out);
    }
    free_pedtrans_data(pd);

    return result;
}
