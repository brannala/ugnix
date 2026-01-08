/*
 * pedsim_vcf_multipop_main.c - CLI for multi-population VCF pipeline
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "pedsim_vcf_multipop.h"

/* Parse comma-separated integer list */
static int* parse_int_list(const char* str, int* count)
{
    if (!str || !count) return NULL;

    int n = 1;
    for (const char* p = str; *p; p++) {
        if (*p == ',') n++;
    }

    int* arr = malloc(n * sizeof(int));
    if (!arr) return NULL;

    char* copy = strdup(str);
    char* tok = strtok(copy, ",");
    int i = 0;
    while (tok && i < n) {
        arr[i++] = atoi(tok);
        tok = strtok(NULL, ",");
    }

    free(copy);
    *count = i;
    return arr;
}

/* Parse migration matrix string */
static double** parse_migration(const char* str, int n_pops)
{
    if (!str || n_pops <= 0) return NULL;

    double** mig = malloc(n_pops * sizeof(double*));
    if (!mig) return NULL;

    for (int i = 0; i < n_pops; i++) {
        mig[i] = calloc(n_pops, sizeof(double));
        if (!mig[i]) {
            for (int j = 0; j < i; j++) free(mig[j]);
            free(mig);
            return NULL;
        }
    }

    char* copy = strdup(str);
    char* row_str = strtok(copy, ";");
    int row = 0;
    while (row_str && row < n_pops) {
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

static void free_migration(double** mig, int n_pops)
{
    if (!mig) return;
    for (int i = 0; i < n_pops; i++) {
        free(mig[i]);
    }
    free(mig);
}

static void print_usage(const char* progname)
{
    fprintf(stderr, "Usage: %s [options] -o <output.vcf>\n", progname);
    fprintf(stderr, "\nMulti-population pedigree-coalescent VCF simulation\n\n");
    fprintf(stderr, "Population options:\n");
    fprintf(stderr, "  -J NUM    Number of populations (default: 1)\n");
    fprintf(stderr, "  -N LIST   Population sizes, comma-separated (e.g., 1000,500,2000)\n");
    fprintf(stderr, "  -n LIST   Sample sizes, comma-separated (e.g., 10,5,20)\n");
    fprintf(stderr, "  -g NUM    Number of generations to trace back (default: 5)\n");
    fprintf(stderr, "\nMigration options:\n");
    fprintf(stderr, "  -m STR    Full migration matrix, semicolon-separated rows\n");
    fprintf(stderr, "  -M NUM    Symmetric island model migration rate\n");
    fprintf(stderr, "\nCoalescent options:\n");
    fprintf(stderr, "  -E NUM    Effective population size for coalescent (default: 10000)\n");
    fprintf(stderr, "  -r NUM    Recombination rate in cM (default: 1.0)\n");
    fprintf(stderr, "  -u NUM    Mutation rate (default: 0.5)\n");
    fprintf(stderr, "  -s NUM    Random seed (default: 0 = random)\n");
    fprintf(stderr, "\nOutput options:\n");
    fprintf(stderr, "  -o FILE   Output VCF file (required)\n");
    fprintf(stderr, "  -S        Output samples only (exclude founders)\n");
    fprintf(stderr, "  -k        Keep temporary files\n");
    fprintf(stderr, "  -v        Verbose output\n");
    fprintf(stderr, "  -h        Show this help message\n");
    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  # Two populations with island model migration\n");
    fprintf(stderr, "  %s -J 2 -N 1000,500 -n 10,5 -M 0.01 -g 5 -E 10000 -o output.vcf\n", progname);
    fprintf(stderr, "\n  # Three populations with custom migration\n");
    fprintf(stderr, "  %s -J 3 -N 1000,500,2000 -n 10,5,20 \\\n", progname);
    fprintf(stderr, "     -m \"0,0.01,0.02;0.01,0,0.01;0.02,0.01,0\" -o output.vcf\n");
}

int main(int argc, char** argv)
{
    multipop_vcf_pipeline_params params = {
        .n_populations = 1,
        .populations = NULL,
        .migration = NULL,
        .island_m = -1.0,
        .n_generations = 5,
        .coal_pop_size = 10000,
        .rec_rate = 1.0,
        .mut_rate = 0.5,
        .seed = 0,
        .output_vcf = NULL,
        .samples_only = 0,
        .keep_temp = 0,
        .verbose = 0
    };

    char* pop_sizes_str = NULL;
    char* sample_sizes_str = NULL;
    char* mig_str = NULL;

    int c;
    while ((c = getopt(argc, argv, "J:N:n:g:m:M:E:r:u:s:o:Skvh")) != -1) {
        switch (c) {
            case 'J':
                params.n_populations = atoi(optarg);
                break;
            case 'N':
                pop_sizes_str = optarg;
                break;
            case 'n':
                sample_sizes_str = optarg;
                break;
            case 'g':
                params.n_generations = atoi(optarg);
                break;
            case 'm':
                mig_str = optarg;
                break;
            case 'M':
                params.island_m = atof(optarg);
                break;
            case 'E':
                params.coal_pop_size = atof(optarg);
                break;
            case 'r':
                params.rec_rate = atof(optarg);
                break;
            case 'u':
                params.mut_rate = atof(optarg);
                break;
            case 's':
                params.seed = atoi(optarg);
                break;
            case 'o':
                params.output_vcf = optarg;
                break;
            case 'S':
                params.samples_only = 1;
                break;
            case 'k':
                params.keep_temp = 1;
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

    /* Validate required options */
    if (!params.output_vcf) {
        fprintf(stderr, "Error: output VCF file required (-o)\n");
        print_usage(argv[0]);
        return 1;
    }

    if (params.n_populations < 1 || params.n_populations > MAX_POPULATIONS) {
        fprintf(stderr, "Error: number of populations must be 1-%d\n",
                MAX_POPULATIONS);
        return 1;
    }

    if (mig_str && params.island_m >= 0) {
        fprintf(stderr, "Error: cannot specify both -m and -M\n");
        return 1;
    }

    /* Parse population sizes */
    int* pop_sizes = NULL;
    int n_pop_sizes = 0;
    if (pop_sizes_str) {
        pop_sizes = parse_int_list(pop_sizes_str, &n_pop_sizes);
        if (!pop_sizes) {
            fprintf(stderr, "Error: failed to parse population sizes\n");
            return 1;
        }
    }

    /* Parse sample sizes */
    int* sample_sizes = NULL;
    int n_sample_sizes = 0;
    if (sample_sizes_str) {
        sample_sizes = parse_int_list(sample_sizes_str, &n_sample_sizes);
        if (!sample_sizes) {
            fprintf(stderr, "Error: failed to parse sample sizes\n");
            free(pop_sizes);
            return 1;
        }
    }

    /* Build population parameters */
    params.populations = malloc(params.n_populations * sizeof(pop_vcf_params));
    if (!params.populations) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(pop_sizes);
        free(sample_sizes);
        return 1;
    }

    for (int i = 0; i < params.n_populations; i++) {
        if (pop_sizes) {
            if (n_pop_sizes == 1) {
                params.populations[i].pop_size = pop_sizes[0];
            } else if (i < n_pop_sizes) {
                params.populations[i].pop_size = pop_sizes[i];
            } else {
                fprintf(stderr, "Error: not enough population sizes specified\n");
                goto cleanup_error;
            }
        } else {
            params.populations[i].pop_size = 1000;
        }

        if (sample_sizes) {
            if (n_sample_sizes == 1) {
                params.populations[i].sample_size = sample_sizes[0];
            } else if (i < n_sample_sizes) {
                params.populations[i].sample_size = sample_sizes[i];
            } else {
                fprintf(stderr, "Error: not enough sample sizes specified\n");
                goto cleanup_error;
            }
        } else {
            params.populations[i].sample_size = 10;
        }
    }

    free(pop_sizes);
    free(sample_sizes);
    pop_sizes = sample_sizes = NULL;

    /* Parse migration matrix if provided */
    if (mig_str) {
        params.migration = parse_migration(mig_str, params.n_populations);
        if (!params.migration) {
            fprintf(stderr, "Error: failed to parse migration matrix\n");
            free(params.populations);
            return 1;
        }
    }

    /* Print parameters */
    fprintf(stderr, "Multi-population Pedigree-Coalescent VCF Simulation\n");
    fprintf(stderr, "===================================================\n");
    fprintf(stderr, "Populations: %d\n", params.n_populations);
    for (int i = 0; i < params.n_populations; i++) {
        fprintf(stderr, "  Pop%d: N=%d, n=%d\n", i,
                params.populations[i].pop_size,
                params.populations[i].sample_size);
    }
    fprintf(stderr, "Generations: %d\n", params.n_generations);

    if (params.migration) {
        fprintf(stderr, "Migration matrix:\n");
        for (int i = 0; i < params.n_populations; i++) {
            fprintf(stderr, "  ");
            for (int j = 0; j < params.n_populations; j++) {
                fprintf(stderr, "%.4f ", params.migration[i][j]);
            }
            fprintf(stderr, "\n");
        }
    } else if (params.island_m >= 0) {
        fprintf(stderr, "Island model migration: %.4f\n", params.island_m);
    } else {
        fprintf(stderr, "Migration: none\n");
    }

    fprintf(stderr, "Coalescent: Ne=%.0f, r=%.2f cM, m=%.2f\n",
            params.coal_pop_size, params.rec_rate, params.mut_rate);

    /* Run pipeline */
    int result = run_multipop_vcf_pipeline(&params);

    /* Cleanup */
    if (params.migration) {
        free_migration(params.migration, params.n_populations);
    }
    free(params.populations);

    return result;

cleanup_error:
    free(pop_sizes);
    free(sample_sizes);
    free(params.populations);
    return 1;
}
