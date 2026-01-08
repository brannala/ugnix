/*
 * pedsim_vcf_multipop.c - Multi-population VCF Simulation Pipeline
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include "pedsim_vcf_multipop.h"

/* Run a command and return exit status */
static int run_command(const char* cmd, int verbose)
{
    if (verbose) {
        fprintf(stderr, "Running: %s\n", cmd);
    }
    int ret = system(cmd);
    if (ret == -1) {
        fprintf(stderr, "Error: failed to execute command\n");
        return -1;
    }
    return WEXITSTATUS(ret);
}

/* Count founders in multi-population pedigree file */
static int count_founders(const char* ped_file)
{
    FILE* fp = fopen(ped_file, "r");
    if (!fp) return -1;

    int count = 0;
    char line[1024];

    while (fgets(line, sizeof(line), fp)) {
        /* Skip comments and empty lines */
        if (line[0] == '#' || line[0] == '\n') continue;

        char id[256], father[256], mother[256];
        int pop_id;

        /* Try 4-column format first (with population ID) */
        int n = sscanf(line, "%s %s %s %d", id, father, mother, &pop_id);
        if (n >= 3) {
            if (strcmp(father, "0") == 0 && strcmp(mother, "0") == 0) {
                count++;
            }
        }
    }

    fclose(fp);
    return count;
}

/* Build migration matrix string for command line */
static char* build_migration_string(double** migration, int n_pops)
{
    /* Estimate size: each entry is at most 10 chars, plus separators */
    int buf_size = n_pops * n_pops * 12 + n_pops * 2 + 10;
    char* str = malloc(buf_size);
    if (!str) return NULL;

    str[0] = '\0';
    char temp[32];

    for (int i = 0; i < n_pops; i++) {
        if (i > 0) strcat(str, ";");
        for (int j = 0; j < n_pops; j++) {
            if (j > 0) strcat(str, ",");
            snprintf(temp, sizeof(temp), "%.6f", migration[i][j]);
            strcat(str, temp);
        }
    }

    return str;
}

/* Build population sizes string */
static char* build_pop_sizes_string(multipop_vcf_pipeline_params* params)
{
    int buf_size = params->n_populations * 16 + 16;
    char* str = malloc(buf_size);
    if (!str) return NULL;

    str[0] = '\0';
    char temp[32];

    for (int i = 0; i < params->n_populations; i++) {
        if (i > 0) strcat(str, ",");
        snprintf(temp, sizeof(temp), "%d", params->populations[i].pop_size);
        strcat(str, temp);
    }

    return str;
}

/* Build sample sizes string */
static char* build_sample_sizes_string(multipop_vcf_pipeline_params* params)
{
    int buf_size = params->n_populations * 16 + 16;
    char* str = malloc(buf_size);
    if (!str) return NULL;

    str[0] = '\0';
    char temp[32];

    for (int i = 0; i < params->n_populations; i++) {
        if (i > 0) strcat(str, ",");
        snprintf(temp, sizeof(temp), "%d", params->populations[i].sample_size);
        strcat(str, temp);
    }

    return str;
}

int run_multipop_vcf_pipeline(multipop_vcf_pipeline_params* params)
{
    char cmd[8192];
    char* temp_dir = NULL;
    int ret = 0;

    /* Create temporary directory */
    char template[] = "/tmp/pedsim_vcf_multipop_XXXXXX";
    temp_dir = mkdtemp(template);
    if (!temp_dir) {
        fprintf(stderr, "Error: failed to create temp directory\n");
        return 1;
    }

    char ped_file[512], vcf_file[512], seg_file[512];
    snprintf(ped_file, sizeof(ped_file), "%s/ped.txt", temp_dir);
    snprintf(vcf_file, sizeof(vcf_file), "%s/founder.vcf", temp_dir);
    snprintf(seg_file, sizeof(seg_file), "%s/segments.txt", temp_dir);

    /* Build parameter strings */
    char* pop_sizes = build_pop_sizes_string(params);
    char* sample_sizes = build_sample_sizes_string(params);
    char* mig_str = NULL;

    if (!pop_sizes || !sample_sizes) {
        fprintf(stderr, "Error: failed to allocate parameter strings\n");
        ret = 1;
        goto cleanup;
    }

    /* Step 1: Run pedsim_multipop to generate pedigree */
    fprintf(stderr, "\n=== Step 1: Generating multi-population pedigree ===\n");

    if (params->migration) {
        mig_str = build_migration_string(params->migration, params->n_populations);
        if (!mig_str) {
            ret = 1;
            goto cleanup;
        }
        snprintf(cmd, sizeof(cmd),
                 "./pedsim_multipop -J %d -N %s -n %s -k %d -m \"%s\" -s %u -o %s",
                 params->n_populations, pop_sizes, sample_sizes,
                 params->n_generations, mig_str, params->seed, ped_file);
    } else if (params->island_m >= 0) {
        snprintf(cmd, sizeof(cmd),
                 "./pedsim_multipop -J %d -N %s -n %s -k %d -M %.6f -s %u -o %s",
                 params->n_populations, pop_sizes, sample_sizes,
                 params->n_generations, params->island_m, params->seed, ped_file);
    } else {
        snprintf(cmd, sizeof(cmd),
                 "./pedsim_multipop -J %d -N %s -n %s -k %d -s %u -o %s",
                 params->n_populations, pop_sizes, sample_sizes,
                 params->n_generations, params->seed, ped_file);
    }

    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: pedsim_multipop failed\n");
        goto cleanup;
    }

    /* Count founders */
    int n_founders = count_founders(ped_file);
    if (n_founders <= 0) {
        fprintf(stderr, "Error: failed to count founders\n");
        ret = 1;
        goto cleanup;
    }
    int n_haplotypes = n_founders * 2;
    fprintf(stderr, "Pedigree has %d founders (%d haplotypes) across %d populations\n",
            n_founders, n_haplotypes, params->n_populations);

    /* Step 2: Run coalsim to simulate founder haplotypes */
    fprintf(stderr, "\n=== Step 2: Simulating founder haplotypes (VCF) ===\n");
    snprintf(cmd, sizeof(cmd),
             "./coalsim -c %d -N %.0f -r %.4f -m %.4f -s %u -u b -V %s",
             n_haplotypes, params->coal_pop_size, params->rec_rate,
             params->mut_rate, params->seed + 1, vcf_file);
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: coalsim failed\n");
        goto cleanup;
    }

    /* Step 3: Run pedtrans to simulate chromosome transmission */
    fprintf(stderr, "\n=== Step 3: Simulating chromosome transmission ===\n");
    snprintf(cmd, sizeof(cmd),
             "./pedtrans -r %.4f -s %u -o %s %s",
             params->rec_rate, params->seed + 2, seg_file, ped_file);
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: pedtrans failed\n");
        goto cleanup;
    }

    /* Step 4: Run vcfassemble to create sample VCF */
    fprintf(stderr, "\n=== Step 4: Assembling sample VCF ===\n");
    snprintf(cmd, sizeof(cmd),
             "./vcfassemble -v %s -p %s -o %s%s",
             vcf_file, seg_file, params->output_vcf,
             params->samples_only ? " -s" : "");
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: vcfassemble failed\n");
        goto cleanup;
    }

    fprintf(stderr, "\n=== Pipeline complete ===\n");
    fprintf(stderr, "Output VCF: %s\n", params->output_vcf);

    /* Count total samples */
    int total_samples = 0;
    for (int i = 0; i < params->n_populations; i++) {
        total_samples += params->populations[i].sample_size;
    }
    fprintf(stderr, "Total samples: %d across %d populations\n",
            total_samples, params->n_populations);

cleanup:
    free(pop_sizes);
    free(sample_sizes);
    free(mig_str);

    /* Clean up temporary files */
    if (!params->keep_temp && temp_dir) {
        snprintf(cmd, sizeof(cmd), "rm -rf %s", temp_dir);
        int cleanup_ret = system(cmd);
        (void)cleanup_ret;
    } else if (params->keep_temp) {
        fprintf(stderr, "Temp files kept in: %s\n", temp_dir);
    }

    return ret;
}
