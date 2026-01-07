/*
 * pedsim_seq.c - Integrated Pedigree-Coalescent Sequence Simulation
 *
 * Uses a pipeline approach: runs pedsim, coalsim, pedtrans, and seqassemble
 * in sequence using temporary files.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include "pedsim_seq.h"

void init_pedsim_seq_params(pedsim_seq_params* params) {
    params->n_samples = 10;
    params->N = 1000;
    params->k = 5;
    params->rec_rate = 0.1;
    params->mut_rate = 0.5;
    params->seq_length = 0;  /* use coalsim default */
    params->pedtrans_rec = 1.0;
    params->output_file = NULL;
    params->format = "fasta";
    params->samples_only = 0;
    params->verbose = 0;
    params->seed = 0;
}

/*
 * Run a command and return exit status.
 */
static int run_command(const char* cmd, int verbose) {
    if (verbose) {
        fprintf(stderr, "  Running: %s\n", cmd);
    }
    int status = system(cmd);
    if (status == -1) {
        fprintf(stderr, "Error: Failed to execute command\n");
        return -1;
    }
    return WEXITSTATUS(status);
}

/*
 * Count founders in pedigree file.
 * Founders are lines ending with "0 0".
 */
static int count_founders(const char* ped_file) {
    FILE* fp = fopen(ped_file, "r");
    if (!fp) return -1;

    int count = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;  /* skip comments */
        /* Check if line ends with " 0 0" (founder pattern) */
        int len = strlen(line);
        /* Remove trailing newline */
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }
        if (len >= 4) {
            /* Check last 4 characters for " 0 0" pattern */
            if (line[len-1] == '0' && line[len-2] == ' ' &&
                line[len-3] == '0' && line[len-4] == ' ') {
                count++;
            }
        }
    }
    fclose(fp);
    return count;
}

/*
 * Create temporary directory and file paths.
 */
static char* create_temp_dir(void) {
    char* template = strdup("/tmp/pedsim_seq_XXXXXX");
    if (!template) return NULL;
    char* dir = mkdtemp(template);
    if (!dir) {
        free(template);
        return NULL;
    }
    return template;  /* caller frees */
}

int run_pedsim_seq(pedsim_seq_params* params) {
    /* Initialize seed if not set */
    if (params->seed == 0) {
        params->seed = time(NULL);
    }

    if (params->verbose) {
        fprintf(stderr, "pedsim_seq - Integrated Pedigree-Coalescent Simulation\n");
        fprintf(stderr, "Random seed: %lu\n\n", params->seed);
    }

    /* Create temporary directory */
    char* temp_dir = create_temp_dir();
    if (!temp_dir) {
        fprintf(stderr, "Error: Failed to create temporary directory\n");
        return 1;
    }

    /* Note: coalsim limits output filename to 30 chars, use short names */
    char ped_file[512], founders_file[512], segments_file[512];
    snprintf(ped_file, sizeof(ped_file), "%s/ped.txt", temp_dir);
    snprintf(founders_file, sizeof(founders_file), "%s/f.fa", temp_dir);
    snprintf(segments_file, sizeof(segments_file), "%s/seg.txt", temp_dir);

    char cmd[2048];
    int status;

    /* Phase 1: Generate pedigree using pedsim */
    if (params->verbose) {
        fprintf(stderr, "Phase 1: Generating pedigree...\n");
        fprintf(stderr, "  Samples: %d, N: %d, Generations: %d\n",
                params->n_samples, params->N, params->k);
    }

    snprintf(cmd, sizeof(cmd),
             "./pedsim -n %d -N %d -k %d -s %lu -o %s 2>/dev/null",
             params->n_samples, params->N, params->k, params->seed, ped_file);

    status = run_command(cmd, params->verbose);
    if (status != 0) {
        fprintf(stderr, "Error: pedsim failed (exit code %d)\n", status);
        goto cleanup;
    }

    /* Count founders */
    int n_founders = count_founders(ped_file);
    if (n_founders <= 0) {
        fprintf(stderr, "Error: Failed to count founders in pedigree\n");
        status = 1;
        goto cleanup;
    }

    if (params->verbose) {
        fprintf(stderr, "  Founders: %d\n\n", n_founders);
    }

    /* Phase 2: Simulate coalescent for founders using coalsim */
    if (params->verbose) {
        fprintf(stderr, "Phase 2: Simulating coalescent...\n");
        fprintf(stderr, "  Founder haplotypes: %d\n", 2 * n_founders);
        fprintf(stderr, "  Recombination rate: %.4f\n", params->rec_rate);
        fprintf(stderr, "  Mutation rate: %.4f\n", params->mut_rate);
    }

    /* Use same seed offset for reproducibility */
    unsigned long coal_seed = params->seed + 1;

    snprintf(cmd, sizeof(cmd),
             "./coalsim -c %d -N %d -r %g -m %g -s %lu -o %s 2>/dev/null",
             2 * n_founders, params->N, params->rec_rate, params->mut_rate,
             coal_seed, founders_file);

    status = run_command(cmd, params->verbose);
    if (status != 0) {
        fprintf(stderr, "Error: coalsim failed (exit code %d)\n", status);
        goto cleanup;
    }

    if (params->verbose) {
        fprintf(stderr, "\n");
    }

    /* Phase 3: Simulate chromosome transmission using pedtrans */
    if (params->verbose) {
        fprintf(stderr, "Phase 3: Simulating chromosome transmission...\n");
        fprintf(stderr, "  Recombination rate: %.4f\n", params->pedtrans_rec);
    }

    /* Use different seed for pedtrans */
    unsigned long pedtrans_seed = params->seed + 2;

    snprintf(cmd, sizeof(cmd),
             "./pedtrans -r %g -s %lu %s -o %s 2>/dev/null",
             params->pedtrans_rec, pedtrans_seed, ped_file, segments_file);

    status = run_command(cmd, params->verbose);
    if (status != 0) {
        fprintf(stderr, "Error: pedtrans failed (exit code %d)\n", status);
        goto cleanup;
    }

    if (params->verbose) {
        fprintf(stderr, "\n");
    }

    /* Phase 4: Assemble sequences using seqassemble */
    if (params->verbose) {
        fprintf(stderr, "Phase 4: Assembling sequences...\n");
    }

    const char* output = params->output_file ? params->output_file : "-";
    const char* samples_flag = params->samples_only ? "-s" : "";

    if (strcmp(output, "-") == 0) {
        /* Output to stdout */
        snprintf(cmd, sizeof(cmd),
                 "./seqassemble -f %s %s %s %s",
                 params->format, samples_flag, segments_file, founders_file);
    } else {
        snprintf(cmd, sizeof(cmd),
                 "./seqassemble -f %s %s -o %s %s %s",
                 params->format, samples_flag, output, segments_file, founders_file);
    }

    status = run_command(cmd, params->verbose);
    if (status != 0) {
        fprintf(stderr, "Error: seqassemble failed (exit code %d)\n", status);
        goto cleanup;
    }

    if (params->verbose) {
        fprintf(stderr, "\nSimulation complete.\n");
        if (params->output_file) {
            fprintf(stderr, "Output written to %s\n", params->output_file);
        }
    }

    status = 0;

cleanup:
    /* Remove temporary files */
    unlink(ped_file);
    unlink(founders_file);
    unlink(segments_file);
    rmdir(temp_dir);
    free(temp_dir);

    return status;
}
