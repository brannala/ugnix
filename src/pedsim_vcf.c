#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include "pedsim_vcf.h"

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

/* Count founders in pedigree file */
static int count_founders(const char* ped_file)
{
    FILE* fp = fopen(ped_file, "r");
    if (!fp) return -1;

    int count = 0;
    char line[1024];

    /* Skip header */
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }

    while (fgets(line, sizeof(line), fp)) {
        char id[256], father[256], mother[256];
        if (sscanf(line, "%s %s %s", id, father, mother) == 3) {
            if (strcmp(father, "0") == 0 && strcmp(mother, "0") == 0) {
                count++;
            }
        }
    }

    fclose(fp);
    return count;
}

/* Run the integrated pipeline */
int run_pedsim_vcf_pipeline(pedsim_vcf_params* params)
{
    char cmd[4096];
    char* temp_dir = NULL;
    int ret = 0;

    /* Create temporary directory */
    char template[] = "/tmp/pedsim_vcf_XXXXXX";
    temp_dir = mkdtemp(template);
    if (!temp_dir) {
        fprintf(stderr, "Error: failed to create temp directory\n");
        return 1;
    }

    char ped_file[512], vcf_file[512], seg_file[512];
    snprintf(ped_file, sizeof(ped_file), "%s/ped.txt", temp_dir);
    snprintf(vcf_file, sizeof(vcf_file), "%s/f.vcf", temp_dir);
    snprintf(seg_file, sizeof(seg_file), "%s/seg.txt", temp_dir);

    /* Step 1: Run pedsim to generate pedigree */
    fprintf(stderr, "\n=== Step 1: Generating pedigree ===\n");
    snprintf(cmd, sizeof(cmd),
             "./pedsim -k %d -n %d -N %.0f -s %u -o %s",
             params->n_generations, params->sample_size, params->ped_pop_size,
             params->seed, ped_file);
    ret = run_command(cmd, params->verbose);
    if (ret != 0) {
        fprintf(stderr, "Error: pedsim failed\n");
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
    fprintf(stderr, "Pedigree has %d founders (%d haplotypes)\n",
            n_founders, n_haplotypes);

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

cleanup:
    /* Clean up temporary files */
    if (!params->keep_temp && temp_dir) {
        snprintf(cmd, sizeof(cmd), "rm -rf %s", temp_dir);
        system(cmd);
    } else if (params->keep_temp) {
        fprintf(stderr, "Temp files kept in: %s\n", temp_dir);
    }

    return ret;
}
