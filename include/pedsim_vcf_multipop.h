/*
 * pedsim_vcf_multipop.h - Multi-population VCF Simulation Pipeline
 *
 * Integrated pipeline for multi-population pedigree-coalescent VCF simulation.
 * 1. pedsim_multipop - generate multi-population pedigree with migration
 * 2. coalsim - simulate founder haplotypes (shared coalescent)
 * 3. pedtrans - simulate chromosome transmission through merged pedigree
 * 4. vcfassemble - assemble sample VCF from founder variants
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#ifndef PEDSIM_VCF_MULTIPOP_H
#define PEDSIM_VCF_MULTIPOP_H

/* Maximum number of populations */
#define MAX_POPULATIONS 100

/* Population-specific parameters */
typedef struct {
    int pop_size;           /* N_j: population size */
    int sample_size;        /* n_j: number of sampled individuals */
} pop_vcf_params;

/* Integrated pipeline parameters */
typedef struct {
    /* Population structure */
    int n_populations;              /* J: number of populations */
    pop_vcf_params* populations;    /* Parameters for each population */
    double** migration;             /* J x J migration matrix */
    double island_m;                /* Island model migration rate (-1 if not used) */

    /* Shared pedigree parameters */
    int n_generations;              /* k: generations to trace back */
    int n_chromosomes;              /* c: number of chromosomes to simulate */

    /* Coalescent parameters */
    double coal_pop_size;           /* Ne: effective population size */
    double rec_rate;                /* r: recombination rate in cM */
    double mut_rate;                /* m: mutation rate */

    /* Random seeds and output */
    unsigned int seed;
    char* output_vcf;
    int samples_only;               /* Exclude founders from output */
    int keep_temp;                  /* Keep temporary files */
    int verbose;                    /* Verbose output */
} multipop_vcf_pipeline_params;

/*
 * Run the integrated multi-population pipeline:
 * 1. pedsim_multipop - generate pedigree with migration
 * 2. coalsim - simulate founder haplotypes (VCF output)
 * 3. pedtrans - simulate chromosome transmission
 * 4. vcfassemble - assemble sample VCF
 *
 * Returns 0 on success, non-zero on error
 */
int run_multipop_vcf_pipeline(multipop_vcf_pipeline_params* params);

#endif /* PEDSIM_VCF_MULTIPOP_H */
