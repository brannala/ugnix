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
    double div_time;                /* T: star-tree divergence time (generations) of the
                                       populations founders are drawn from; < 0 = off
                                       (all founders from one panmictic coalescent) */

    /* Migration before sampling and output thinning */
    int migrate_before_sampling;    /* 1: sampled individuals can be first-gen migrants */
    int snps_per_chrom;             /* keep at most this many SNPs per chromosome
                                       (polymorphic in the output samples); 0 = all */

    /* Random seeds and output */
    unsigned int seed;
    char* output_vcf;
    int samples_only;               /* Exclude founders from output */
    int keep_temp;                  /* Keep temporary files */
    int verbose;                    /* Verbose output */
} multipop_vcf_pipeline_params;

/*
 * Run the integrated multi-population pipeline:
 * 1. pedsim_multipop - generate pedigree with migration (+ truth file)
 * 2. pedtrans - simulate chromosome transmission
 * 3. founder haplotypes: coalsim (one panmictic population), or coalsim_msc
 *    on a star tree with divergence time T (founders drawn from their own
 *    population), one run per chromosome
 * 4. vcfassemble - assemble sample VCF
 * 5. optional thinning to snps_per_chrom SNPs per chromosome
 * Also writes <output>.meta.txt (sample, population; for BA3 -M) and
 * <output>.truth.tsv (migrant ancestry of the sampled individuals).
 *
 * Returns 0 on success, non-zero on error
 */
int run_multipop_vcf_pipeline(multipop_vcf_pipeline_params* params);

#endif /* PEDSIM_VCF_MULTIPOP_H */
