#ifndef PEDSIM_VCF_H
#define PEDSIM_VCF_H

/* Integrated pipeline parameters */
typedef struct {
    /* Pedigree parameters */
    int n_generations;      /* number of generations to trace back */
    int sample_size;        /* number of individuals at present */
    double ped_pop_size;    /* population size for pedigree simulation */

    /* Coalescent parameters */
    double coal_pop_size;   /* effective population size for coalescent */
    double rec_rate;        /* recombination rate in cM */
    double mut_rate;        /* mutation rate */
    unsigned int seed;

    /* Output control */
    char* output_vcf;
    int samples_only;
    int keep_temp;
    int verbose;
} pedsim_vcf_params;

/* Run the integrated pipeline:
 * 1. pedsim - generate pedigree
 * 2. coalsim - simulate founder haplotypes (VCF output)
 * 3. pedtrans - simulate chromosome transmission
 * 4. vcfassemble - assemble sample VCF
 */
int run_pedsim_vcf_pipeline(pedsim_vcf_params* params);

#endif /* PEDSIM_VCF_H */
