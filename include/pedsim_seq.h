/*
 * pedsim_seq.h - Integrated Pedigree-Coalescent Sequence Simulation
 *
 * Combines pedsim, coalsim, pedtrans, and seqassemble into a single tool
 * for simulating genetic sequences through pedigrees.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#ifndef PEDSIM_SEQ_H
#define PEDSIM_SEQ_H

#include <stdio.h>

/*
 * Parameters for the integrated simulation.
 */
typedef struct {
    /* Pedigree generation (pedsim) */
    int n_samples;          /* number of sample individuals */
    int N;                  /* effective population size */
    int k;                  /* number of generations */

    /* Coalescent simulation (coalsim) */
    double rec_rate;        /* recombination rate per chromosome */
    double mut_rate;        /* mutation rate per chromosome */
    long seq_length;        /* sequence length in bases (0 = default) */

    /* Chromosome transmission (pedtrans) */
    double pedtrans_rec;    /* recombination rate for pedtrans */

    /* Output options */
    const char* output_file;
    const char* format;     /* "fasta" or "vcf" */
    int samples_only;       /* exclude founders from output */
    int verbose;

    /* Random seed */
    unsigned long seed;
} pedsim_seq_params;

/*
 * Initialize default parameters.
 */
void init_pedsim_seq_params(pedsim_seq_params* params);

/*
 * Run the integrated simulation.
 * Returns 0 on success, non-zero on error.
 */
int run_pedsim_seq(pedsim_seq_params* params);

#endif /* PEDSIM_SEQ_H */
