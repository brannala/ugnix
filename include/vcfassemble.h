#ifndef VCFASSEMBLE_H
#define VCFASSEMBLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NAME_LENGTH 256
#define MAX_LINE_LENGTH 65536
#define INITIAL_CAPACITY 1024

/* A single mutation from founder VCF */
typedef struct {
    int chrom_id;           /* chromosome ID (1-based) */
    long position;          /* 1-based position */
    char ref;               /* reference allele */
    char alt;               /* alternate allele */
    int* founder_geno;      /* genotype for each founder haplotype (0 or 1) */
    int n_founders;         /* number of founder haplotypes */
} vcf_mutation;

/* Collection of mutations from founder VCF */
typedef struct {
    vcf_mutation* mutations;
    int n_mutations;
    int capacity;
    int n_chromosomes;      /* number of chromosomes */
    long* chrom_lengths;    /* length of each chromosome (1-indexed) */
    char** founder_names;   /* sample names from VCF header */
    int n_founder_samples;  /* number of founder samples (haplotypes = 2x this) */
} founder_vcf;

/* Segment from pedtrans output */
typedef struct {
    double start;           /* start position (0-1 scale) */
    double end;             /* end position (0-1 scale) */
    char founder[MAX_NAME_LENGTH];  /* founder:homolog string (e.g., "I35:pat") */
    int homolog;            /* VCF haplotype index (computed from founder:homolog) */
} segment;

/* Sample's chromosome inheritance */
typedef struct {
    char name[MAX_NAME_LENGTH];
    int chrom_id;           /* chromosome ID (1-based) */
    int homolog;            /* which homolog (0=paternal, 1=maternal) */
    segment* segments;
    int n_segments;
    int capacity;
} sample_chrom;

/* Forward declaration for internal founder map */
struct founder_map_internal;

/* All sample chromosomes from pedtrans */
typedef struct {
    sample_chrom* chroms;
    int n_chroms;
    int capacity;
    int n_chromosomes;      /* number of distinct chromosomes */
    char** sample_names;    /* unique sample names */
    int n_samples;
    void* founder_map;      /* internal founder name -> VCF index mapping */
} pedtrans_segments;

/* Read founder VCF file */
founder_vcf* read_founder_vcf(const char* filename);

/* Free founder VCF */
void free_founder_vcf(founder_vcf* fvcf);

/* Read pedtrans segment output */
pedtrans_segments* read_pedtrans_segments(const char* filename);

/* Free pedtrans segments */
void free_pedtrans_segments(pedtrans_segments* pts);

/* Get founder haplotype index from name and homolog */
int get_founder_haplotype_idx(founder_vcf* fvcf, const char* founder, int homolog);

/* Assemble sample VCF from founder VCF and pedtrans segments */
int assemble_sample_vcf(founder_vcf* fvcf, pedtrans_segments* pts,
                        FILE* out, int samples_only);

#endif /* VCFASSEMBLE_H */
