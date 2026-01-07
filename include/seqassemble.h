/*
 * seqassemble.h - Sequence Assembly from Pedigree Segments
 *
 * Combines founder sequences (from coalsim) with segment transmission data
 * (from pedtrans) to assemble final sample sequences.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#ifndef SEQASSEMBLE_H
#define SEQASSEMBLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Maximum line length for parsing */
#define MAX_LINE_LENGTH 4096
#define MAX_NAME_LENGTH 256

/*
 * Founder haplotype sequence.
 * Each founder has two haplotypes (paternal and maternal).
 */
typedef struct {
    char* name;           /* founder name from pedigree */
    char* pat_seq;        /* paternal haplotype sequence */
    char* mat_seq;        /* maternal haplotype sequence */
    long seq_length;      /* sequence length in bases */
} founder_haplotypes;

/*
 * Collection of all founder sequences.
 */
typedef struct {
    founder_haplotypes* founders;
    int n_founders;
    int capacity;
    long seq_length;      /* common sequence length */
} founder_sequences;

/*
 * Segment from pedtrans output.
 */
typedef struct {
    double start;         /* start position (0-1) */
    double end;           /* end position (0-1) */
    char founder[MAX_NAME_LENGTH];  /* founder name */
    int homolog;          /* 0=paternal, 1=maternal */
} parsed_segment;

/*
 * Chromosome (list of segments).
 */
typedef struct {
    parsed_segment* segments;
    int n_segments;
    int capacity;
} parsed_chromosome;

/*
 * Individual with two chromosomes.
 */
typedef struct {
    char name[MAX_NAME_LENGTH];
    parsed_chromosome paternal;
    parsed_chromosome maternal;
    int is_founder;
} parsed_individual;

/*
 * Parsed pedtrans output.
 */
typedef struct {
    parsed_individual* individuals;
    int n_individuals;
    int capacity;
    char** founder_names;  /* array of founder names */
    int n_founders;
} pedtrans_data;

/*
 * ============================================================================
 * FASTA Index for Streaming (Memory-Efficient Mode)
 * ============================================================================
 */

/*
 * Index entry for one haplotype sequence in FASTA file.
 */
typedef struct {
    char founder[MAX_NAME_LENGTH];  /* founder name */
    int homolog;                     /* 0=paternal, 1=maternal */
    long file_offset;               /* byte offset to start of sequence data */
    long seq_length;                /* sequence length in bases */
    int line_width;                 /* characters per line (for seeking) */
} fasta_index_entry;

/*
 * Index of all sequences in a FASTA file.
 */
typedef struct {
    fasta_index_entry* entries;
    int n_entries;
    int capacity;
    FILE* fp;                       /* open file handle for reading */
    char* filename;                 /* filename for reopening if needed */
    long seq_length;                /* common sequence length */
} fasta_index;

/*
 * Create FASTA index by scanning file headers.
 * Maps founder names to file offsets for random access.
 */
fasta_index* create_fasta_index(const char* filename,
                                 const char** founder_names, int n_founders);

/*
 * Free FASTA index and close file.
 */
void free_fasta_index(fasta_index* idx);

/*
 * Read a range of bases from indexed FASTA file.
 * Returns number of bases read, or -1 on error.
 * Caller provides buffer of at least (end - start) bytes.
 */
int read_sequence_range(fasta_index* idx, const char* founder,
                        int homolog, long start, long end, char* buffer);

/*
 * ============================================================================
 * Founder Sequence Functions (In-Memory Mode)
 * ============================================================================
 */

/*
 * Create an empty founder sequences container.
 */
founder_sequences* create_founder_sequences(void);

/*
 * Free founder sequences.
 */
void free_founder_sequences(founder_sequences* fs);

/*
 * Read founder sequences from FASTA file.
 * Expects sequences labeled as FounderName:pat and FounderName:mat
 * or as sample0, sample1, ... (mapping to founders in order)
 *
 * If use_sample_names is 1, expects sample0, sample1 format where
 * even indices are paternal, odd indices are maternal.
 */
int read_founder_fasta(founder_sequences* fs, const char* filename,
                       const char** founder_names, int n_founders);

/*
 * Get sequence for a founder haplotype by name.
 * homolog: 0=paternal, 1=maternal
 * Returns pointer to sequence or NULL if not found.
 */
const char* get_founder_sequence(founder_sequences* fs,
                                  const char* founder_name, int homolog);

/*
 * ============================================================================
 * Pedtrans Data Functions
 * ============================================================================
 */

/*
 * Create empty pedtrans data container.
 */
pedtrans_data* create_pedtrans_data(void);

/*
 * Free pedtrans data.
 */
void free_pedtrans_data(pedtrans_data* pd);

/*
 * Parse pedtrans output file.
 */
int parse_pedtrans_output(pedtrans_data* pd, const char* filename);

/*
 * ============================================================================
 * Sequence Assembly Functions
 * ============================================================================
 */

/*
 * Assemble sequence for one chromosome from segments.
 * Returns allocated sequence string (caller must free).
 */
char* assemble_chromosome_sequence(parsed_chromosome* chr,
                                    founder_sequences* fs,
                                    long seq_length);

/*
 * Assemble all sequences for all individuals.
 * Writes output in FASTA format.
 */
int assemble_all_sequences(pedtrans_data* pd, founder_sequences* fs,
                           FILE* out, int samples_only);

/*
 * Streaming version - assembles directly from indexed FASTA file.
 * Uses minimal memory by reading sequences on demand.
 */
int assemble_all_sequences_streaming(pedtrans_data* pd, fasta_index* idx,
                                      FILE* out, int samples_only);

/*
 * ============================================================================
 * VCF Output Functions
 * ============================================================================
 */

/*
 * Write assembled sequences in VCF format.
 * Only outputs polymorphic sites.
 */
int write_vcf_output(pedtrans_data* pd, founder_sequences* fs,
                     FILE* out, int samples_only);

#endif /* SEQASSEMBLE_H */
