#ifndef BITARRAY_H
#define BITARRAY_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Bitarray: flexible bit array for tracking ancestry in coalescent simulations.
 * Supports arbitrary numbers of samples (not limited to 32 or 64).
 *
 * Each bitarray stores bits in an array of unsigned long (64-bit) words.
 * Bit i is stored in word[i/64], bit position (i % 64).
 */

#define BITS_PER_WORD 64

typedef struct {
    unsigned long *bits;  /* array of 64-bit words */
    int nwords;           /* number of words allocated */
    int nbits;            /* number of bits (sample size) */
    int is_zero;          /* cached: 1 if all bits are zero, 0 otherwise */
} bitarray;

/* Calculate number of words needed for n bits */
static inline int bitarray_nwords(int nbits) {
    return (nbits + BITS_PER_WORD - 1) / BITS_PER_WORD;
}

/* Create a new bitarray with all bits cleared */
bitarray* bitarray_create(int nbits);

/* Free a bitarray */
void bitarray_free(bitarray *ba);

/* Create a deep copy of a bitarray */
bitarray* bitarray_copy(const bitarray *ba);

/* Set bit at position i */
static inline void bitarray_set(bitarray *ba, int i) {
    ba->bits[i / BITS_PER_WORD] |= (1UL << (i % BITS_PER_WORD));
    ba->is_zero = 0;  /* now has at least one bit set */
}

/* Clear bit at position i - may need to recompute is_zero */
static inline void bitarray_clear(bitarray *ba, int i) {
    ba->bits[i / BITS_PER_WORD] &= ~(1UL << (i % BITS_PER_WORD));
    /* is_zero status unknown - mark as needing recompute (-1) */
    ba->is_zero = -1;
}

/* Test if bit at position i is set */
static inline int bitarray_test(const bitarray *ba, int i) {
    return (ba->bits[i / BITS_PER_WORD] & (1UL << (i % BITS_PER_WORD))) != 0;
}

/* Set all bits to zero */
void bitarray_clear_all(bitarray *ba);

/* dst = dst | src (union/OR) */
void bitarray_union(bitarray *dst, const bitarray *src);

/* dst = a | b (creates union without modifying a or b) */
void bitarray_union_into(bitarray *dst, const bitarray *a, const bitarray *b);

/* dst = dst & src (intersection/AND) */
void bitarray_intersect(bitarray *dst, const bitarray *src);

/* dst = ~src (complement, only within nbits) */
void bitarray_complement(bitarray *dst, const bitarray *src);

/* Test if two bitarrays are equal */
int bitarray_equal(const bitarray *a, const bitarray *b);

/* Test if bitarray is all zeros */
int bitarray_is_zero(const bitarray *ba);

/* Test if exactly one bit is set (singleton) */
int bitarray_is_singleton(const bitarray *ba);

/* Count number of bits set */
int bitarray_popcount(const bitarray *ba);

/* Test if (a & b) is non-zero */
int bitarray_intersects(const bitarray *a, const bitarray *b);

/* Copy src to dst (dst must already be allocated with same size) */
void bitarray_copy_to(bitarray *dst, const bitarray *src);

/* Print bitarray as string of 0s and 1s (for nbits bits) */
void bitarray_print(const bitarray *ba, FILE *out);

/* Print bitarray to string buffer (must have space for nbits+1 chars) */
void bitarray_to_string(const bitarray *ba, char *buf);

/* Create bitarray with single bit i set */
bitarray* bitarray_singleton(int nbits, int i);

/* Create bitarray with all bits 0..nbits-1 set (full MRCA) */
bitarray* bitarray_full(int nbits);

/*
 * Memory pool for bitarray allocation.
 * Significantly reduces malloc/free overhead in coalescent simulations
 * where millions of bitarrays of the same size are created and destroyed.
 */

/* Initialize the pool for bitarrays of given size. Call once before simulation. */
void bitarray_pool_init(int nbits);

/* Destroy the pool and free all memory. Call after simulation. */
void bitarray_pool_destroy(void);

/* Get pool statistics (for debugging/profiling) */
void bitarray_pool_stats(int *pool_size, int *total_allocs, int *pool_hits);

#endif /* BITARRAY_H */
