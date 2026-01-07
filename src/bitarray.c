#include "bitarray.h"

/*
 * Bitarray implementation for coalescent ancestry tracking.
 */

/*
 * Memory pool for bitarray allocation.
 * Uses a free list to recycle bitarrays instead of malloc/free.
 */
typedef struct pool_entry {
    bitarray ba;
    struct pool_entry *next;
} pool_entry;

static struct {
    pool_entry *free_list;    /* linked list of free bitarrays */
    int pool_nbits;           /* size of bitarrays in pool */
    int pool_nwords;          /* words per bitarray */
    int pool_size;            /* current number of entries in free list */
    int total_allocs;         /* total allocations requested */
    int pool_hits;            /* allocations served from pool */
    int initialized;          /* pool initialized flag */
} g_pool = {NULL, 0, 0, 0, 0, 0, 0};

void bitarray_pool_init(int nbits) {
    g_pool.pool_nbits = nbits;
    g_pool.pool_nwords = bitarray_nwords(nbits);
    g_pool.free_list = NULL;
    g_pool.pool_size = 0;
    g_pool.total_allocs = 0;
    g_pool.pool_hits = 0;
    g_pool.initialized = 1;
}

void bitarray_pool_destroy(void) {
    pool_entry *curr = g_pool.free_list;
    while (curr) {
        pool_entry *next = curr->next;
        free(curr->ba.bits);
        free(curr);
        curr = next;
    }
    g_pool.free_list = NULL;
    g_pool.pool_size = 0;
    g_pool.initialized = 0;
}

void bitarray_pool_stats(int *pool_size, int *total_allocs, int *pool_hits) {
    if (pool_size) *pool_size = g_pool.pool_size;
    if (total_allocs) *total_allocs = g_pool.total_allocs;
    if (pool_hits) *pool_hits = g_pool.pool_hits;
}

bitarray* bitarray_create(int nbits) {
    g_pool.total_allocs++;

    /* Try to get from pool if initialized and size matches */
    if (g_pool.initialized && nbits == g_pool.pool_nbits && g_pool.free_list) {
        pool_entry *entry = g_pool.free_list;
        g_pool.free_list = entry->next;
        g_pool.pool_size--;
        g_pool.pool_hits++;

        /* Clear the bits array and return */
        bitarray *ba = &entry->ba;
        memset(ba->bits, 0, ba->nwords * sizeof(unsigned long));
        ba->is_zero = 1;  /* freshly cleared */
        return ba;
    }

    /* Fall back to malloc */
    pool_entry *entry = malloc(sizeof(pool_entry));
    if (!entry) return NULL;

    entry->ba.nbits = nbits;
    entry->ba.nwords = bitarray_nwords(nbits);
    entry->ba.bits = calloc(entry->ba.nwords, sizeof(unsigned long));
    if (!entry->ba.bits) {
        free(entry);
        return NULL;
    }
    entry->ba.is_zero = 1;  /* calloc zeros the bits */
    entry->next = NULL;
    return &entry->ba;
}

void bitarray_free(bitarray *ba) {
    if (!ba) return;

    /* Get the containing pool_entry */
    pool_entry *entry = (pool_entry *)ba;

    /* Return to pool if initialized and size matches */
    if (g_pool.initialized && ba->nbits == g_pool.pool_nbits) {
        entry->next = g_pool.free_list;
        g_pool.free_list = entry;
        g_pool.pool_size++;
        return;
    }

    /* Fall back to free */
    free(ba->bits);
    free(entry);
}

bitarray* bitarray_copy(const bitarray *ba) {
    if (!ba) return NULL;

    bitarray *copy = bitarray_create(ba->nbits);
    if (!copy) return NULL;

    memcpy(copy->bits, ba->bits, ba->nwords * sizeof(unsigned long));
    copy->is_zero = ba->is_zero;  /* copy the cached flag */
    return copy;
}

void bitarray_clear_all(bitarray *ba) {
    memset(ba->bits, 0, ba->nwords * sizeof(unsigned long));
    ba->is_zero = 1;  /* now all zeros */
}

void bitarray_union(bitarray *dst, const bitarray *src) {
    for (int i = 0; i < dst->nwords; i++) {
        dst->bits[i] |= src->bits[i];
    }
    /* Result is zero only if both were zero */
    if (dst->is_zero == 1 && src->is_zero == 1)
        dst->is_zero = 1;
    else if (src->is_zero == 0)
        dst->is_zero = 0;  /* src had bits, so dst now has bits */
    /* else dst->is_zero unchanged (src was zero) */
}

void bitarray_union_into(bitarray *dst, const bitarray *a, const bitarray *b) {
    for (int i = 0; i < dst->nwords; i++) {
        dst->bits[i] = a->bits[i] | b->bits[i];
    }
    /* Result is zero only if both inputs are zero */
    dst->is_zero = (a->is_zero == 1 && b->is_zero == 1) ? 1 : 0;
}

void bitarray_intersect(bitarray *dst, const bitarray *src) {
    for (int i = 0; i < dst->nwords; i++) {
        dst->bits[i] &= src->bits[i];
    }
    /* Either input being zero means result is zero */
    if (dst->is_zero == 1 || src->is_zero == 1)
        dst->is_zero = 1;
    else
        dst->is_zero = -1;  /* unknown, need to recompute if checked */
}

void bitarray_complement(bitarray *dst, const bitarray *src) {
    for (int i = 0; i < dst->nwords; i++) {
        dst->bits[i] = ~src->bits[i];
    }
    /* Clear bits beyond nbits in the last word */
    int extra_bits = dst->nbits % BITS_PER_WORD;
    if (extra_bits > 0) {
        unsigned long mask = (1UL << extra_bits) - 1;
        dst->bits[dst->nwords - 1] &= mask;
    }
    /* Complement of zero is non-zero (if nbits > 0) */
    dst->is_zero = (src->is_zero == 1 && dst->nbits > 0) ? 0 : -1;
}

int bitarray_equal(const bitarray *a, const bitarray *b) {
    if (a->nwords != b->nwords) return 0;
    for (int i = 0; i < a->nwords; i++) {
        if (a->bits[i] != b->bits[i]) return 0;
    }
    return 1;
}

int bitarray_is_zero(const bitarray *ba) {
    /* Use cached value if available */
    if (ba->is_zero >= 0) return ba->is_zero;

    /* Recompute if unknown (-1) */
    for (int i = 0; i < ba->nwords; i++) {
        if (ba->bits[i] != 0) {
            ((bitarray*)ba)->is_zero = 0;  /* cache the result */
            return 0;
        }
    }
    ((bitarray*)ba)->is_zero = 1;  /* cache the result */
    return 1;
}

int bitarray_is_singleton(const bitarray *ba) {
    int count = 0;
    for (int i = 0; i < ba->nwords; i++) {
        count += __builtin_popcountl(ba->bits[i]);
        if (count > 1) return 0;
    }
    return count == 1;
}

int bitarray_popcount(const bitarray *ba) {
    int count = 0;
    for (int i = 0; i < ba->nwords; i++) {
        count += __builtin_popcountl(ba->bits[i]);
    }
    return count;
}

int bitarray_intersects(const bitarray *a, const bitarray *b) {
    for (int i = 0; i < a->nwords; i++) {
        if (a->bits[i] & b->bits[i]) return 1;
    }
    return 0;
}

void bitarray_copy_to(bitarray *dst, const bitarray *src) {
    memcpy(dst->bits, src->bits, dst->nwords * sizeof(unsigned long));
    dst->is_zero = src->is_zero;  /* copy the cached flag */
}

void bitarray_print(const bitarray *ba, FILE *out) {
    for (int i = 0; i < ba->nbits; i++) {
        fprintf(out, "%c", bitarray_test(ba, i) ? '1' : '0');
    }
}

void bitarray_to_string(const bitarray *ba, char *buf) {
    for (int i = 0; i < ba->nbits; i++) {
        buf[i] = bitarray_test(ba, i) ? '1' : '0';
    }
    buf[ba->nbits] = '\0';
}

bitarray* bitarray_singleton(int nbits, int i) {
    bitarray *ba = bitarray_create(nbits);
    if (ba && i >= 0 && i < nbits) {
        bitarray_set(ba, i);  /* bitarray_set already sets is_zero = 0 */
    }
    return ba;
}

bitarray* bitarray_full(int nbits) {
    bitarray *ba = bitarray_create(nbits);
    if (!ba) return NULL;

    /* Set all bits in complete words */
    int full_words = nbits / BITS_PER_WORD;
    for (int i = 0; i < full_words; i++) {
        ba->bits[i] = ~0UL;
    }

    /* Set remaining bits in last word */
    int extra_bits = nbits % BITS_PER_WORD;
    if (extra_bits > 0) {
        ba->bits[ba->nwords - 1] = (1UL << extra_bits) - 1;
    }

    ba->is_zero = (nbits > 0) ? 0 : 1;  /* full bitarray is non-zero if nbits > 0 */
    return ba;
}
