/*
 * fenwick.c - Fenwick tree (Binary Indexed Tree) implementation
 *
 * Provides O(log n) prefix sum queries, point updates, and weighted
 * random selection. Used to accelerate chromosome lookup during
 * recombination and mutation event placement in the coalescent simulator.
 */

#include <stdlib.h>
#include <string.h>
#include "fenwick.h"

fenwick_tree* fenwick_create(int capacity)
{
    fenwick_tree* ft = malloc(sizeof(fenwick_tree));
    if (!ft) return NULL;

    ft->tree = calloc(capacity + 1, sizeof(double));    /* 1-indexed */
    ft->weights = calloc(capacity, sizeof(double));      /* 0-indexed */
    if (!ft->tree || !ft->weights) {
        free(ft->tree);
        free(ft->weights);
        free(ft);
        return NULL;
    }

    ft->n = 0;
    ft->capacity = capacity;
    return ft;
}

void fenwick_free(fenwick_tree* ft)
{
    if (!ft) return;
    free(ft->tree);
    free(ft->weights);
    free(ft);
}

void fenwick_update(fenwick_tree* ft, int i, double delta)
{
    ft->weights[i] += delta;
    for (i += 1; i <= ft->n; i += i & (-i))
        ft->tree[i] += delta;
}

void fenwick_set(fenwick_tree* ft, int i, double new_weight)
{
    double delta = new_weight - ft->weights[i];
    ft->weights[i] = new_weight;
    for (i += 1; i <= ft->n; i += i & (-i))
        ft->tree[i] += delta;
}

double fenwick_query(fenwick_tree* ft, int i)
{
    double sum = 0;
    for (i += 1; i > 0; i -= i & (-i))
        sum += ft->tree[i];
    return sum;
}

double fenwick_total(fenwick_tree* ft)
{
    if (ft->n <= 0) return 0;
    return fenwick_query(ft, ft->n - 1);
}

double fenwick_point_value(fenwick_tree* ft, int i)
{
    return ft->weights[i];
}

int fenwick_find(fenwick_tree* ft, double target)
{
    int pos = 0;

    /* Find highest power of 2 <= n */
    int bitmask = 1;
    while (bitmask <= ft->n) bitmask <<= 1;
    bitmask >>= 1;

    while (bitmask > 0) {
        int next = pos + bitmask;
        if (next <= ft->n && ft->tree[next] <= target) {
            target -= ft->tree[next];
            pos = next;
        }
        bitmask >>= 1;
    }

    /* pos is 1-indexed count of elements whose prefix sum <= target.
     * This equals the 0-indexed chromosome index. */
    return pos;
}

void fenwick_rebuild(fenwick_tree* ft, int n)
{
    int old_n = ft->n;
    ft->n = n;

    /* Zero the full tree array to clear stale values from previous states
     * where ft->n may have been larger. */
    int clear_to = (old_n > n ? old_n : n);
    memset(ft->tree, 0, (clear_to + 1) * sizeof(double));

    /* Zero stale weights beyond new n */
    if (old_n > n)
        memset(ft->weights + n, 0, (old_n - n) * sizeof(double));

    /* Copy weights into 1-indexed tree positions */
    for (int i = 1; i <= n; i++)
        ft->tree[i] = ft->weights[i - 1];

    /* In-place linear-time construction: propagate each node to its parent */
    for (int i = 1; i <= n; i++) {
        int parent = i + (i & (-i));
        if (parent <= n)
            ft->tree[parent] += ft->tree[i];
    }
}

void fenwick_grow(fenwick_tree* ft, int new_capacity)
{
    if (new_capacity <= ft->capacity) return;

    ft->tree = realloc(ft->tree, (new_capacity + 1) * sizeof(double));
    ft->weights = realloc(ft->weights, new_capacity * sizeof(double));

    /* Zero-fill new weight slots */
    memset(ft->weights + ft->capacity, 0,
           (new_capacity - ft->capacity) * sizeof(double));

    ft->capacity = new_capacity;

    /* Rebuild tree to ensure internal consistency after realloc */
    memset(ft->tree, 0, (new_capacity + 1) * sizeof(double));
    fenwick_rebuild(ft, ft->n);
}
