#ifndef FENWICK_H
#define FENWICK_H

/*
 * Fenwick tree (Binary Indexed Tree) for O(log n) prefix sum queries
 * and point updates. Used for weighted random selection of chromosomes
 * during recombination and mutation event placement.
 *
 * Stores a parallel weights[] array for O(1) point value queries,
 * which is needed for the swap-and-pop deletion pattern.
 */

typedef struct {
    double* tree;     /* 1-indexed Fenwick array, size capacity+1 */
    double* weights;  /* 0-indexed raw weights, size capacity */
    int n;            /* current number of active elements */
    int capacity;     /* allocated capacity */
} fenwick_tree;

/* Create a Fenwick tree with given capacity. All weights initialized to 0. */
fenwick_tree* fenwick_create(int capacity);

/* Free a Fenwick tree */
void fenwick_free(fenwick_tree* ft);

/* Add delta to element at index i (0-indexed). O(log n). */
void fenwick_update(fenwick_tree* ft, int i, double delta);

/* Set element at index i to new_weight (computes delta internally). O(log n). */
void fenwick_set(fenwick_tree* ft, int i, double new_weight);

/* Prefix sum of elements [0, i] (0-indexed). O(log n). */
double fenwick_query(fenwick_tree* ft, int i);

/* Total sum of all n elements. O(log n). */
double fenwick_total(fenwick_tree* ft);

/* Get the weight at a single index. O(1) via cached weights array. */
double fenwick_point_value(fenwick_tree* ft, int i);

/*
 * Find the 0-indexed chromosome where a cumulative weight target falls.
 * Returns the largest index k such that prefix_sum(0..k-1) <= target.
 * Equivalently: the chromosome that contains position 'target' in the
 * cumulative weight space. Single O(log n) top-down walk.
 */
int fenwick_find(fenwick_tree* ft, double target);

/*
 * Rebuild the Fenwick tree from weights[0..n-1]. O(n).
 * Call this after bulk weight changes (e.g., after coalescence).
 */
void fenwick_rebuild(fenwick_tree* ft, int n);

/*
 * Grow the Fenwick tree to new_capacity. Rebuilds internal structure.
 * No-op if new_capacity <= current capacity.
 */
void fenwick_grow(fenwick_tree* ft, int new_capacity);

#endif /* FENWICK_H */
