/*
 * test_pedtrans.c - Unit tests for pedigree chromosome transmission simulator
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pedtrans.h"
#include "unity.h"

char version[] = "test_pedtrans";

/* Test fixtures */
pedigree* test_ped = NULL;
ped_chromosome* test_chr = NULL;
gsl_rng* test_rng = NULL;

void setUp(void) {
    /* Initialize RNG with fixed seed for reproducibility */
    test_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(test_rng, 12345);
}

void tearDown(void) {
    if (test_ped) {
        free_pedigree(test_ped);
        test_ped = NULL;
    }
    if (test_chr) {
        free_chromosome(test_chr);
        test_chr = NULL;
    }
    if (test_rng) {
        gsl_rng_free(test_rng);
        test_rng = NULL;
    }
}

/* Test pedigree creation */
void test_create_pedigree(void) {
    test_ped = create_pedigree(100);
    TEST_ASSERT_NOT_NULL(test_ped);
    TEST_ASSERT_NOT_NULL(test_ped->individuals);
    TEST_ASSERT_NOT_NULL(test_ped->name_to_id);
    TEST_ASSERT_NOT_NULL(test_ped->id_to_name);
    TEST_ASSERT_EQUAL_INT(0, test_ped->n_individuals);
    TEST_ASSERT_EQUAL_INT(0, test_ped->n_founders);
    TEST_ASSERT_EQUAL_INT(100, test_ped->capacity);
}

/* Test adding individuals to pedigree */
void test_pedigree_add_individual(void) {
    test_ped = create_pedigree(10);

    /* Add founder (parents = "0") */
    int id1 = pedigree_add_individual(test_ped, "F1", "0", "0");
    TEST_ASSERT_EQUAL_INT(0, id1);
    TEST_ASSERT_EQUAL_INT(1, test_ped->n_individuals);
    TEST_ASSERT_EQUAL_INT(1, test_ped->n_founders);
    TEST_ASSERT_EQUAL_INT(1, test_ped->individuals[id1].is_founder);
    TEST_ASSERT_EQUAL_INT(-1, test_ped->individuals[id1].father_id);
    TEST_ASSERT_EQUAL_INT(-1, test_ped->individuals[id1].mother_id);

    /* Add another founder */
    int id2 = pedigree_add_individual(test_ped, "F2", "0", "0");
    TEST_ASSERT_EQUAL_INT(1, id2);
    TEST_ASSERT_EQUAL_INT(2, test_ped->n_individuals);
    TEST_ASSERT_EQUAL_INT(2, test_ped->n_founders);

    /* Add non-founder child */
    int id3 = pedigree_add_individual(test_ped, "C1", "F1", "F2");
    TEST_ASSERT_EQUAL_INT(2, id3);
    TEST_ASSERT_EQUAL_INT(3, test_ped->n_individuals);
    TEST_ASSERT_EQUAL_INT(2, test_ped->n_founders);  /* Still 2 founders */
    TEST_ASSERT_EQUAL_INT(0, test_ped->individuals[id3].is_founder);
    TEST_ASSERT_EQUAL_INT(0, test_ped->individuals[id3].father_id);  /* F1 */
    TEST_ASSERT_EQUAL_INT(1, test_ped->individuals[id3].mother_id);  /* F2 */
}

/* Test name lookup */
void test_pedigree_name_lookup(void) {
    test_ped = create_pedigree(10);
    pedigree_add_individual(test_ped, "Alpha", "0", "0");
    pedigree_add_individual(test_ped, "Beta", "0", "0");
    pedigree_add_individual(test_ped, "Gamma", "Alpha", "Beta");

    TEST_ASSERT_EQUAL_INT(0, pedigree_get_id(test_ped, "Alpha"));
    TEST_ASSERT_EQUAL_INT(1, pedigree_get_id(test_ped, "Beta"));
    TEST_ASSERT_EQUAL_INT(2, pedigree_get_id(test_ped, "Gamma"));
    TEST_ASSERT_EQUAL_INT(-1, pedigree_get_id(test_ped, "NotFound"));

    TEST_ASSERT_EQUAL_STRING("Alpha", pedigree_get_name(test_ped, 0));
    TEST_ASSERT_EQUAL_STRING("Beta", pedigree_get_name(test_ped, 1));
    TEST_ASSERT_EQUAL_STRING("Gamma", pedigree_get_name(test_ped, 2));
    TEST_ASSERT_NULL(pedigree_get_name(test_ped, 99));
}

/* Test topological sort with simple trio */
void test_topological_sort_trio(void) {
    test_ped = create_pedigree(10);
    pedigree_add_individual(test_ped, "child", "dad", "mom");
    pedigree_add_individual(test_ped, "dad", "0", "0");
    pedigree_add_individual(test_ped, "mom", "0", "0");

    int result = topological_sort(test_ped);
    TEST_ASSERT_EQUAL_INT(0, result);  /* Success */

    /* Founders should come before child */
    int child_pos = -1, dad_pos = -1, mom_pos = -1;
    for (int i = 0; i < test_ped->n_individuals; i++) {
        int id = test_ped->topo_order[i];
        const char* name = pedigree_get_name(test_ped, id);
        if (strcmp(name, "child") == 0) child_pos = i;
        else if (strcmp(name, "dad") == 0) dad_pos = i;
        else if (strcmp(name, "mom") == 0) mom_pos = i;
    }

    TEST_ASSERT_TRUE(dad_pos < child_pos);
    TEST_ASSERT_TRUE(mom_pos < child_pos);
}

/* Test topological sort with multi-generation pedigree */
void test_topological_sort_multigenerational(void) {
    test_ped = create_pedigree(20);

    /* Generation 0: founders */
    pedigree_add_individual(test_ped, "GF1", "0", "0");
    pedigree_add_individual(test_ped, "GM1", "0", "0");
    pedigree_add_individual(test_ped, "GF2", "0", "0");
    pedigree_add_individual(test_ped, "GM2", "0", "0");

    /* Generation 1: children of founders */
    pedigree_add_individual(test_ped, "P1", "GF1", "GM1");
    pedigree_add_individual(test_ped, "P2", "GF2", "GM2");

    /* Generation 2: grandchildren */
    pedigree_add_individual(test_ped, "GC", "P1", "P2");

    int result = topological_sort(test_ped);
    TEST_ASSERT_EQUAL_INT(0, result);

    /* Check order constraints */
    int pos[7];
    const char* names[] = {"GF1", "GM1", "GF2", "GM2", "P1", "P2", "GC"};
    for (int i = 0; i < test_ped->n_individuals; i++) {
        int id = test_ped->topo_order[i];
        const char* name = pedigree_get_name(test_ped, id);
        for (int j = 0; j < 7; j++) {
            if (strcmp(name, names[j]) == 0) pos[j] = i;
        }
    }

    /* Parents before children */
    TEST_ASSERT_TRUE(pos[0] < pos[4]);  /* GF1 < P1 */
    TEST_ASSERT_TRUE(pos[1] < pos[4]);  /* GM1 < P1 */
    TEST_ASSERT_TRUE(pos[2] < pos[5]);  /* GF2 < P2 */
    TEST_ASSERT_TRUE(pos[3] < pos[5]);  /* GM2 < P2 */
    TEST_ASSERT_TRUE(pos[4] < pos[6]);  /* P1 < GC */
    TEST_ASSERT_TRUE(pos[5] < pos[6]);  /* P2 < GC */
}

/* Test chromosome creation */
void test_create_chromosome(void) {
    test_chr = create_chromosome(10);
    TEST_ASSERT_NOT_NULL(test_chr);
    TEST_ASSERT_NOT_NULL(test_chr->segments);
    TEST_ASSERT_EQUAL_INT(0, test_chr->n_segments);
    TEST_ASSERT_EQUAL_INT(10, test_chr->capacity);
}

/* Test founder chromosome creation */
void test_create_founder_chromosome(void) {
    test_chr = create_founder_chromosome(5, 1);
    TEST_ASSERT_NOT_NULL(test_chr);
    TEST_ASSERT_EQUAL_INT(1, test_chr->n_segments);
    TEST_ASSERT_EQUAL_FLOAT(0.0, test_chr->segments[0].start);
    TEST_ASSERT_EQUAL_FLOAT(1.0, test_chr->segments[0].end);
    TEST_ASSERT_EQUAL_INT(5, test_chr->segments[0].origin.founder_id);
    TEST_ASSERT_EQUAL_INT(1, test_chr->segments[0].origin.homolog);
}

/* Test adding segments */
void test_chromosome_add_segment(void) {
    test_chr = create_chromosome(10);
    founder_origin origin1 = {0, 0};
    founder_origin origin2 = {1, 1};

    chromosome_add_segment(test_chr, 0.0, 0.5, origin1);
    TEST_ASSERT_EQUAL_INT(1, test_chr->n_segments);

    chromosome_add_segment(test_chr, 0.5, 1.0, origin2);
    TEST_ASSERT_EQUAL_INT(2, test_chr->n_segments);

    TEST_ASSERT_EQUAL_FLOAT(0.0, test_chr->segments[0].start);
    TEST_ASSERT_EQUAL_FLOAT(0.5, test_chr->segments[0].end);
    TEST_ASSERT_EQUAL_INT(0, test_chr->segments[0].origin.founder_id);

    TEST_ASSERT_EQUAL_FLOAT(0.5, test_chr->segments[1].start);
    TEST_ASSERT_EQUAL_FLOAT(1.0, test_chr->segments[1].end);
    TEST_ASSERT_EQUAL_INT(1, test_chr->segments[1].origin.founder_id);
}

/* Test segment merging */
void test_merge_adjacent_segments(void) {
    test_chr = create_chromosome(10);
    founder_origin same_origin = {3, 0};
    founder_origin diff_origin = {4, 1};

    /* Add segments: same, same, different, same */
    chromosome_add_segment(test_chr, 0.0, 0.25, same_origin);
    chromosome_add_segment(test_chr, 0.25, 0.5, same_origin);  /* Should merge with prev */
    chromosome_add_segment(test_chr, 0.5, 0.75, diff_origin);  /* Different origin */
    chromosome_add_segment(test_chr, 0.75, 1.0, diff_origin);  /* Should merge with prev */

    TEST_ASSERT_EQUAL_INT(4, test_chr->n_segments);  /* Before merge */

    merge_adjacent_segments(test_chr);

    TEST_ASSERT_EQUAL_INT(2, test_chr->n_segments);  /* After merge */
    TEST_ASSERT_EQUAL_FLOAT(0.0, test_chr->segments[0].start);
    TEST_ASSERT_EQUAL_FLOAT(0.5, test_chr->segments[0].end);
    TEST_ASSERT_EQUAL_INT(3, test_chr->segments[0].origin.founder_id);
    TEST_ASSERT_EQUAL_FLOAT(0.5, test_chr->segments[1].start);
    TEST_ASSERT_EQUAL_FLOAT(1.0, test_chr->segments[1].end);
    TEST_ASSERT_EQUAL_INT(4, test_chr->segments[1].origin.founder_id);
}

/* Test copy_segments_in_range */
void test_copy_segments_in_range(void) {
    /* Source chromosome: [0, 0.3) origin 0, [0.3, 0.7) origin 1, [0.7, 1.0] origin 2 */
    ped_chromosome* src = create_chromosome(10);
    founder_origin o0 = {0, 0};
    founder_origin o1 = {1, 0};
    founder_origin o2 = {2, 0};
    chromosome_add_segment(src, 0.0, 0.3, o0);
    chromosome_add_segment(src, 0.3, 0.7, o1);
    chromosome_add_segment(src, 0.7, 1.0, o2);

    /* Copy range [0.2, 0.8) - should clip boundaries */
    test_chr = create_chromosome(10);
    copy_segments_in_range(test_chr, src, 0.2, 0.8);

    /* Should have 3 segments after clipping */
    TEST_ASSERT_EQUAL_INT(3, test_chr->n_segments);

    /* First segment: clipped [0.2, 0.3) from origin 0 */
    TEST_ASSERT_EQUAL_FLOAT(0.2, test_chr->segments[0].start);
    TEST_ASSERT_EQUAL_FLOAT(0.3, test_chr->segments[0].end);
    TEST_ASSERT_EQUAL_INT(0, test_chr->segments[0].origin.founder_id);

    /* Second segment: full [0.3, 0.7) from origin 1 */
    TEST_ASSERT_EQUAL_FLOAT(0.3, test_chr->segments[1].start);
    TEST_ASSERT_EQUAL_FLOAT(0.7, test_chr->segments[1].end);
    TEST_ASSERT_EQUAL_INT(1, test_chr->segments[1].origin.founder_id);

    /* Third segment: clipped [0.7, 0.8) from origin 2 */
    TEST_ASSERT_EQUAL_FLOAT(0.7, test_chr->segments[2].start);
    TEST_ASSERT_EQUAL_FLOAT(0.8, test_chr->segments[2].end);
    TEST_ASSERT_EQUAL_INT(2, test_chr->segments[2].origin.founder_id);

    free_chromosome(src);
}

/* Test meiosis with no recombination (rec_rate = 0) */
void test_meiosis_no_recombination(void) {
    ped_chromosome* pat = create_founder_chromosome(0, 0);  /* Founder 0, paternal */
    ped_chromosome* mat = create_founder_chromosome(0, 1);  /* Founder 0, maternal */

    /* With rec_rate = 0, should get exactly one segment from one parent */
    test_chr = meiosis(pat, mat, 0.0, test_rng);

    TEST_ASSERT_NOT_NULL(test_chr);
    TEST_ASSERT_EQUAL_INT(1, test_chr->n_segments);
    TEST_ASSERT_EQUAL_FLOAT(0.0, test_chr->segments[0].start);
    TEST_ASSERT_EQUAL_FLOAT(1.0, test_chr->segments[0].end);
    TEST_ASSERT_EQUAL_INT(0, test_chr->segments[0].origin.founder_id);
    /* Homolog should be either 0 or 1 (random) */
    TEST_ASSERT_TRUE(test_chr->segments[0].origin.homolog == 0 ||
                     test_chr->segments[0].origin.homolog == 1);

    free_chromosome(pat);
    free_chromosome(mat);
}

/* Test meiosis reproducibility with same seed */
void test_meiosis_reproducibility(void) {
    ped_chromosome* pat = create_founder_chromosome(0, 0);
    ped_chromosome* mat = create_founder_chromosome(0, 1);

    /* Run meiosis with specific seed */
    gsl_rng_set(test_rng, 99999);
    ped_chromosome* chr1 = meiosis(pat, mat, 2.0, test_rng);

    /* Reset and run again with same seed */
    gsl_rng_set(test_rng, 99999);
    ped_chromosome* chr2 = meiosis(pat, mat, 2.0, test_rng);

    /* Results should be identical */
    TEST_ASSERT_EQUAL_INT(chr1->n_segments, chr2->n_segments);
    for (int i = 0; i < chr1->n_segments; i++) {
        TEST_ASSERT_EQUAL_FLOAT(chr1->segments[i].start, chr2->segments[i].start);
        TEST_ASSERT_EQUAL_FLOAT(chr1->segments[i].end, chr2->segments[i].end);
        TEST_ASSERT_EQUAL_INT(chr1->segments[i].origin.founder_id,
                             chr2->segments[i].origin.founder_id);
        TEST_ASSERT_EQUAL_INT(chr1->segments[i].origin.homolog,
                             chr2->segments[i].origin.homolog);
    }

    free_chromosome(pat);
    free_chromosome(mat);
    free_chromosome(chr1);
    free_chromosome(chr2);
    test_chr = NULL;  /* Prevent double-free in tearDown */
}

/* Test chromosome copy */
void test_copy_chromosome(void) {
    ped_chromosome* orig = create_chromosome(10);
    founder_origin o1 = {1, 0};
    founder_origin o2 = {2, 1};
    chromosome_add_segment(orig, 0.0, 0.4, o1);
    chromosome_add_segment(orig, 0.4, 1.0, o2);

    test_chr = copy_chromosome(orig);

    TEST_ASSERT_NOT_NULL(test_chr);
    TEST_ASSERT_EQUAL_INT(orig->n_segments, test_chr->n_segments);
    TEST_ASSERT_NOT_EQUAL(orig->segments, test_chr->segments);  /* Different memory */

    for (int i = 0; i < orig->n_segments; i++) {
        TEST_ASSERT_EQUAL_FLOAT(orig->segments[i].start, test_chr->segments[i].start);
        TEST_ASSERT_EQUAL_FLOAT(orig->segments[i].end, test_chr->segments[i].end);
        TEST_ASSERT_EQUAL_INT(orig->segments[i].origin.founder_id,
                             test_chr->segments[i].origin.founder_id);
    }

    free_chromosome(orig);
}

/* Test simple trio simulation from file */
void test_simulate_trio(void) {
    /* Create temporary pedigree file */
    const char* tmp_file = "/tmp/test_trio.ped";
    FILE* fp = fopen(tmp_file, "w");
    TEST_ASSERT_NOT_NULL(fp);
    fprintf(fp, "dad 0 0\n");
    fprintf(fp, "mom 0 0\n");
    fprintf(fp, "child dad mom\n");
    fclose(fp);

    /* Run simulation */
    ped_simulation* sim = simulate_pedigree(tmp_file, 1.0, 12345);
    TEST_ASSERT_NOT_NULL(sim);
    TEST_ASSERT_NOT_NULL(sim->ped);
    TEST_ASSERT_NOT_NULL(sim->chromosomes);
    TEST_ASSERT_EQUAL_INT(3, sim->ped->n_individuals);
    TEST_ASSERT_EQUAL_INT(2, sim->ped->n_founders);

    /* Check founders have single-segment chromosomes */
    int dad_id = pedigree_get_id(sim->ped, "dad");
    int mom_id = pedigree_get_id(sim->ped, "mom");
    int child_id = pedigree_get_id(sim->ped, "child");

    TEST_ASSERT_EQUAL_INT(1, sim->chromosomes[dad_id].paternal->n_segments);
    TEST_ASSERT_EQUAL_INT(1, sim->chromosomes[dad_id].maternal->n_segments);
    TEST_ASSERT_EQUAL_INT(1, sim->chromosomes[mom_id].paternal->n_segments);
    TEST_ASSERT_EQUAL_INT(1, sim->chromosomes[mom_id].maternal->n_segments);

    /* Child's chromosomes should have segments traceable to founders */
    TEST_ASSERT_NOT_NULL(sim->chromosomes[child_id].paternal);
    TEST_ASSERT_NOT_NULL(sim->chromosomes[child_id].maternal);
    TEST_ASSERT_TRUE(sim->chromosomes[child_id].paternal->n_segments >= 1);
    TEST_ASSERT_TRUE(sim->chromosomes[child_id].maternal->n_segments >= 1);

    /* Verify child's paternal chromosome comes from dad's founder origins */
    for (int i = 0; i < sim->chromosomes[child_id].paternal->n_segments; i++) {
        int fid = sim->chromosomes[child_id].paternal->segments[i].origin.founder_id;
        TEST_ASSERT_EQUAL_INT(dad_id, fid);  /* Should trace to dad */
    }

    /* Verify child's maternal chromosome comes from mom's founder origins */
    for (int i = 0; i < sim->chromosomes[child_id].maternal->n_segments; i++) {
        int fid = sim->chromosomes[child_id].maternal->segments[i].origin.founder_id;
        TEST_ASSERT_EQUAL_INT(mom_id, fid);  /* Should trace to mom */
    }

    free_simulation(sim);
    remove(tmp_file);
}

/* Test origin_equal utility */
void test_origin_equal(void) {
    founder_origin a = {1, 0};
    founder_origin b = {1, 0};
    founder_origin c = {1, 1};
    founder_origin d = {2, 0};

    TEST_ASSERT_TRUE(origin_equal(a, b));
    TEST_ASSERT_FALSE(origin_equal(a, c));
    TEST_ASSERT_FALSE(origin_equal(a, d));
}

int main(void) {
    UNITY_BEGIN();

    /* Pedigree tests */
    RUN_TEST(test_create_pedigree);
    RUN_TEST(test_pedigree_add_individual);
    RUN_TEST(test_pedigree_name_lookup);
    RUN_TEST(test_topological_sort_trio);
    RUN_TEST(test_topological_sort_multigenerational);

    /* Chromosome/segment tests */
    RUN_TEST(test_create_chromosome);
    RUN_TEST(test_create_founder_chromosome);
    RUN_TEST(test_chromosome_add_segment);
    RUN_TEST(test_merge_adjacent_segments);
    RUN_TEST(test_copy_segments_in_range);
    RUN_TEST(test_copy_chromosome);

    /* Meiosis tests */
    RUN_TEST(test_meiosis_no_recombination);
    RUN_TEST(test_meiosis_reproducibility);

    /* Simulation tests */
    RUN_TEST(test_simulate_trio);

    /* Utility tests */
    RUN_TEST(test_origin_equal);

    return UNITY_END();
}
