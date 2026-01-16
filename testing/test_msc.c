/*
 * test_msc.c - Unit tests for Multispecies Coalescent with Migration (MSC-M)
 *
 * Tests cover:
 * 1. Zero migration rate (regression tests - should match original MSC behavior)
 * 2. Migration rate calculation
 * 3. Migration event handling
 * 4. Migration stops at divergence time
 * 5. Ancestral migration between internal nodes
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "unity.h"
#include "msc.h"
#include "species_tree.h"
#include "coalescent.h"
#include "bitarray.h"

char version[] = "test_msc";

/* Global test fixtures */
static gsl_rng* rng = NULL;

void setUp(void) {
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 12345);
    g_noSamples = 10;
    bitarray_pool_init(10);
}

void tearDown(void) {
    if (rng) {
        gsl_rng_free(rng);
        rng = NULL;
    }
    bitarray_pool_destroy();
}

/* ============== Newick Parser Tests ============== */

void test_parse_newick_with_internal_names(void) {
    const char* newick = "((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;";
    species_tree* tree = species_tree_parse_newick(newick);

    TEST_ASSERT_NOT_NULL(tree);
    TEST_ASSERT_EQUAL_INT(5, tree->n_nodes);
    TEST_ASSERT_EQUAL_INT(3, tree->n_tips);

    /* Check root has name "ROOT" */
    TEST_ASSERT_EQUAL_STRING("ROOT", tree->root->name);

    /* Check internal node has name "AB" */
    species_node* ab = species_tree_find_node_by_name(tree, "AB");
    TEST_ASSERT_NOT_NULL(ab);
    TEST_ASSERT_EQUAL_INT(0, ab->is_tip);

    /* Check tip nodes */
    species_node* a = species_tree_find_node_by_name(tree, "A");
    species_node* b = species_tree_find_node_by_name(tree, "B");
    species_node* c = species_tree_find_node_by_name(tree, "C");
    TEST_ASSERT_NOT_NULL(a);
    TEST_ASSERT_NOT_NULL(b);
    TEST_ASSERT_NOT_NULL(c);
    TEST_ASSERT_EQUAL_INT(1, a->is_tip);
    TEST_ASSERT_EQUAL_INT(1, b->is_tip);
    TEST_ASSERT_EQUAL_INT(1, c->is_tip);

    species_tree_free(tree);
}

void test_parse_newick_without_internal_names(void) {
    /* Original format without internal node names should still work */
    const char* newick = "((A:300#100,B:300#100):200#100,C:500#100)#100;";
    species_tree* tree = species_tree_parse_newick(newick);

    TEST_ASSERT_NOT_NULL(tree);
    TEST_ASSERT_EQUAL_INT(5, tree->n_nodes);
    TEST_ASSERT_EQUAL_INT(3, tree->n_tips);

    /* Internal nodes should have empty names */
    TEST_ASSERT_EQUAL_STRING("", tree->root->name);
    TEST_ASSERT_EQUAL_STRING("", tree->root->left->name);

    species_tree_free(tree);
}

void test_find_node_by_name(void) {
    const char* newick = "((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;";
    species_tree* tree = species_tree_parse_newick(newick);
    TEST_ASSERT_NOT_NULL(tree);

    /* Find tip nodes */
    TEST_ASSERT_NOT_NULL(species_tree_find_node_by_name(tree, "A"));
    TEST_ASSERT_NOT_NULL(species_tree_find_node_by_name(tree, "B"));
    TEST_ASSERT_NOT_NULL(species_tree_find_node_by_name(tree, "C"));

    /* Find internal nodes */
    TEST_ASSERT_NOT_NULL(species_tree_find_node_by_name(tree, "AB"));
    TEST_ASSERT_NOT_NULL(species_tree_find_node_by_name(tree, "ROOT"));

    /* Non-existent node */
    TEST_ASSERT_NULL(species_tree_find_node_by_name(tree, "XYZ"));

    species_tree_free(tree);
}

/* ============== Migration Parsing Tests ============== */

void test_parse_control_file_with_migration(void) {
    /* Create a temporary control file */
    const char* filename = "/tmp/test_msc_control.ctl";
    FILE* f = fopen(filename, "w");
    TEST_ASSERT_NOT_NULL(f);
    fprintf(f, "species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:2,B:2,C:2\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:0.5,B->A:0.3\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_NOT_NULL(params->migrations);
    TEST_ASSERT_EQUAL_INT(2, params->migrations->n_bands);

    /* Check migration band values */
    int found_a_to_b = 0, found_b_to_a = 0;
    for (int i = 0; i < params->migrations->n_bands; i++) {
        migration_band* band = &params->migrations->bands[i];
        species_node* from = NULL;
        species_node* to = NULL;
        for (int j = 0; j < params->tree->n_nodes; j++) {
            if (params->tree->nodes[j]->id == band->from_pop)
                from = params->tree->nodes[j];
            if (params->tree->nodes[j]->id == band->to_pop)
                to = params->tree->nodes[j];
        }
        if (from && to && strcmp(from->name, "A") == 0 && strcmp(to->name, "B") == 0) {
            TEST_ASSERT_EQUAL_FLOAT(0.5, band->M);
            found_a_to_b = 1;
        }
        if (from && to && strcmp(from->name, "B") == 0 && strcmp(to->name, "A") == 0) {
            TEST_ASSERT_EQUAL_FLOAT(0.3, band->M);
            found_b_to_a = 1;
        }
    }
    TEST_ASSERT_EQUAL_INT(1, found_a_to_b);
    TEST_ASSERT_EQUAL_INT(1, found_b_to_a);

    msc_free_params(params);
    remove(filename);
}

void test_parse_control_file_without_migration(void) {
    /* MSC without migration should still work */
    const char* filename = "/tmp/test_msc_no_mig.ctl";
    FILE* f = fopen(filename, "w");
    TEST_ASSERT_NOT_NULL(f);
    fprintf(f, "species_tree: ((A:300#100,B:300#100):200#100,C:500#100)#100;\n");
    fprintf(f, "samples: A:2,B:2,C:2\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_NULL(params->migrations);

    msc_free_params(params);
    remove(filename);
}

void test_parse_ancestral_migration(void) {
    /* Test migration between ancestral populations */
    const char* filename = "/tmp/test_msc_ancestral.ctl";
    FILE* f = fopen(filename, "w");
    TEST_ASSERT_NOT_NULL(f);
    fprintf(f, "species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:2,B:2,C:2\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: AB->C:0.2,C->AB:0.1\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_NOT_NULL(params->migrations);
    TEST_ASSERT_EQUAL_INT(2, params->migrations->n_bands);

    msc_free_params(params);
    remove(filename);
}

/* ============== Population Existence Tests ============== */

void test_population_exists_at_time(void) {
    const char* newick = "((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;";
    species_tree* tree = species_tree_parse_newick(newick);
    TEST_ASSERT_NOT_NULL(tree);
    species_tree_compute_divergence_times(tree);

    species_node* a = species_tree_find_node_by_name(tree, "A");
    species_node* c = species_tree_find_node_by_name(tree, "C");
    species_node* ab = species_tree_find_node_by_name(tree, "AB");
    species_node* root = species_tree_find_node_by_name(tree, "ROOT");

    /* Tips A, B, C exist from t=0 until their parent's divergence time */
    /* A, B parent is AB which has divergence_time=300 */
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(a, 0));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(a, 100));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(a, 299));
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(a, 300));
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(a, 400));

    /* C's parent is ROOT which has divergence_time=500 */
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(c, 0));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(c, 300));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(c, 499));
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(c, 500));

    /* AB exists from t=300 until ROOT's divergence_time=500 */
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(ab, 0));
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(ab, 299));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(ab, 300));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(ab, 400));
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(ab, 500));

    /* ROOT exists from t=500 onwards (no parent) */
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(root, 0));
    TEST_ASSERT_EQUAL_INT(0, msc_population_exists_at_time(root, 499));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(root, 500));
    TEST_ASSERT_EQUAL_INT(1, msc_population_exists_at_time(root, 1000));

    species_tree_free(tree);
}

/* ============== Rate Calculation Tests ============== */

void test_rate_calculation_no_migration(void) {
    const char* filename = "/tmp/test_rate_no_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    species_tree_compute_divergence_times(params->tree);

    g_noSamples = params->total_samples;
    bitarray_pool_destroy();
    bitarray_pool_init(params->total_samples);
    bitarray* mrca = bitarray_full(params->total_samples);

    chrsample* sample = msc_create_sample(params->tree, params->total_samples);
    TEST_ASSERT_NOT_NULL(sample);

    msc_rates* rates = msc_calculate_rates(sample, params->tree, params, 0.0, mrca);
    TEST_ASSERT_NOT_NULL(rates);

    /* With no migration, total_mig_rate should be 0 */
    TEST_ASSERT_EQUAL_FLOAT(0.0, rates->total_mig_rate);
    TEST_ASSERT_EQUAL_INT(0, rates->n_mig_pairs);

    /* Coalescence rate: n(n-1)/(2*theta) for each pop with n=3, theta=100 */
    /* = 3*2/(2*100) = 0.03 per population, 0.06 total */
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.06, rates->total_coal_rate);

    msc_free_rates(rates);
    delete_sample(sample);
    bitarray_free(mrca);
    msc_free_params(params);
    remove(filename);
}

void test_rate_calculation_with_migration(void) {
    const char* filename = "/tmp/test_rate_with_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:1.0,B->A:2.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    species_tree_compute_divergence_times(params->tree);

    g_noSamples = params->total_samples;
    bitarray_pool_destroy();
    bitarray_pool_init(params->total_samples);
    bitarray* mrca = bitarray_full(params->total_samples);

    chrsample* sample = msc_create_sample(params->tree, params->total_samples);
    TEST_ASSERT_NOT_NULL(sample);

    msc_rates* rates = msc_calculate_rates(sample, params->tree, params, 0.0, mrca);
    TEST_ASSERT_NOT_NULL(rates);

    /* Migration rate = n_to * M / 2 for each band */
    /* A->B: n_B=3, M=1.0 => rate = 3 * 1.0 / 2 = 1.5 */
    /* B->A: n_A=3, M=2.0 => rate = 3 * 2.0 / 2 = 3.0 */
    /* Total: 4.5 */
    TEST_ASSERT_EQUAL_INT(2, rates->n_mig_pairs);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 4.5, rates->total_mig_rate);

    msc_free_rates(rates);
    delete_sample(sample);
    bitarray_free(mrca);
    msc_free_params(params);
    remove(filename);
}

void test_migration_rate_zero_after_divergence(void) {
    const char* filename = "/tmp/test_mig_after_div.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:1.0,B->A:1.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    species_tree_compute_divergence_times(params->tree);

    g_noSamples = params->total_samples;
    bitarray_pool_destroy();
    bitarray_pool_init(params->total_samples);
    bitarray* mrca = bitarray_full(params->total_samples);

    chrsample* sample = msc_create_sample(params->tree, params->total_samples);
    TEST_ASSERT_NOT_NULL(sample);

    /* Before divergence (t=0), migration should be active */
    msc_rates* rates = msc_calculate_rates(sample, params->tree, params, 0.0, mrca);
    TEST_ASSERT_NOT_NULL(rates);
    TEST_ASSERT_GREATER_THAN(0.0, rates->total_mig_rate);
    msc_free_rates(rates);

    /* After divergence (t=500), populations A and B no longer exist */
    /* All lineages should be in ROOT, so migration between A and B is invalid */
    rates = msc_calculate_rates(sample, params->tree, params, 500.0, mrca);
    TEST_ASSERT_NOT_NULL(rates);
    TEST_ASSERT_EQUAL_FLOAT(0.0, rates->total_mig_rate);
    msc_free_rates(rates);

    delete_sample(sample);
    bitarray_free(mrca);
    msc_free_params(params);
    remove(filename);
}

/* ============== Migration Selection Tests ============== */

void test_select_migration_pair(void) {
    /* Create mock rates structure */
    msc_rates rates;
    rates.mig_rates = malloc(3 * sizeof(mig_rate_info));
    rates.n_mig_pairs = 3;
    rates.mig_rates[0].from_pop = 0;
    rates.mig_rates[0].to_pop = 1;
    rates.mig_rates[0].rate = 1.0;
    rates.mig_rates[1].from_pop = 1;
    rates.mig_rates[1].to_pop = 0;
    rates.mig_rates[1].rate = 2.0;
    rates.mig_rates[2].from_pop = 0;
    rates.mig_rates[2].to_pop = 2;
    rates.mig_rates[2].rate = 3.0;
    rates.total_mig_rate = 6.0;

    /* Run many selections and verify distribution */
    int counts[3] = {0, 0, 0};
    int n_trials = 10000;
    for (int i = 0; i < n_trials; i++) {
        int idx = msc_select_migration_pair(rng, &rates);
        TEST_ASSERT_GREATER_OR_EQUAL(0, idx);
        TEST_ASSERT_LESS_THAN(3, idx);
        counts[idx]++;
    }

    /* Check approximate distribution (1:2:3 ratio) */
    /* With 10000 trials, expected: ~1667, ~3333, ~5000 */
    /* Allow 10% tolerance */
    TEST_ASSERT_INT_WITHIN(500, 1667, counts[0]);
    TEST_ASSERT_INT_WITHIN(500, 3333, counts[1]);
    TEST_ASSERT_INT_WITHIN(500, 5000, counts[2]);

    free(rates.mig_rates);
}

void test_select_lineage_in_pop(void) {
    const char* filename = "/tmp/test_select_lineage.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:5,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);

    g_noSamples = params->total_samples;
    bitarray_pool_destroy();
    bitarray_pool_init(params->total_samples);

    chrsample* sample = msc_create_sample(params->tree, params->total_samples);
    TEST_ASSERT_NOT_NULL(sample);

    species_node* a = species_tree_find_species(params->tree, "A");
    species_node* b = species_tree_find_species(params->tree, "B");
    TEST_ASSERT_NOT_NULL(a);
    TEST_ASSERT_NOT_NULL(b);

    /* Select lineages from population A (should have 5) */
    for (int i = 0; i < 100; i++) {
        int idx = msc_select_lineage_in_pop(rng, sample, a->id);
        TEST_ASSERT_GREATER_OR_EQUAL(0, idx);
        TEST_ASSERT_LESS_THAN(sample->count, idx);
        TEST_ASSERT_EQUAL_INT(a->id, sample->chrs[idx]->population_id);
    }

    /* Select lineages from population B (should have 3) */
    for (int i = 0; i < 100; i++) {
        int idx = msc_select_lineage_in_pop(rng, sample, b->id);
        TEST_ASSERT_GREATER_OR_EQUAL(0, idx);
        TEST_ASSERT_LESS_THAN(sample->count, idx);
        TEST_ASSERT_EQUAL_INT(b->id, sample->chrs[idx]->population_id);
    }

    /* Select from non-existent population should return -1 */
    int idx = msc_select_lineage_in_pop(rng, sample, 999);
    TEST_ASSERT_EQUAL_INT(-1, idx);

    delete_sample(sample);
    msc_free_params(params);
    remove(filename);
}

/* ============== Full Simulation Tests ============== */

void test_simulation_no_migration_completes(void) {
    /* Verify MSC without migration completes successfully */
    const char* filename = "/tmp/test_sim_no_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_EQUAL_INT(0, msc_validate_params(params));

    int result = msc_simulate(params);
    TEST_ASSERT_EQUAL_INT(0, result);

    msc_free_params(params);
    remove(filename);
}

void test_simulation_with_migration_completes(void) {
    /* Verify MSC with migration completes successfully */
    const char* filename = "/tmp/test_sim_with_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:1.0,B->A:1.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_EQUAL_INT(0, msc_validate_params(params));

    int result = msc_simulate(params);
    TEST_ASSERT_EQUAL_INT(0, result);

    msc_free_params(params);
    remove(filename);
}

void test_simulation_ancestral_migration_completes(void) {
    /* Verify MSC with ancestral migration completes successfully */
    const char* filename = "/tmp/test_sim_anc_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:2,B:2,C:2\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:0.5,B->A:0.5,AB->C:0.3,C->AB:0.3\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_EQUAL_INT(0, msc_validate_params(params));

    int result = msc_simulate(params);
    TEST_ASSERT_EQUAL_INT(0, result);

    msc_free_params(params);
    remove(filename);
}

void test_simulation_high_migration_rate(void) {
    /* With very high migration rate, should behave like panmictic population */
    const char* filename = "/tmp/test_sim_high_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:1000#100,B:1000#100)ROOT#100;\n");
    fprintf(f, "samples: A:4,B:4\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:100.0,B->A:100.0\n");  /* Very high M */
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_EQUAL_INT(0, msc_validate_params(params));

    int result = msc_simulate(params);
    TEST_ASSERT_EQUAL_INT(0, result);

    msc_free_params(params);
    remove(filename);
}

void test_simulation_asymmetric_migration(void) {
    /* Test asymmetric migration (different rates A->B vs B->A) */
    const char* filename = "/tmp/test_sim_asym_mig.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:2.0,B->A:0.5\n");  /* Asymmetric */
    fprintf(f, "seed: 54321\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_EQUAL_INT(0, msc_validate_params(params));

    int result = msc_simulate(params);
    TEST_ASSERT_EQUAL_INT(0, result);

    msc_free_params(params);
    remove(filename);
}

void test_simulation_zero_migration_matches_original(void) {
    /*
     * Test that MSC with zero migration rates behaves identically
     * to MSC without migration specification (regression test).
     * Both should produce the same results with the same seed.
     */
    const char* filename1 = "/tmp/test_no_mig_spec.ctl";
    const char* filename2 = "/tmp/test_zero_mig_spec.ctl";

    /* First: no migration specification */
    FILE* f = fopen(filename1, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "seed: 99999\n");
    fclose(f);

    /* Second: migration with zero rates */
    f = fopen(filename2, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:3,B:3\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:0.0,B->A:0.0\n");
    fprintf(f, "seed: 99999\n");
    fclose(f);

    /* Both should complete successfully */
    msc_params* params1 = msc_parse_control_file(filename1);
    TEST_ASSERT_NOT_NULL(params1);
    TEST_ASSERT_EQUAL_INT(0, msc_simulate(params1));

    msc_params* params2 = msc_parse_control_file(filename2);
    TEST_ASSERT_NOT_NULL(params2);
    TEST_ASSERT_EQUAL_INT(0, msc_simulate(params2));

    msc_free_params(params1);
    msc_free_params(params2);
    remove(filename1);
    remove(filename2);
}

/* ============== Memory Management Tests ============== */

void test_free_params_with_migration(void) {
    const char* filename = "/tmp/test_free_params.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:2,B:2\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:1.0,B->A:1.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    TEST_ASSERT_NOT_NULL(params->migrations);

    /* Should not crash or leak memory */
    msc_free_params(params);

    remove(filename);
}

void test_free_rates_with_migration(void) {
    const char* filename = "/tmp/test_free_rates.ctl";
    FILE* f = fopen(filename, "w");
    fprintf(f, "species_tree: (A:500#100,B:500#100)ROOT#100;\n");
    fprintf(f, "samples: A:2,B:2\n");
    fprintf(f, "recombination_rate: 0.0\n");
    fprintf(f, "mutation_rate: 0.0\n");
    fprintf(f, "migration: A->B:1.0,B->A:1.0\n");
    fprintf(f, "seed: 12345\n");
    fclose(f);

    msc_params* params = msc_parse_control_file(filename);
    TEST_ASSERT_NOT_NULL(params);
    species_tree_compute_divergence_times(params->tree);

    g_noSamples = params->total_samples;
    bitarray_pool_destroy();
    bitarray_pool_init(params->total_samples);
    bitarray* mrca = bitarray_full(params->total_samples);

    chrsample* sample = msc_create_sample(params->tree, params->total_samples);
    TEST_ASSERT_NOT_NULL(sample);

    msc_rates* rates = msc_calculate_rates(sample, params->tree, params, 0.0, mrca);
    TEST_ASSERT_NOT_NULL(rates);
    TEST_ASSERT_NOT_NULL(rates->mig_rates);

    /* Should not crash or leak memory */
    msc_free_rates(rates);

    delete_sample(sample);
    bitarray_free(mrca);
    msc_free_params(params);
    remove(filename);
}

/* ============== Main ============== */

int main(void) {
    UNITY_BEGIN();

    /* Newick parser tests */
    RUN_TEST(test_parse_newick_with_internal_names);
    RUN_TEST(test_parse_newick_without_internal_names);
    RUN_TEST(test_find_node_by_name);

    /* Migration parsing tests */
    RUN_TEST(test_parse_control_file_with_migration);
    RUN_TEST(test_parse_control_file_without_migration);
    RUN_TEST(test_parse_ancestral_migration);

    /* Population existence tests */
    RUN_TEST(test_population_exists_at_time);

    /* Rate calculation tests */
    RUN_TEST(test_rate_calculation_no_migration);
    RUN_TEST(test_rate_calculation_with_migration);
    RUN_TEST(test_migration_rate_zero_after_divergence);

    /* Migration selection tests */
    RUN_TEST(test_select_migration_pair);
    RUN_TEST(test_select_lineage_in_pop);

    /* Full simulation tests */
    RUN_TEST(test_simulation_no_migration_completes);
    RUN_TEST(test_simulation_with_migration_completes);
    RUN_TEST(test_simulation_ancestral_migration_completes);
    RUN_TEST(test_simulation_high_migration_rate);
    RUN_TEST(test_simulation_asymmetric_migration);
    RUN_TEST(test_simulation_zero_migration_matches_original);

    /* Memory management tests */
    RUN_TEST(test_free_params_with_migration);
    RUN_TEST(test_free_rates_with_migration);

    return UNITY_END();
}
