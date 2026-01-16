/*
 * msc.c - Multispecies Coalescent Simulation
 *
 * Implements the multispecies coalescent model where lineages can only
 * coalesce within the same population, and populations merge at divergence
 * times according to the species tree.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "msc.h"
#include "species_tree.h"
#include "coalescent.h"
#include "bitarray.h"

/* Global sample size for bitarray allocation (from coalescent.c) */
extern int g_noSamples;

/* Forward declarations for migration helpers */
static species_node* find_node_by_id(species_tree* tree, int id);
static int populations_contemporary(species_tree* tree, int pop1_id, int pop2_id, double current_time);

/* ============== Control File Parsing ============== */

/* Trim leading/trailing whitespace in place */
static char* trim(char* str) {
    while (isspace(*str)) str++;
    if (*str == '\0') return str;

    char* end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) end--;
    *(end + 1) = '\0';

    return str;
}

/* Parse "key: value" line, returns 1 if matched */
static int parse_key_value(char* line, const char* key, char** value) {
    size_t keylen = strlen(key);
    if (strncmp(line, key, keylen) != 0) return 0;

    char* p = line + keylen;
    while (*p && (isspace(*p) || *p == ':')) p++;

    *value = p;
    return 1;
}

/* Parse "species:count,species:count,..." sample specification */
static int parse_samples(const char* str, species_tree* tree) {
    char* copy = strdup(str);
    if (!copy) return -1;

    char* token = strtok(copy, ",");
    while (token) {
        token = trim(token);

        char* colon = strchr(token, ':');
        if (!colon) {
            fprintf(stderr, "Error: Invalid sample format '%s', expected 'species:count'\n", token);
            free(copy);
            return -1;
        }

        *colon = '\0';
        char* species_name = trim(token);
        int count = atoi(colon + 1);

        if (count <= 0) {
            fprintf(stderr, "Error: Invalid sample count for species '%s'\n", species_name);
            free(copy);
            return -1;
        }

        species_node* sp = species_tree_find_species(tree, species_name);
        if (!sp) {
            fprintf(stderr, "Error: Species '%s' not found in tree\n", species_name);
            free(copy);
            return -1;
        }

        sp->sample_size = count;
        token = strtok(NULL, ",");
    }

    free(copy);
    return 0;
}

/* Parse boolean value */
static int parse_bool(const char* str) {
    if (strcasecmp(str, "true") == 0 || strcmp(str, "1") == 0 ||
        strcasecmp(str, "yes") == 0) {
        return 1;
    }
    return 0;
}

/*
 * Parse migration specification: "A->B:0.5,B->A:0.3,AB->C:0.1,..."
 * Returns 0 on success, -1 on error.
 */
static int parse_migrations(const char* str, msc_params* params) {
    if (!str || !*str) return 0;  /* No migration is okay */

    char* copy = strdup(str);
    if (!copy) return -1;

    /* Count migration bands */
    int n_bands = 0;
    for (const char* p = str; *p; p++) {
        if (*p == ',') n_bands++;
    }
    n_bands++;  /* One more than number of commas */

    /* Allocate migration list */
    params->migrations = malloc(sizeof(migration_band_list));
    if (!params->migrations) {
        free(copy);
        return -1;
    }
    params->migrations->bands = malloc(n_bands * sizeof(migration_band));
    if (!params->migrations->bands) {
        free(params->migrations);
        params->migrations = NULL;
        free(copy);
        return -1;
    }
    params->migrations->n_bands = 0;

    /* Parse each band: "source->dest:M" */
    char* token = strtok(copy, ",");
    while (token) {
        token = trim(token);

        /* Find "->" separator */
        char* arrow = strstr(token, "->");
        if (!arrow) {
            fprintf(stderr, "Error: Invalid migration format '%s', expected 'source->dest:M'\n", token);
            free(copy);
            return -1;
        }

        /* Find ":" separator */
        char* colon = strchr(arrow, ':');
        if (!colon) {
            fprintf(stderr, "Error: Invalid migration format '%s', missing rate after ':'\n", token);
            free(copy);
            return -1;
        }

        /* Extract source, dest, and rate */
        *arrow = '\0';
        *colon = '\0';
        char* source_name = trim(token);
        char* dest_name = trim(arrow + 2);
        double M = atof(colon + 1);

        if (M < 0) {
            fprintf(stderr, "Error: Negative migration rate for %s->%s\n", source_name, dest_name);
            free(copy);
            return -1;
        }

        /* Look up source and destination nodes (can be tips or internal) */
        species_node* source = species_tree_find_node_by_name(params->tree, source_name);
        species_node* dest = species_tree_find_node_by_name(params->tree, dest_name);

        if (!source) {
            fprintf(stderr, "Error: Migration source '%s' not found in tree\n", source_name);
            free(copy);
            return -1;
        }
        if (!dest) {
            fprintf(stderr, "Error: Migration destination '%s' not found in tree\n", dest_name);
            free(copy);
            return -1;
        }

        if (source->id == dest->id) {
            fprintf(stderr, "Error: Cannot have migration from population to itself (%s)\n", source_name);
            free(copy);
            return -1;
        }

        /* Add the migration band */
        migration_band* band = &params->migrations->bands[params->migrations->n_bands];
        band->from_pop = source->id;
        band->to_pop = dest->id;
        band->M = M;
        params->migrations->n_bands++;

        token = strtok(NULL, ",");
    }

    free(copy);
    return 0;
}

msc_params* msc_parse_control_file(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open control file '%s'\n", filename);
        return NULL;
    }

    msc_params* params = calloc(1, sizeof(msc_params));
    if (!params) {
        fclose(f);
        return NULL;
    }

    /* Set defaults */
    params->recombination_rate = 0.0;
    params->mutation_rate = 0.0;
    params->seed = 0;
    params->output_gene_trees = 0;
    params->output_vcf = 0;
    params->vcf_filename = NULL;
    params->verbose = 0;
    params->migrations = NULL;

    char line[4096];
    char* value;
    int line_num = 0;

    while (fgets(line, sizeof(line), f)) {
        line_num++;

        /* Remove newline */
        char* nl = strchr(line, '\n');
        if (nl) *nl = '\0';

        /* Skip comments and empty lines */
        char* trimmed = trim(line);
        if (*trimmed == '\0' || *trimmed == '#') continue;

        if (parse_key_value(trimmed, "species_tree", &value)) {
            value = trim(value);
            params->tree = species_tree_parse_newick(value);
            if (!params->tree) {
                fprintf(stderr, "Error parsing species tree at line %d\n", line_num);
                msc_free_params(params);
                fclose(f);
                return NULL;
            }
        }
        else if (parse_key_value(trimmed, "samples", &value)) {
            if (!params->tree) {
                fprintf(stderr, "Error: 'samples' must come after 'species_tree' (line %d)\n", line_num);
                msc_free_params(params);
                fclose(f);
                return NULL;
            }
            if (parse_samples(trim(value), params->tree) < 0) {
                msc_free_params(params);
                fclose(f);
                return NULL;
            }
        }
        else if (parse_key_value(trimmed, "recombination_rate", &value)) {
            params->recombination_rate = atof(trim(value));
        }
        else if (parse_key_value(trimmed, "mutation_rate", &value)) {
            params->mutation_rate = atof(trim(value));
        }
        else if (parse_key_value(trimmed, "seed", &value)) {
            params->seed = (unsigned long)atol(trim(value));
        }
        else if (parse_key_value(trimmed, "output_gene_trees", &value)) {
            params->output_gene_trees = parse_bool(trim(value));
        }
        else if (parse_key_value(trimmed, "output_vcf", &value)) {
            params->output_vcf = parse_bool(trim(value));
        }
        else if (parse_key_value(trimmed, "vcf_file", &value)) {
            params->vcf_filename = strdup(trim(value));
        }
        else if (parse_key_value(trimmed, "verbose", &value)) {
            params->verbose = parse_bool(trim(value));
        }
        else if (parse_key_value(trimmed, "migration", &value)) {
            if (!params->tree) {
                fprintf(stderr, "Error: 'migration' must come after 'species_tree' (line %d)\n", line_num);
                msc_free_params(params);
                fclose(f);
                return NULL;
            }
            if (parse_migrations(trim(value), params) < 0) {
                msc_free_params(params);
                fclose(f);
                return NULL;
            }
        }
        else {
            fprintf(stderr, "Warning: Unknown parameter at line %d: %s\n", line_num, trimmed);
        }
    }

    fclose(f);

    /* Compute derived parameters */
    if (params->tree) {
        species_tree_compute_divergence_times(params->tree);
        params->total_samples = species_tree_total_samples(params->tree);
    }

    return params;
}

int msc_validate_params(msc_params* params) {
    if (!params) {
        fprintf(stderr, "Error: NULL parameters\n");
        return -1;
    }

    if (!params->tree) {
        fprintf(stderr, "Error: species_tree is required\n");
        return -1;
    }

    if (params->total_samples < 2) {
        fprintf(stderr, "Error: Total sample size must be at least 2\n");
        return -1;
    }

    /* Check all tips have sample sizes */
    for (int i = 0; i < params->tree->n_tips; i++) {
        if (params->tree->tips[i]->sample_size <= 0) {
            fprintf(stderr, "Error: No sample size specified for species '%s'\n",
                    params->tree->tips[i]->name);
            return -1;
        }
    }

    /* Check all nodes have theta values */
    for (int i = 0; i < params->tree->n_nodes; i++) {
        if (params->tree->nodes[i]->theta <= 0) {
            fprintf(stderr, "Error: Invalid or missing theta for node %d\n",
                    params->tree->nodes[i]->id);
            return -1;
        }
    }

    if (params->output_vcf && !params->vcf_filename) {
        fprintf(stderr, "Error: vcf_file is required when output_vcf is true\n");
        return -1;
    }

    return 0;
}

void msc_free_params(msc_params* params) {
    if (!params) return;
    if (params->tree) species_tree_free(params->tree);
    if (params->vcf_filename) free(params->vcf_filename);
    if (params->migrations) msc_free_migrations(params->migrations);
    free(params);
}

/* ============== Sample Creation ============== */

chrsample* msc_create_sample(species_tree* tree, int total_samples) {
    chrsample* sample = malloc(sizeof(chrsample));
    if (!sample) return NULL;

    sample->capacity = total_samples * 4;  /* Room for recombination */
    sample->chrs = malloc(sample->capacity * sizeof(chromosome*));
    if (!sample->chrs) {
        free(sample);
        return NULL;
    }

    sample->count = 0;
    sample->ancLength = (double)total_samples;
    sample->activeAncLength = (double)total_samples;
    sample->activeAncValid = 0;

    /* Create chromosomes for each species */
    int global_idx = 0;
    for (int i = 0; i < tree->n_tips; i++) {
        species_node* sp = tree->tips[i];
        for (int j = 0; j < sp->sample_size; j++) {
            chromosome* chr = malloc(sizeof(chromosome));
            if (!chr) {
                /* Cleanup on failure */
                for (int k = 0; k < sample->count; k++) {
                    delete_anc(sample->chrs[k]->anc);
                    free(sample->chrs[k]);
                }
                free(sample->chrs);
                free(sample);
                return NULL;
            }

            chr->anc = calloc(1, sizeof(ancestry));
            chr->anc->abits = bitarray_singleton(total_samples, global_idx);
            chr->anc->position = 1.0;
            chr->anc->next = NULL;
            chr->anc->is_mrca = 0;
            chr->ancLen = 1.0;
            chr->activeLen = 1.0;
            chr->activeLenValid = 0;
            chr->population_id = sp->id;  /* Tag with population */

            sample->chrs[sample->count++] = chr;
            global_idx++;
        }
        sp->active_lineages = sp->sample_size;
    }

    return sample;
}

/* ============== Rate Calculation ============== */

msc_rates* msc_calculate_rates(chrsample* sample, species_tree* tree,
                                msc_params* params, double current_time,
                                const bitarray* mrca) {
    msc_rates* rates = calloc(1, sizeof(msc_rates));
    if (!rates) return NULL;

    /* First count lineages per population */
    msc_count_lineages(sample, tree);

    /* Allocate space for per-population rates */
    rates->pop_rates = calloc(tree->n_nodes, sizeof(pop_rate_info));
    if (!rates->pop_rates) {
        free(rates);
        return NULL;
    }

    rates->n_active_pops = 0;
    rates->total_coal_rate = 0;

    /* Calculate coalescence rate for each population with >= 2 lineages */
    for (int i = 0; i < tree->n_nodes; i++) {
        species_node* node = tree->nodes[i];
        int n = node->active_lineages;

        if (n >= 2) {
            /* Coalescence rate = n(n-1) / (2 * theta) */
            double coal_rate = (double)n * (n - 1) / (2.0 * node->theta);

            rates->pop_rates[rates->n_active_pops].pop_id = node->id;
            rates->pop_rates[rates->n_active_pops].n_lineages = n;
            rates->pop_rates[rates->n_active_pops].coal_rate = coal_rate;
            rates->n_active_pops++;

            rates->total_coal_rate += coal_rate;
        }
    }

    /* Calculate migration rates if migration is configured */
    rates->mig_rates = NULL;
    rates->n_mig_pairs = 0;
    rates->total_mig_rate = 0;

    if (params->migrations && params->migrations->n_bands > 0) {
        rates->mig_rates = calloc(params->migrations->n_bands, sizeof(mig_rate_info));
        if (!rates->mig_rates) {
            free(rates->pop_rates);
            free(rates);
            return NULL;
        }

        for (int b = 0; b < params->migrations->n_bands; b++) {
            migration_band* band = &params->migrations->bands[b];

            /* Find lineage count in "to" population (source going backwards) */
            species_node* to_node = find_node_by_id(tree, band->to_pop);
            int n_to = to_node ? to_node->active_lineages : 0;

            /* Only add if destination has lineages AND both pops are contemporary */
            if (n_to > 0 && populations_contemporary(tree, band->from_pop, band->to_pop, current_time)) {
                double mig_rate = n_to * band->M / 2.0;
                rates->mig_rates[rates->n_mig_pairs].from_pop = band->from_pop;
                rates->mig_rates[rates->n_mig_pairs].to_pop = band->to_pop;
                rates->mig_rates[rates->n_mig_pairs].rate = mig_rate;
                rates->total_mig_rate += mig_rate;
                rates->n_mig_pairs++;
            }
        }
    }

    /* Recombination and mutation rates scale with active ancestry length */
    double anc_len = getActiveAncLength(sample, mrca);
    rates->rec_rate = params->recombination_rate * anc_len;
    rates->mut_rate = params->mutation_rate * anc_len;
    rates->total_rate = rates->total_coal_rate + rates->total_mig_rate +
                        rates->rec_rate + rates->mut_rate;

    return rates;
}

void msc_free_rates(msc_rates* rates) {
    if (!rates) return;
    free(rates->pop_rates);
    free(rates->mig_rates);
    free(rates);
}

/* ============== Population Selection ============== */

int msc_select_coal_population(gsl_rng* r, msc_rates* rates) {
    if (rates->n_active_pops == 0) return -1;

    double u = gsl_rng_uniform_pos(r) * rates->total_coal_rate;
    double cumsum = 0;

    for (int i = 0; i < rates->n_active_pops; i++) {
        cumsum += rates->pop_rates[i].coal_rate;
        if (u <= cumsum) {
            return rates->pop_rates[i].pop_id;
        }
    }

    /* Fallback (should not happen) */
    return rates->pop_rates[rates->n_active_pops - 1].pop_id;
}

int msc_get_coal_pair_in_pop(gsl_rng* r, chrsample* sample,
                              int pop_id, coalescent_pair* pair) {
    /* Build list of chromosome indices in this population */
    int* indices = malloc(sample->count * sizeof(int));
    if (!indices) return -1;

    int n = 0;
    for (int i = 0; i < sample->count; i++) {
        if (sample->chrs[i]->population_id == pop_id) {
            indices[n++] = i;
        }
    }

    if (n < 2) {
        free(indices);
        return -1;  /* Cannot coalesce */
    }

    /* Pick two distinct indices */
    int i1 = gsl_rng_uniform_int(r, n);
    int i2 = gsl_rng_uniform_int(r, n - 1);
    if (i2 >= i1) i2++;

    pair->chr1 = indices[i1];
    pair->chr2 = indices[i2];

    free(indices);
    return 0;
}

/* ============== Divergence Handling ============== */

void msc_process_divergence(chrsample* sample, divergence_event* event) {
    species_node* parent = event->parent;
    species_node* c1 = event->child1;
    species_node* c2 = event->child2;

    /* Re-assign all chromosomes from children to parent population */
    int count = 0;
    for (int i = 0; i < sample->count; i++) {
        chromosome* chr = sample->chrs[i];
        if (chr->population_id == c1->id || chr->population_id == c2->id) {
            chr->population_id = parent->id;
            count++;
        }
    }

    /* Update lineage counts */
    parent->active_lineages = count;
    c1->active_lineages = 0;
    c2->active_lineages = 0;
}

/* ============== Lineage Counting ============== */

void msc_count_lineages(chrsample* sample, species_tree* tree) {
    /* Reset all counts */
    for (int i = 0; i < tree->n_nodes; i++) {
        tree->nodes[i]->active_lineages = 0;
    }

    /* Count lineages per population */
    for (int i = 0; i < sample->count; i++) {
        int pop_id = sample->chrs[i]->population_id;
        /* Find the node with this ID */
        for (int j = 0; j < tree->n_nodes; j++) {
            if (tree->nodes[j]->id == pop_id) {
                tree->nodes[j]->active_lineages++;
                break;
            }
        }
    }
}

/* ============== Migration Functions ============== */

/* Find species node by ID */
static species_node* find_node_by_id(species_tree* tree, int id) {
    for (int i = 0; i < tree->n_nodes; i++) {
        if (tree->nodes[i]->id == id) {
            return tree->nodes[i];
        }
    }
    return NULL;
}

/*
 * Check if a population exists at the given time.
 * A population exists at time t if:
 * - For tips: always exists from t=0 until parent's divergence_time
 * - For internal nodes: exists from its divergence_time until parent's divergence_time
 *
 * In our simulation going backwards in time:
 * - Tips exist when current_time < tip's parent divergence_time
 * - Internal nodes exist when current_time >= node's divergence_time AND
 *   current_time < parent's divergence_time (or parent is NULL for root)
 */
int msc_population_exists_at_time(species_node* node, double current_time) {
    if (!node) return 0;

    if (node->is_tip) {
        /* Tips exist from present (t=0) until their parent's divergence time */
        if (!node->parent) return 1;  /* Should not happen for tips */
        return current_time < node->parent->divergence_time;
    } else {
        /* Internal nodes exist from their divergence time onwards */
        if (current_time < node->divergence_time) return 0;

        /* Until parent's divergence time (or forever if root) */
        if (!node->parent) return 1;
        return current_time < node->parent->divergence_time;
    }
}

/*
 * Check if two populations are contemporary (both exist at current time).
 * Migration can only occur between contemporary populations.
 */
static int populations_contemporary(species_tree* tree, int pop1_id, int pop2_id, double current_time) {
    species_node* node1 = find_node_by_id(tree, pop1_id);
    species_node* node2 = find_node_by_id(tree, pop2_id);

    if (!node1 || !node2) return 0;

    return msc_population_exists_at_time(node1, current_time) &&
           msc_population_exists_at_time(node2, current_time);
}

int msc_select_migration_pair(gsl_rng* r, msc_rates* rates) {
    if (rates->n_mig_pairs == 0) return -1;

    double u = gsl_rng_uniform_pos(r) * rates->total_mig_rate;
    double cumsum = 0;

    for (int i = 0; i < rates->n_mig_pairs; i++) {
        cumsum += rates->mig_rates[i].rate;
        if (u <= cumsum) {
            return i;
        }
    }

    /* Fallback (should not happen) */
    return rates->n_mig_pairs - 1;
}

int msc_select_lineage_in_pop(gsl_rng* r, chrsample* sample, int pop_id) {
    /* Count lineages in this population */
    int count = 0;
    for (int i = 0; i < sample->count; i++) {
        if (sample->chrs[i]->population_id == pop_id) {
            count++;
        }
    }

    if (count == 0) return -1;

    /* Select random one */
    int target = gsl_rng_uniform_int(r, count);
    int seen = 0;
    for (int i = 0; i < sample->count; i++) {
        if (sample->chrs[i]->population_id == pop_id) {
            if (seen == target) {
                return i;
            }
            seen++;
        }
    }

    return -1;  /* Should not happen */
}

void msc_free_migrations(migration_band_list* migrations) {
    if (!migrations) return;
    free(migrations->bands);
    free(migrations);
}

/* ============== Main Simulation ============== */

int msc_simulate(msc_params* params) {
    /* Initialize RNG */
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    if (params->seed == 0) {
        params->seed = (unsigned long)time(NULL);
    }
    gsl_rng_set(r, params->seed);

    /* Set global sample size for bitarray allocation */
    g_noSamples = params->total_samples;

    /* Initialize bitarray pool */
    bitarray_pool_init(params->total_samples);

    /* Create MRCA bitarray (all bits set) */
    bitarray* mrca = bitarray_full(params->total_samples);

    /* Create initial sample */
    chrsample* sample = msc_create_sample(params->tree, params->total_samples);
    if (!sample) {
        fprintf(stderr, "Error: Failed to create sample\n");
        bitarray_free(mrca);
        gsl_rng_free(r);
        return -1;
    }

    /* Create divergence event queue */
    divergence_event_list* div_events = species_tree_create_divergence_events(params->tree);
    divergence_event* next_div = div_events ? div_events->head : NULL;

    /* Data structures for tracking events */
    struct coalescent_events* coalescent_list = NULL;
    struct coalescent_events* coalescent_list_tail = NULL;
    mutation* mutation_list = NULL;
    mutation* mutation_list_tail = NULL;

    double current_time = 0;
    unsigned int n_chrom = sample->count;
    int n_recombinations = 0;
    int n_mutations = 0;
    int n_coalescences = 0;
    int n_migrations = 0;

    if (params->verbose) {
        fprintf(stderr, "Starting MSC simulation:\n");
        fprintf(stderr, "  Total samples: %d\n", params->total_samples);
        fprintf(stderr, "  Species: %d\n", params->tree->n_tips);
        fprintf(stderr, "  Recombination rate: %.4f\n", params->recombination_rate);
        fprintf(stderr, "  Mutation rate: %.4f\n", params->mutation_rate);
        if (params->migrations) {
            fprintf(stderr, "  Migration bands: %d\n", params->migrations->n_bands);
        }
        fprintf(stderr, "  Seed: %lu\n", params->seed);
        species_tree_print(params->tree, stderr);
    }

    /* Main simulation loop */
    while (n_chrom > 1) {
        /* Calculate rates */
        msc_rates* rates = msc_calculate_rates(sample, params->tree,
                                                params, current_time, mrca);

        if (!rates || rates->total_rate <= 0) {
            /* No coalescence/recombination/mutation possible right now.
             * But if there are pending divergence events, we should wait
             * for them to merge populations, enabling further coalescence. */
            if (next_div) {
                /* Jump to next divergence event */
                current_time = next_div->time;
                if (params->verbose) {
                    fprintf(stderr, "Time %.1f: Divergence - merging populations %d and %d into %d\n",
                            current_time, next_div->child1->id, next_div->child2->id,
                            next_div->parent->id);
                }
                msc_process_divergence(sample, next_div);
                next_div = next_div->next;
                if (rates) msc_free_rates(rates);
                continue;
            }
            /* No divergence events and no events possible - truly done */
            if (rates) msc_free_rates(rates);
            break;
        }

        /* Time to next event */
        double dt = gsl_ran_exponential(r, 1.0 / rates->total_rate);
        double event_time = current_time + dt;

        /* Check if divergence event happens first */
        if (next_div && event_time > next_div->time) {
            /* Process divergence */
            current_time = next_div->time;
            if (params->verbose) {
                fprintf(stderr, "Time %.1f: Divergence - merging populations %d and %d into %d\n",
                        current_time, next_div->child1->id, next_div->child2->id,
                        next_div->parent->id);
            }
            msc_process_divergence(sample, next_div);
            next_div = next_div->next;
            msc_free_rates(rates);
            continue;
        }

        current_time = event_time;

        /* Determine event type */
        double u = gsl_rng_uniform_pos(r) * rates->total_rate;

        if (u < rates->total_coal_rate) {
            /* Coalescence event */
            int pop_id = msc_select_coal_population(r, rates);
            coalescent_pair pair;

            if (msc_get_coal_pair_in_pop(r, sample, pop_id, &pair) == 0) {
                /* Use existing coalescence function */
                coalescence(pair, &n_chrom, sample, mrca);
                n_coalescences++;

                /* Invalidate active length cache */
                invalidateActiveAncLength(sample);

                /* Track for gene tree output */
                if (params->output_gene_trees) {
                    updateCoalescentEvents(&coalescent_list, &coalescent_list_tail,
                                           sample, current_time);
                }
            }
        }
        else if (u < rates->total_coal_rate + rates->total_mig_rate) {
            /* Migration event */
            int pair_idx = msc_select_migration_pair(r, rates);
            if (pair_idx >= 0) {
                int to_pop = rates->mig_rates[pair_idx].to_pop;
                int from_pop = rates->mig_rates[pair_idx].from_pop;

                /* Select random lineage in to_pop and move it to from_pop */
                int chr_idx = msc_select_lineage_in_pop(r, sample, to_pop);
                if (chr_idx >= 0) {
                    sample->chrs[chr_idx]->population_id = from_pop;
                    n_migrations++;

                    if (params->verbose) {
                        /* Find node names for nice output */
                        species_node* from_node = find_node_by_id(params->tree, from_pop);
                        species_node* to_node = find_node_by_id(params->tree, to_pop);
                        const char* from_name = (from_node && from_node->name[0]) ? from_node->name : "?";
                        const char* to_name = (to_node && to_node->name[0]) ? to_node->name : "?";
                        fprintf(stderr, "Time %.1f: Migration - lineage from %s to %s\n",
                                current_time, to_name, from_name);
                    }

                    /* Update lineage counts (no need to invalidate active length) */
                    msc_count_lineages(sample, params->tree);
                }
            }
        }
        else if (u < rates->total_coal_rate + rates->total_mig_rate + rates->rec_rate) {
            /* Recombination event */
            double rec_pos = (u - rates->total_coal_rate - rates->total_mig_rate) / params->recombination_rate;
            recombination_event recEv;
            getRecEventActive(sample, rec_pos, &recEv, mrca);

            /* recombination() in coalescent.c now inherits population_id automatically */
            recombination(&n_chrom, recEv, sample);
            n_recombinations++;
        }
        else {
            /* Mutation event */
            double mut_pos = (u - rates->total_coal_rate - rates->total_mig_rate - rates->rec_rate) / params->mutation_rate;
            mutation* mutEv = malloc(sizeof(mutation));
            getMutEventActive(sample, mut_pos, mutEv, current_time, mrca);

            /* Add to mutation list */
            mutEv->next = NULL;
            if (!mutation_list) {
                mutation_list = mutEv;
                mutation_list_tail = mutEv;
            } else {
                mutation_list_tail->next = mutEv;
                mutation_list_tail = mutEv;
            }
            n_mutations++;
        }

        msc_free_rates(rates);

        /* In MSC, we don't check for MRCA early termination because
         * populations are isolated until divergence times. The loop
         * terminates when n_chrom == 1 (all lineages coalesced) or
         * when total_rate == 0 (only MRCA segments remain within
         * isolated populations that can't coalesce further). */
    }

    /* Print summary */
    printf("MSC Simulation Complete:\n");
    printf("  Coalescences: %d\n", n_coalescences);
    printf("  Migrations: %d\n", n_migrations);
    printf("  Recombinations: %d\n", n_recombinations);
    printf("  Mutations: %d\n", n_mutations);
    printf("  Final time: %.2f generations\n", current_time);

    /* Output VCF if requested */
    if (params->output_vcf && params->vcf_filename) {
        FILE* vcf_out = fopen(params->vcf_filename, "w");
        if (vcf_out) {
            /* Assume 1Mb chromosome for now - this should be configurable */
            long chrom_bases = 1000000;
            writeVCF(mutation_list, chrom_bases, params->total_samples, mrca, vcf_out, r);
            fclose(vcf_out);
            printf("  VCF written to: %s\n", params->vcf_filename);
        } else {
            fprintf(stderr, "Warning: Could not open VCF file '%s'\n", params->vcf_filename);
        }
    }

    /* Cleanup */
    delete_sample(sample);
    bitarray_free(mrca);
    gsl_rng_free(r);
    if (div_events) divergence_event_list_free(div_events);

    /* Free mutation list */
    mutation* mut = mutation_list;
    while (mut) {
        mutation* next = mut->next;
        if (mut->abits) bitarray_free(mut->abits);
        free(mut);
        mut = next;
    }

    /* Free coalescent event list */
    struct coalescent_events* ce = coalescent_list;
    while (ce) {
        struct coalescent_events* next = ce->next;
        /* Note: chr is owned by sample, already freed */
        free(ce);
        ce = next;
    }

    bitarray_pool_destroy();

    return 0;
}

/* ============== Utility Functions ============== */

void msc_print_summary(msc_params* params, FILE* out) {
    if (!params) return;

    fprintf(out, "MSC Parameters:\n");
    fprintf(out, "  Total samples: %d\n", params->total_samples);
    fprintf(out, "  Recombination rate: %.6f\n", params->recombination_rate);
    fprintf(out, "  Mutation rate: %.6f\n", params->mutation_rate);
    fprintf(out, "  Seed: %lu\n", params->seed);
    fprintf(out, "  Output gene trees: %s\n", params->output_gene_trees ? "yes" : "no");
    fprintf(out, "  Output VCF: %s\n", params->output_vcf ? "yes" : "no");

    if (params->tree) {
        fprintf(out, "\nSpecies Tree:\n");
        species_tree_print(params->tree, out);
    }
}
