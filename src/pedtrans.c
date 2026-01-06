/*
 * pedtrans.c - Pedigree Chromosome Transmission Simulator
 *
 * Implementation of efficient pedigree parsing, meiosis simulation,
 * and chromosome segment tracking.
 *
 * Copyright (C) 2025 Bruce Rannala
 * GNU Affero General Public License v3
 */

#include "pedtrans.h"
#include <time.h>
#include <math.h>

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

int compare_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* Ensure array has capacity for at least one more element */
static int ensure_capacity_individuals(pedigree* ped) {
    if (ped->n_individuals >= ped->capacity) {
        int new_capacity = ped->capacity * 2;
        ped_indiv* new_indivs = realloc(ped->individuals,
                                         new_capacity * sizeof(ped_indiv));
        if (!new_indivs) return -1;
        ped->individuals = new_indivs;

        char** new_names = realloc(ped->id_to_name,
                                    new_capacity * sizeof(char*));
        if (!new_names) return -1;
        ped->id_to_name = new_names;

        ped->capacity = new_capacity;
    }
    return 0;
}

static int ensure_capacity_segments(ped_chromosome* chr) {
    if (chr->n_segments >= chr->capacity) {
        int new_capacity = chr->capacity * 2;
        ped_segment* new_segs = realloc(chr->segments,
                                         new_capacity * sizeof(ped_segment));
        if (!new_segs) return -1;
        chr->segments = new_segs;
        chr->capacity = new_capacity;
    }
    return 0;
}

/* ============================================================================
 * PEDIGREE FUNCTIONS
 * ============================================================================ */

pedigree* create_pedigree(int initial_capacity) {
    pedigree* ped = malloc(sizeof(pedigree));
    if (!ped) return NULL;

    ped->individuals = malloc(initial_capacity * sizeof(ped_indiv));
    ped->id_to_name = malloc(initial_capacity * sizeof(char*));
    ped->name_to_id = g_hash_table_new_full(g_str_hash, g_str_equal, NULL, g_free);

    if (!ped->individuals || !ped->id_to_name || !ped->name_to_id) {
        free(ped->individuals);
        free(ped->id_to_name);
        if (ped->name_to_id) g_hash_table_destroy(ped->name_to_id);
        free(ped);
        return NULL;
    }

    ped->topo_order = NULL;
    ped->children = NULL;
    ped->child_offsets = NULL;
    ped->n_children = NULL;
    ped->n_individuals = 0;
    ped->n_founders = 0;
    ped->capacity = initial_capacity;

    return ped;
}

void free_pedigree(pedigree* ped) {
    if (!ped) return;

    /* Free individual names */
    for (int i = 0; i < ped->n_individuals; i++) {
        free(ped->id_to_name[i]);
    }
    free(ped->id_to_name);
    free(ped->individuals);
    free(ped->topo_order);
    free(ped->children);
    free(ped->child_offsets);
    free(ped->n_children);

    if (ped->name_to_id) {
        g_hash_table_destroy(ped->name_to_id);
    }

    free(ped);
}

int pedigree_get_id(pedigree* ped, const char* name) {
    gint* id_ptr = g_hash_table_lookup(ped->name_to_id, name);
    if (id_ptr) return *id_ptr;
    return -1;
}

const char* pedigree_get_name(pedigree* ped, int id) {
    if (id < 0 || id >= ped->n_individuals) return NULL;
    return ped->id_to_name[id];
}

/*
 * Get or create ID for a name
 * If name is "0", returns -1 (represents missing/founder parent)
 */
static int get_or_create_id(pedigree* ped, const char* name) {
    if (strcmp(name, "0") == 0) {
        return -1;
    }

    int existing_id = pedigree_get_id(ped, name);
    if (existing_id >= 0) {
        return existing_id;
    }

    /* Need to create new entry */
    if (ensure_capacity_individuals(ped) < 0) {
        return -2; /* Error */
    }

    int new_id = ped->n_individuals;

    /* Copy name */
    char* name_copy = strdup(name);
    if (!name_copy) return -2;

    /* Add to hash table */
    gint* id_ptr = g_new(gint, 1);
    *id_ptr = new_id;
    g_hash_table_insert(ped->name_to_id, name_copy, id_ptr);

    /* Store name reference */
    ped->id_to_name[new_id] = name_copy;

    /* Initialize individual (will be filled in later) */
    ped->individuals[new_id].id = new_id;
    ped->individuals[new_id].father_id = -1;
    ped->individuals[new_id].mother_id = -1;
    ped->individuals[new_id].is_founder = 0;
    ped->individuals[new_id].in_degree = 0;

    ped->n_individuals++;
    return new_id;
}

int pedigree_add_individual(pedigree* ped, const char* name,
                            const char* father, const char* mother) {
    /* Get or create ID for individual */
    int id = get_or_create_id(ped, name);
    if (id < 0) return -1;

    /* Get parent IDs */
    int father_id = get_or_create_id(ped, father);
    int mother_id = get_or_create_id(ped, mother);

    if (father_id == -2 || mother_id == -2) {
        return -1; /* Allocation error */
    }

    /* Update individual info */
    ped->individuals[id].father_id = father_id;
    ped->individuals[id].mother_id = mother_id;

    /* Determine if founder */
    if (father_id == -1 && mother_id == -1) {
        ped->individuals[id].is_founder = 1;
        ped->individuals[id].in_degree = 0;
        ped->n_founders++;
    } else {
        ped->individuals[id].is_founder = 0;
        /* in_degree = number of parents (can be 1 or 2 if one parent is "0") */
        ped->individuals[id].in_degree = 0;
        if (father_id >= 0) ped->individuals[id].in_degree++;
        if (mother_id >= 0) ped->individuals[id].in_degree++;
    }

    return id;
}

pedigree* parse_pedigree(FILE* fp) {
    pedigree* ped = create_pedigree(PED_INITIAL_CAPACITY);
    if (!ped) return NULL;

    char line[1024];
    char indiv[PED_MAX_NAME], father[PED_MAX_NAME], mother[PED_MAX_NAME];

    while (fgets(line, sizeof(line), fp)) {
        /* Skip empty lines and comments */
        char* p = line;
        while (*p && (*p == ' ' || *p == '\t')) p++;
        if (*p == '\0' || *p == '\n' || *p == '#') continue;

        /* Parse three fields */
        if (sscanf(line, "%63s %63s %63s", indiv, father, mother) != 3) {
            fprintf(stderr, "Error: malformed line: %s", line);
            free_pedigree(ped);
            return NULL;
        }

        if (pedigree_add_individual(ped, indiv, father, mother) < 0) {
            fprintf(stderr, "Error: failed to add individual %s\n", indiv);
            free_pedigree(ped);
            return NULL;
        }
    }

    /* Perform topological sort */
    if (topological_sort(ped) < 0) {
        fprintf(stderr, "Error: pedigree contains cycle or invalid structure\n");
        free_pedigree(ped);
        return NULL;
    }

    return ped;
}

pedigree* parse_pedigree_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        return NULL;
    }

    pedigree* ped = parse_pedigree(fp);
    fclose(fp);
    return ped;
}

void build_child_lists(pedigree* ped) {
    int n = ped->n_individuals;

    /* Free existing child lists if any */
    free(ped->children);
    free(ped->child_offsets);
    free(ped->n_children);

    /* Count children for each individual */
    ped->n_children = calloc(n, sizeof(int));
    if (!ped->n_children) return;

    for (int i = 0; i < n; i++) {
        int fid = ped->individuals[i].father_id;
        int mid = ped->individuals[i].mother_id;
        if (fid >= 0) ped->n_children[fid]++;
        if (mid >= 0) ped->n_children[mid]++;
    }

    /* Build offset array */
    ped->child_offsets = malloc(n * sizeof(int));
    if (!ped->child_offsets) return;

    int total_edges = 0;
    for (int i = 0; i < n; i++) {
        ped->child_offsets[i] = total_edges;
        total_edges += ped->n_children[i];
    }

    /* Allocate flat children array */
    ped->children = malloc(total_edges * sizeof(int));
    if (!ped->children) return;

    /* Fill children array */
    int* current_offset = calloc(n, sizeof(int));
    if (!current_offset) return;

    for (int i = 0; i < n; i++) {
        int fid = ped->individuals[i].father_id;
        int mid = ped->individuals[i].mother_id;
        if (fid >= 0) {
            int pos = ped->child_offsets[fid] + current_offset[fid];
            ped->children[pos] = i;
            current_offset[fid]++;
        }
        if (mid >= 0) {
            int pos = ped->child_offsets[mid] + current_offset[mid];
            ped->children[pos] = i;
            current_offset[mid]++;
        }
    }

    free(current_offset);
}

int topological_sort(pedigree* ped) {
    int n = ped->n_individuals;
    if (n == 0) return 0;

    /* Build child lists for traversal */
    build_child_lists(ped);
    if (!ped->children && n > 0 && ped->n_founders < n) {
        return -1; /* Allocation failed */
    }

    /* Allocate result array */
    free(ped->topo_order);
    ped->topo_order = malloc(n * sizeof(int));
    if (!ped->topo_order) return -1;

    /* Copy in_degree values (we'll decrement them) */
    int* in_degree = malloc(n * sizeof(int));
    if (!in_degree) return -1;

    for (int i = 0; i < n; i++) {
        in_degree[i] = ped->individuals[i].in_degree;
    }

    /* Queue for individuals with in_degree = 0 */
    int* queue = malloc(n * sizeof(int));
    if (!queue) {
        free(in_degree);
        return -1;
    }

    int queue_front = 0, queue_back = 0;

    /* Enqueue all founders */
    for (int i = 0; i < n; i++) {
        if (in_degree[i] == 0) {
            queue[queue_back++] = i;
        }
    }

    /* Kahn's algorithm */
    int order_idx = 0;
    while (queue_front < queue_back) {
        int i = queue[queue_front++];
        ped->topo_order[order_idx++] = i;

        /* Decrement in_degree of all children */
        for (int j = 0; j < ped->n_children[i]; j++) {
            int child = ped->children[ped->child_offsets[i] + j];
            in_degree[child]--;
            if (in_degree[child] == 0) {
                queue[queue_back++] = child;
            }
        }
    }

    free(in_degree);
    free(queue);

    /* Check if all individuals were processed */
    if (order_idx != n) {
        return -1; /* Cycle detected */
    }

    return 0;
}

/* ============================================================================
 * CHROMOSOME FUNCTIONS
 * ============================================================================ */

ped_chromosome* create_chromosome(int initial_capacity) {
    ped_chromosome* chr = malloc(sizeof(ped_chromosome));
    if (!chr) return NULL;

    chr->segments = malloc(initial_capacity * sizeof(ped_segment));
    if (!chr->segments) {
        free(chr);
        return NULL;
    }

    chr->n_segments = 0;
    chr->capacity = initial_capacity;
    return chr;
}

ped_chromosome* create_founder_chromosome(int founder_id, int homolog) {
    ped_chromosome* chr = create_chromosome(4);
    if (!chr) return NULL;

    founder_origin origin = {founder_id, homolog};
    chromosome_add_segment(chr, 0.0, 1.0, origin);
    return chr;
}

void free_chromosome(ped_chromosome* chr) {
    if (!chr) return;
    free(chr->segments);
    free(chr);
}

void chromosome_add_segment(ped_chromosome* chr, double start, double end,
                            founder_origin origin) {
    if (ensure_capacity_segments(chr) < 0) return;

    chr->segments[chr->n_segments].start = start;
    chr->segments[chr->n_segments].end = end;
    chr->segments[chr->n_segments].origin = origin;
    chr->n_segments++;
}

void copy_segments_in_range(ped_chromosome* dest, ped_chromosome* source,
                            double range_start, double range_end) {
    for (int i = 0; i < source->n_segments; i++) {
        ped_segment* seg = &source->segments[i];

        /* Skip segments entirely outside range */
        if (seg->end <= range_start || seg->start >= range_end) {
            continue;
        }

        /* Clip segment to range */
        double new_start = (seg->start > range_start) ? seg->start : range_start;
        double new_end = (seg->end < range_end) ? seg->end : range_end;

        chromosome_add_segment(dest, new_start, new_end, seg->origin);
    }
}

void merge_adjacent_segments(ped_chromosome* chr) {
    if (chr->n_segments <= 1) return;

    int write_idx = 0;
    for (int read_idx = 1; read_idx < chr->n_segments; read_idx++) {
        if (origin_equal(chr->segments[write_idx].origin,
                         chr->segments[read_idx].origin)) {
            /* Merge: extend current segment */
            chr->segments[write_idx].end = chr->segments[read_idx].end;
        } else {
            /* Move to next position */
            write_idx++;
            if (write_idx != read_idx) {
                chr->segments[write_idx] = chr->segments[read_idx];
            }
        }
    }
    chr->n_segments = write_idx + 1;
}

ped_chromosome* copy_chromosome(ped_chromosome* src) {
    if (!src) return NULL;

    ped_chromosome* dest = create_chromosome(src->n_segments + 4);
    if (!dest) return NULL;

    for (int i = 0; i < src->n_segments; i++) {
        chromosome_add_segment(dest, src->segments[i].start,
                               src->segments[i].end, src->segments[i].origin);
    }
    return dest;
}

/* ============================================================================
 * MEIOSIS AND SIMULATION FUNCTIONS
 * ============================================================================ */

ped_chromosome* meiosis(ped_chromosome* parent_pat, ped_chromosome* parent_mat,
                        double rec_rate, gsl_rng* rng) {
    /* Generate number of crossovers from Poisson distribution */
    unsigned int n_crossovers = gsl_ran_poisson(rng, rec_rate);

    /* Allocate breakpoints array: [0, crossover positions..., 1] */
    double* breakpoints = malloc((n_crossovers + 2) * sizeof(double));
    if (!breakpoints) return NULL;

    breakpoints[0] = 0.0;
    breakpoints[n_crossovers + 1] = 1.0;

    /* Generate uniform random crossover positions */
    for (unsigned int i = 0; i < n_crossovers; i++) {
        breakpoints[i + 1] = gsl_rng_uniform(rng);
    }

    /* Sort crossover positions */
    if (n_crossovers > 1) {
        qsort(breakpoints + 1, n_crossovers, sizeof(double), compare_double);
    }

    /* Choose starting chromosome (0 = paternal, 1 = maternal) */
    int current_source = gsl_rng_uniform_int(rng, 2);

    /* Create gamete chromosome */
    int estimated_segments = parent_pat->n_segments + parent_mat->n_segments + n_crossovers;
    ped_chromosome* gamete = create_chromosome(estimated_segments);
    if (!gamete) {
        free(breakpoints);
        return NULL;
    }

    /* Build gamete by copying segments from alternating parents */
    for (unsigned int i = 0; i <= n_crossovers; i++) {
        double seg_start = breakpoints[i];
        double seg_end = breakpoints[i + 1];

        ped_chromosome* source = (current_source == 0) ? parent_pat : parent_mat;
        copy_segments_in_range(gamete, source, seg_start, seg_end);

        /* Alternate at each crossover */
        current_source = 1 - current_source;
    }

    free(breakpoints);

    /* Merge adjacent segments with identical origin */
    merge_adjacent_segments(gamete);

    return gamete;
}

ped_simulation* create_simulation(void) {
    ped_simulation* sim = malloc(sizeof(ped_simulation));
    if (!sim) return NULL;

    sim->ped = NULL;
    sim->chromosomes = NULL;
    sim->rec_rate = 1.0;
    sim->rng = NULL;
    sim->seed = 0;

    return sim;
}

void free_simulation(ped_simulation* sim) {
    if (!sim) return;

    if (sim->chromosomes && sim->ped) {
        for (int i = 0; i < sim->ped->n_individuals; i++) {
            free_chromosome(sim->chromosomes[i].paternal);
            free_chromosome(sim->chromosomes[i].maternal);
        }
        free(sim->chromosomes);
    }

    free_pedigree(sim->ped);

    if (sim->rng) {
        gsl_rng_free(sim->rng);
    }

    free(sim);
}

ped_simulation* simulate_pedigree(const char* ped_file, double rec_rate,
                                   unsigned long seed) {
    ped_simulation* sim = create_simulation();
    if (!sim) return NULL;

    /* Parse pedigree */
    sim->ped = parse_pedigree_file(ped_file);
    if (!sim->ped) {
        free_simulation(sim);
        return NULL;
    }

    sim->rec_rate = rec_rate;

    /* Initialize RNG */
    gsl_rng_env_setup();
    sim->rng = gsl_rng_alloc(gsl_rng_default);
    if (!sim->rng) {
        free_simulation(sim);
        return NULL;
    }

    if (seed == 0) {
        seed = (unsigned long)time(NULL);
    }
    sim->seed = seed;
    gsl_rng_set(sim->rng, seed);

    /* Allocate chromosome arrays */
    int n = sim->ped->n_individuals;
    sim->chromosomes = calloc(n, sizeof(diploid_chromosomes));
    if (!sim->chromosomes) {
        free_simulation(sim);
        return NULL;
    }

    /* Track founder index for chromosome labeling */
    int founder_idx = 0;

    /* Process individuals in topological order */
    for (int i = 0; i < n; i++) {
        int id = sim->ped->topo_order[i];
        ped_indiv* indiv = &sim->ped->individuals[id];

        if (indiv->is_founder) {
            /* Founder: create unique chromosomes */
            sim->chromosomes[id].paternal = create_founder_chromosome(founder_idx, 0);
            sim->chromosomes[id].maternal = create_founder_chromosome(founder_idx, 1);
            if (!sim->chromosomes[id].paternal || !sim->chromosomes[id].maternal) {
                free_simulation(sim);
                return NULL;
            }
            founder_idx++;
        } else {
            /* Non-founder: perform meiosis from each parent */
            int fid = indiv->father_id;
            int mid = indiv->mother_id;

            /* Paternal chromosome from father's meiosis */
            if (fid >= 0) {
                sim->chromosomes[id].paternal = meiosis(
                    sim->chromosomes[fid].paternal,
                    sim->chromosomes[fid].maternal,
                    rec_rate, sim->rng);
            } else {
                /* Father unknown - shouldn't happen in valid pedigree */
                sim->chromosomes[id].paternal = NULL;
            }

            /* Maternal chromosome from mother's meiosis */
            if (mid >= 0) {
                sim->chromosomes[id].maternal = meiosis(
                    sim->chromosomes[mid].paternal,
                    sim->chromosomes[mid].maternal,
                    rec_rate, sim->rng);
            } else {
                /* Mother unknown - shouldn't happen in valid pedigree */
                sim->chromosomes[id].maternal = NULL;
            }

            if (!sim->chromosomes[id].paternal || !sim->chromosomes[id].maternal) {
                free_simulation(sim);
                return NULL;
            }
        }
    }

    return sim;
}

/* ============================================================================
 * OUTPUT FUNCTIONS
 * ============================================================================ */

static void print_chromosome(ped_simulation* sim, ped_chromosome* chr,
                             const char* label, FILE* out) {
    fprintf(out, "  %s:\n", label);
    for (int i = 0; i < chr->n_segments; i++) {
        ped_segment* seg = &chr->segments[i];
        int fid = seg->origin.founder_id;

        /* Find founder's name */
        const char* founder_name = NULL;
        int founder_count = 0;
        for (int j = 0; j < sim->ped->n_individuals && founder_count <= fid; j++) {
            if (sim->ped->individuals[j].is_founder) {
                if (founder_count == fid) {
                    founder_name = sim->ped->id_to_name[j];
                    break;
                }
                founder_count++;
            }
        }

        const char* homolog_str = (seg->origin.homolog == 0) ? "pat" : "mat";

        if (founder_name) {
            fprintf(out, "    %.6f %.6f %s:%s\n",
                    seg->start, seg->end, founder_name, homolog_str);
        } else {
            fprintf(out, "    %.6f %.6f F%d:%s\n",
                    seg->start, seg->end, fid, homolog_str);
        }
    }
}

void print_individual_segments(ped_simulation* sim, int indiv_id, FILE* out) {
    if (indiv_id < 0 || indiv_id >= sim->ped->n_individuals) return;

    const char* name = pedigree_get_name(sim->ped, indiv_id);
    fprintf(out, "Individual: %s\n", name ? name : "(unknown)");

    diploid_chromosomes* chr = &sim->chromosomes[indiv_id];
    if (chr->paternal) {
        print_chromosome(sim, chr->paternal, "Paternal", out);
    }
    if (chr->maternal) {
        print_chromosome(sim, chr->maternal, "Maternal", out);
    }
}

void print_simulation_header(ped_simulation* sim, FILE* out) {
    fprintf(out, "# pedtrans output\n");
    fprintf(out, "# RecRate: %.6f Seed: %lu\n", sim->rec_rate, sim->seed);
    fprintf(out, "# Individuals: %d Founders: %d\n",
            sim->ped->n_individuals, sim->ped->n_founders);

    /* List founders */
    fprintf(out, "# Founders:");
    for (int i = 0; i < sim->ped->n_individuals; i++) {
        if (sim->ped->individuals[i].is_founder) {
            fprintf(out, " %s", sim->ped->id_to_name[i]);
        }
    }
    fprintf(out, "\n#\n");
}

void print_all_segments(ped_simulation* sim, FILE* out) {
    print_simulation_header(sim, out);

    for (int i = 0; i < sim->ped->n_individuals; i++) {
        int id = sim->ped->topo_order[i];
        print_individual_segments(sim, id, out);
        fprintf(out, "\n");
    }
}
