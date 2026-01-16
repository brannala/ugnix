/*
 * species_tree.c - Species tree parsing and operations for multispecies coalescent
 *
 * Parses extended Newick format with theta annotations:
 *   ((A:1000#0.001,B:1000#0.002):500#0.003,C:1500#0.001)#0.002;
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "species_tree.h"

/* ============== Internal Parser State ============== */

typedef struct {
    const char* input;
    const char* pos;
    int next_id;
    char error_msg[256];
} parser_state;

/* Forward declarations */
static species_node* parse_subtree(parser_state* state);
static void skip_whitespace(parser_state* state);
static int parse_name(parser_state* state, char* name, int max_len);
static double parse_number(parser_state* state);
static double parse_branch_length(parser_state* state);
static double parse_theta_annotation(parser_state* state);

/* ============== Parser Implementation ============== */

static void skip_whitespace(parser_state* state) {
    while (*state->pos && isspace(*state->pos)) {
        state->pos++;
    }
}

static int parse_name(parser_state* state, char* name, int max_len) {
    skip_whitespace(state);
    int i = 0;
    while (*state->pos && i < max_len - 1) {
        char c = *state->pos;
        /* Stop at Newick delimiters */
        if (c == ':' || c == '#' || c == ',' || c == ')' || c == '(' || c == ';') {
            break;
        }
        if (isspace(c)) {
            break;
        }
        name[i++] = c;
        state->pos++;
    }
    name[i] = '\0';
    return i;
}

static double parse_number(parser_state* state) {
    skip_whitespace(state);
    char* end;
    double val = strtod(state->pos, &end);
    if (end == state->pos) {
        snprintf(state->error_msg, sizeof(state->error_msg),
                 "Expected number at position %ld", state->pos - state->input);
        return -1;
    }
    state->pos = end;
    return val;
}

static double parse_branch_length(parser_state* state) {
    skip_whitespace(state);
    if (*state->pos != ':') {
        return 0.0;  /* No branch length specified */
    }
    state->pos++;  /* Skip ':' */
    return parse_number(state);
}

static double parse_theta_annotation(parser_state* state) {
    skip_whitespace(state);
    if (*state->pos != '#') {
        return -1.0;  /* No theta specified (error for MSC) */
    }
    state->pos++;  /* Skip '#' */
    return parse_number(state);
}

static species_node* create_node(parser_state* state) {
    species_node* node = calloc(1, sizeof(species_node));
    if (!node) return NULL;
    node->id = state->next_id++;
    node->theta = -1.0;  /* Mark as unset */
    return node;
}

/*
 * Parse a subtree (recursive descent)
 *
 * subtree -> leaf | internal
 * leaf -> name:length#theta
 * internal -> (subtree,subtree):length#theta
 */
static species_node* parse_subtree(parser_state* state) {
    skip_whitespace(state);

    species_node* node = create_node(state);
    if (!node) {
        snprintf(state->error_msg, sizeof(state->error_msg), "Memory allocation failed");
        return NULL;
    }

    if (*state->pos == '(') {
        /* Internal node */
        state->pos++;  /* Skip '(' */
        node->is_tip = 0;

        /* Parse left child */
        node->left = parse_subtree(state);
        if (!node->left) {
            free(node);
            return NULL;
        }
        node->left->parent = node;

        skip_whitespace(state);
        if (*state->pos != ',') {
            snprintf(state->error_msg, sizeof(state->error_msg),
                     "Expected ',' at position %ld", state->pos - state->input);
            species_tree_free((species_tree*)&node);  /* Cleanup partial tree */
            return NULL;
        }
        state->pos++;  /* Skip ',' */

        /* Parse right child */
        node->right = parse_subtree(state);
        if (!node->right) {
            free(node->left);
            free(node);
            return NULL;
        }
        node->right->parent = node;

        skip_whitespace(state);
        if (*state->pos != ')') {
            snprintf(state->error_msg, sizeof(state->error_msg),
                     "Expected ')' at position %ld", state->pos - state->input);
            free(node->left);
            free(node->right);
            free(node);
            return NULL;
        }
        state->pos++;  /* Skip ')' */

        /* Parse optional internal node name (e.g., ")AB:" or ")ROOT#") */
        skip_whitespace(state);
        if (isalpha(*state->pos)) {
            parse_name(state, node->name, MSC_MAX_SPECIES_NAME);
        } else {
            node->name[0] = '\0';
        }
    } else {
        /* Tip node - parse name */
        node->is_tip = 1;
        parse_name(state, node->name, MSC_MAX_SPECIES_NAME);
        if (node->name[0] == '\0') {
            snprintf(state->error_msg, sizeof(state->error_msg),
                     "Expected species name at position %ld", state->pos - state->input);
            free(node);
            return NULL;
        }
    }

    /* Parse branch length (optional) */
    node->branch_length = parse_branch_length(state);

    /* Parse theta annotation (required for MSC) */
    node->theta = parse_theta_annotation(state);

    return node;
}

/* Count nodes in tree (recursive) */
static int count_nodes(species_node* node) {
    if (!node) return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

/* Count tips in tree (recursive) */
static int count_tips(species_node* node) {
    if (!node) return 0;
    if (node->is_tip) return 1;
    return count_tips(node->left) + count_tips(node->right);
}

/* Collect nodes into array (recursive) */
static void collect_nodes(species_node* node, species_node** nodes, int* idx) {
    if (!node) return;
    nodes[(*idx)++] = node;
    collect_nodes(node->left, nodes, idx);
    collect_nodes(node->right, nodes, idx);
}

/* Collect tips into array (recursive) */
static void collect_tips(species_node* node, species_node** tips, int* idx) {
    if (!node) return;
    if (node->is_tip) {
        tips[(*idx)++] = node;
    }
    collect_tips(node->left, tips, idx);
    collect_tips(node->right, tips, idx);
}

species_tree* species_tree_parse_newick(const char* newick) {
    if (!newick || !*newick) {
        fprintf(stderr, "Error: Empty Newick string\n");
        return NULL;
    }

    parser_state state = {
        .input = newick,
        .pos = newick,
        .next_id = 0,
        .error_msg = ""
    };

    species_node* root = parse_subtree(&state);
    if (!root) {
        fprintf(stderr, "Error parsing species tree: %s\n", state.error_msg);
        return NULL;
    }

    /* Skip optional trailing semicolon */
    skip_whitespace(&state);
    if (*state.pos == ';') {
        state.pos++;
    }

    /* Check for trailing garbage */
    skip_whitespace(&state);
    if (*state.pos != '\0') {
        fprintf(stderr, "Warning: Unexpected characters after tree: '%s'\n", state.pos);
    }

    /* Build species_tree structure */
    species_tree* tree = calloc(1, sizeof(species_tree));
    if (!tree) {
        fprintf(stderr, "Error: Memory allocation failed for species tree\n");
        return NULL;
    }

    tree->root = root;
    tree->n_nodes = count_nodes(root);
    tree->n_tips = count_tips(root);

    /* Allocate and fill node arrays */
    tree->nodes = malloc(tree->n_nodes * sizeof(species_node*));
    tree->tips = malloc(tree->n_tips * sizeof(species_node*));

    if (!tree->nodes || !tree->tips) {
        fprintf(stderr, "Error: Memory allocation failed for node arrays\n");
        species_tree_free(tree);
        return NULL;
    }

    int idx = 0;
    collect_nodes(root, tree->nodes, &idx);
    idx = 0;
    collect_tips(root, tree->tips, &idx);

    return tree;
}

/* ============== Divergence Time Computation ============== */

/*
 * Compute divergence times recursively.
 * For tips: divergence_time = 0 (present)
 * For internal: divergence_time = max(child times) + branch_length of children to this node
 *
 * Actually, more precisely:
 * - Tips are at time 0
 * - Internal node time = child_time + child_branch_length (same for both children in ultrametric tree)
 */
static double compute_time_recursive(species_node* node) {
    if (!node) return 0;

    if (node->is_tip) {
        node->divergence_time = 0.0;
        return node->branch_length;  /* Return time to parent */
    }

    /* Internal node */
    double left_time = compute_time_recursive(node->left);
    double right_time = compute_time_recursive(node->right);

    /* In ultrametric tree, these should be equal */
    /* Use max for robustness */
    double child_time = (left_time > right_time) ? left_time : right_time;

    node->divergence_time = child_time;
    return child_time + node->branch_length;
}

void species_tree_compute_divergence_times(species_tree* tree) {
    if (!tree || !tree->root) return;
    compute_time_recursive(tree->root);
    /* Root divergence time is when all lineages merge */
    /* It's already set by the recursive call */
}

/* ============== Divergence Event List ============== */

/* Insert event in sorted order (by time ascending) */
static void insert_event_sorted(divergence_event_list* list, divergence_event* event) {
    if (!list->head || event->time < list->head->time) {
        event->next = list->head;
        list->head = event;
    } else {
        divergence_event* curr = list->head;
        while (curr->next && curr->next->time < event->time) {
            curr = curr->next;
        }
        event->next = curr->next;
        curr->next = event;
    }
    list->n_events++;
}

/* Create divergence events recursively */
static void create_events_recursive(species_node* node, divergence_event_list* list) {
    if (!node || node->is_tip) return;

    /* Create event for this internal node */
    divergence_event* event = malloc(sizeof(divergence_event));
    if (!event) return;

    event->time = node->divergence_time;
    event->child1 = node->left;
    event->child2 = node->right;
    event->parent = node;
    event->next = NULL;

    insert_event_sorted(list, event);

    /* Recurse to children */
    create_events_recursive(node->left, list);
    create_events_recursive(node->right, list);
}

divergence_event_list* species_tree_create_divergence_events(species_tree* tree) {
    if (!tree) return NULL;

    divergence_event_list* list = calloc(1, sizeof(divergence_event_list));
    if (!list) return NULL;

    create_events_recursive(tree->root, list);
    return list;
}

/* ============== Tree Operations ============== */

species_node* species_tree_find_species(species_tree* tree, const char* name) {
    if (!tree || !name) return NULL;

    for (int i = 0; i < tree->n_tips; i++) {
        if (strcmp(tree->tips[i]->name, name) == 0) {
            return tree->tips[i];
        }
    }
    return NULL;
}

species_node* species_tree_find_node_by_name(species_tree* tree, const char* name) {
    if (!tree || !name || !*name) return NULL;

    for (int i = 0; i < tree->n_nodes; i++) {
        if (tree->nodes[i]->name[0] && strcmp(tree->nodes[i]->name, name) == 0) {
            return tree->nodes[i];
        }
    }
    return NULL;
}

int species_tree_total_samples(species_tree* tree) {
    if (!tree) return 0;

    int total = 0;
    for (int i = 0; i < tree->n_tips; i++) {
        total += tree->tips[i]->sample_size;
    }
    return total;
}

/* ============== Printing ============== */

static void print_node_recursive(species_node* node, FILE* out, int depth) {
    if (!node) return;

    for (int i = 0; i < depth; i++) fprintf(out, "  ");

    if (node->is_tip) {
        fprintf(out, "Tip[%d]: %s, theta=%.6f, branch=%.1f, samples=%d\n",
                node->id, node->name, node->theta, node->branch_length, node->sample_size);
    } else {
        if (node->name[0]) {
            fprintf(out, "Internal[%d]: %s, theta=%.6f, branch=%.1f, div_time=%.1f\n",
                    node->id, node->name, node->theta, node->branch_length, node->divergence_time);
        } else {
            fprintf(out, "Internal[%d]: theta=%.6f, branch=%.1f, div_time=%.1f\n",
                    node->id, node->theta, node->branch_length, node->divergence_time);
        }
    }

    print_node_recursive(node->left, out, depth + 1);
    print_node_recursive(node->right, out, depth + 1);
}

void species_tree_print(species_tree* tree, FILE* out) {
    if (!tree) {
        fprintf(out, "NULL tree\n");
        return;
    }

    fprintf(out, "Species tree: %d nodes, %d tips\n", tree->n_nodes, tree->n_tips);
    print_node_recursive(tree->root, out, 0);
}

/* Print in Newick format (recursive helper) */
static void print_newick_recursive(species_node* node, FILE* out) {
    if (!node) return;

    if (node->is_tip) {
        fprintf(out, "%s", node->name);
    } else {
        fprintf(out, "(");
        print_newick_recursive(node->left, out);
        fprintf(out, ",");
        print_newick_recursive(node->right, out);
        fprintf(out, ")");
        /* Print internal node name if present */
        if (node->name[0]) {
            fprintf(out, "%s", node->name);
        }
    }

    if (node->branch_length > 0) {
        fprintf(out, ":%.1f", node->branch_length);
    }
    if (node->theta >= 0) {
        fprintf(out, "#%.6f", node->theta);
    }
}

void species_tree_print_newick(species_tree* tree, FILE* out) {
    if (!tree || !tree->root) return;
    print_newick_recursive(tree->root, out);
    fprintf(out, ";\n");
}

/* ============== Memory Management ============== */

static void free_node_recursive(species_node* node) {
    if (!node) return;
    free_node_recursive(node->left);
    free_node_recursive(node->right);
    free(node);
}

void species_tree_free(species_tree* tree) {
    if (!tree) return;
    free_node_recursive(tree->root);
    free(tree->nodes);
    free(tree->tips);
    free(tree);
}

void divergence_event_list_free(divergence_event_list* list) {
    if (!list) return;
    divergence_event* curr = list->head;
    while (curr) {
        divergence_event* next = curr->next;
        free(curr);
        curr = next;
    }
    free(list);
}
