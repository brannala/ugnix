# MSC-M: Multispecies Coalescent with Migration

## Overview

This document describes the implementation of migration in the multispecies coalescent (MSC) simulator. The extension adds an Isolation-with-Migration (IM) model where lineages can migrate between contemporary populations before their divergence.

## Model Description

### Biological Interpretation (backwards in time)

In coalescent simulation, we trace lineages backwards in time. Migration events represent a lineage in population B having an ancestor that was in population A. The migration rate M = 4Nm is the population-scaled migration rate, where N is the effective population size and m is the per-generation migration probability.

For a lineage in population B:
- Migration rate to population A = M_AB / 2 per lineage
- This represents (forwards in time) individuals moving from A to B

### Key Properties

1. **Asymmetric migration**: Different rates A→B vs B→A are supported
2. **Contemporary populations only**: Migration only occurs between populations that exist at the current simulation time
3. **Stops at divergence**: When populations merge (going backwards), migration between them ceases
4. **Ancestral migration**: Internal nodes can be named to specify migration between ancestral populations

## Implementation Details

### Data Structures

```c
// Single migration band (in include/msc.h)
typedef struct {
    int from_pop;    // Source population ID (backwards: where lineage goes)
    int to_pop;      // Dest population ID (backwards: where lineage comes from)
    double M;        // Migration rate M = 4Nm
} migration_band;

// List of migration bands
typedef struct {
    migration_band* bands;
    int n_bands;
} migration_band_list;

// Per-band rate during simulation
typedef struct {
    int from_pop;
    int to_pop;
    double rate;     // n_to * M / 2
} mig_rate_info;
```

### Rate Calculation

For each migration band with rate M from population i to j:
```
migration_rate = n_j × M / 2
```
where n_j is the number of lineages currently in population j.

Total event rate:
```
total_rate = coal_rate + mig_rate + rec_rate + mut_rate
```

### Event Selection

Events are selected proportionally to their rates:
1. If u < coal_rate → coalescence event
2. Else if u < coal_rate + mig_rate → migration event
3. Else if u < coal_rate + mig_rate + rec_rate → recombination event
4. Else → mutation event

### Population Existence Logic

A population exists at time t if:
- **Tips**: Always exist from t=0 until parent's divergence_time
- **Internal nodes**: Exist from their divergence_time until parent's divergence_time (or forever for root)

Migration between populations A and B is only valid when both exist at current simulation time.

### Files Modified

| File | Changes |
|------|---------|
| `include/msc.h` | Added migration_band, migration_band_list, mig_rate_info structs; extended msc_params and msc_rates |
| `src/msc.c` | Added parse_migrations(), migration rate calculation, migration event handling, helper functions |
| `include/species_tree.h` | Added species_tree_find_node_by_name() declaration |
| `src/species_tree.c` | Extended Newick parser for internal node names; added species_tree_find_node_by_name() |
| `src/msc_main.c` | Updated help text |

### Key Functions

```c
// Parse migration specification "A->B:0.5,B->A:0.3"
static int parse_migrations(const char* str, msc_params* params);

// Check if population exists at given time
int msc_population_exists_at_time(species_node* node, double current_time);

// Check if two populations are contemporary
static int populations_contemporary(species_tree* tree, int pop1_id, int pop2_id, double current_time);

// Select migration band weighted by rate
int msc_select_migration_pair(gsl_rng* r, msc_rates* rates);

// Select random lineage from a population
int msc_select_lineage_in_pop(gsl_rng* r, chrsample* sample, int pop_id);
```

## Control File Format

### Extended Newick Format

Internal nodes can now have names:
```
(children)NAME:branch_length#theta
```

Example:
```
((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;
```

### Migration Specification

```
migration: source->dest:M,source->dest:M,...
```

Where:
- `source`, `dest` are population names (tips or internal nodes)
- `M` is the migration rate (4Nm)

### Complete Example

```
# Three species with tip and ancestral migration
species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;
samples: A:10,B:10,C:10
recombination_rate: 1.0
mutation_rate: 0.5
migration: A->B:0.5,B->A:0.3,AB->C:0.15,C->AB:0.1
seed: 12345
```

Timeline (going backwards):
- t=0 to t=300: A, B, C exist → A↔B migration active
- t=300 to t=500: AB, C exist → AB↔C migration active
- t>500: Only ROOT exists → no migration

## Testing

Test file: `testing/test_msc.c` (20 tests)

### Test Categories

1. **Newick parser tests**: Internal node name parsing
2. **Migration parsing tests**: Control file parsing
3. **Population existence tests**: Time-dependent population logic
4. **Rate calculation tests**: Migration rate formula verification
5. **Migration selection tests**: Weighted random selection
6. **Full simulation tests**: End-to-end verification
7. **Memory management tests**: Cleanup without leaks

### Running Tests

```bash
make test_msc
./test_msc
```

### Key Regression Tests

- `test_simulation_no_migration_completes`: MSC without migration still works
- `test_simulation_zero_migration_matches_original`: M=0 behaves same as no migration
- `test_migration_rate_zero_after_divergence`: Migration stops at divergence time

## Usage Notes

### High Migration Rates

With very high M (e.g., M=100), populations behave nearly as a single panmictic population. The TMRCA distribution approaches that of a single-population coalescent.

### Asymmetric Migration

Higher A→B migration (going forward) results in more lineages coalescing in population A (going backward), as lineages tend to accumulate in the source population.

### Performance Considerations

Migration events are relatively cheap (just changing a population_id), but high migration rates generate many events. For M=100 with 8 samples over 1000 generations, expect ~60,000+ migration events.

## Future Extensions

Potential enhancements:
1. Time-varying migration rates
2. Migration matrices for >2 populations
3. Stepping-stone models
4. Migration rate inference from data
