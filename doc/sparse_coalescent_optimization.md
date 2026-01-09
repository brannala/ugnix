# Sparse Coalescent Simulation Optimization

## Overview

This document describes optimizations to `coalsim` for simulating coalescent genealogies when only specific target regions of a chromosome are of interest (sparse simulation). These optimizations dramatically improve performance for moderate chromosome lengths (r <= 1-2 cM) but have diminishing returns for very long chromosomes.

## Problem Statement

When simulating genetic data for hybrid classification (e.g., with mongrail2), we often need markers spread across a chromosome but don't need the full genealogy everywhere. The standard coalescent simulation tracks ancestry across the entire chromosome, which becomes expensive for long chromosomes due to:

1. Many recombination events creating new chromosomes
2. Large numbers of ancestry segments to track
3. Expensive bitarray operations for each segment

## Target Region Specification

Target regions are specified with the `-T` option:
```
coalsim -T n,length,start,spacing ...
```
- `n`: Number of target regions
- `length`: Length of each region (in [0,1] scale, e.g., 0.01 = 1% of chromosome)
- `start`: Start position of first region
- `spacing`: Gap between regions (0 = evenly distributed)

Example: `-T 10,0.01,0,0` creates 10 regions of length 0.01 evenly spaced across [0,1].

## Optimizations Implemented

### 1. MRCA Check Optimization

**Change:** Only check for MRCA after coalescence events, not after recombination or mutation events.

**Rationale:** MRCA status can only change when lineages coalesce. Checking after every event was wasteful.

**Impact:** ~2x reduction in MRCA checks.

### 2. Sparse Ancestry Tracking

**Change:** Skip bitarray operations (union, copy) for ancestry segments that don't overlap any target region.

**Implementation:**
- `set_sparse_target_regions()` - Sets global target regions for sparse tracking
- `segment_overlaps_targets()` - Checks if a segment overlaps any target
- Modified `mergeChr()` - Sets `abits = NULL` for non-target segments
- Modified `copy_chrom()` - Skips bitarray copy for non-target segments
- Modified `unionAnc()` - Handles NULL inputs gracefully

**Rationale:** Segments outside target regions don't contribute to the genealogy of interest. We still track segment positions (for chromosome structure) but skip expensive bitarray operations.

**Impact:** ~50x reduction in bitarray allocations for 10% target coverage.

### 3. Chromosome Pruning

**Change:** Remove chromosomes that have no segments overlapping any target region.

**Implementation:**
- `chromosome_overlaps_targets()` - Checks if any segment overlaps targets
- Modified `coalescence()` - Prunes merged chromosome if no target overlap
- Modified `recombination()` - Prunes resulting chromosomes if no target overlap

**Rationale:** Chromosomes with no target-overlapping segments cannot contribute to the genealogy of target regions and can be safely removed.

**Impact:** Keeps chromosome count low, reducing coalescence rate calculations.

### 4. Mutation Filtering

**Change:** Discard mutations that occur outside target regions.

**Implementation:** In coalsim.c main loop, `position_in_target_regions()` filters mutations.

**Rationale:** Only mutations in target regions are needed for downstream analysis.

## Why Single-Region Chromosomes Cannot Be Pruned

A natural question is whether chromosomes overlapping only ONE target region can also be pruned. This is **incorrect** because:

1. Such chromosomes carry active ancestry for samples in that target region
2. They must coalesce with other chromosomes for that region to reach MRCA
3. Pruning them would lose ancestry and produce incorrect genealogies

The key distinction:
- **No target overlap** → Can prune (doesn't contribute to any target genealogy)
- **One target overlap** → Must keep (contributes to that region's genealogy)
- **Multiple target overlaps** → Must keep (connects regions, determines when they become independent)

## Performance Results

### r=1.0 (10 target regions, 10% coverage)

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Total events | 742,862 | 77,089 | ~10x fewer |
| Bitarray allocs | 544M | 10.8M | ~50x fewer |
| Time | ~60s | ~1s | ~60x faster |

### r=10.0 (10 target regions, 10% coverage)

The optimizations help but the simulation remains slow:
- Still ~200 chromosomes active (spanning target regions)
- Must process recombination events to determine when regions separate
- ~770K events in 60s, TMRCA not yet reached

## Fundamental Limitation

For long chromosomes (r >> 1), performance is limited by the number of recombination events needed to separate target regions. This is inherent to the biology:

1. Target regions start connected on each sample's chromosome
2. Recombination events (anywhere on chromosome) can separate them
3. Many events are needed before all regions become independent
4. We cannot skip these events without changing what we're simulating

The recombination events in non-target regions are **necessary** because they determine when target regions become independent. Two regions remain correlated until recombination separates them.

## Possible Future Optimizations

1. **Per-region tracking:** Once a region becomes independent (no chromosomes span it and another region), simulate it separately as a smaller coalescent.

2. **Connectivity tracking:** Track only which target regions are connected, not full segment details in non-target regions. Complex to implement correctly.

3. **Approximate methods:** For very long chromosomes, consider SMC' or other approximations that assume regions become independent quickly.

## Files Modified

- `src/coalescent.c` - Sparse ancestry tracking, chromosome pruning
- `src/coalsim.c` - Target region setup, mutation filtering, MRCA check optimization
- `include/coalescent.h` - New function declarations

## Usage Example

```bash
# Sparse simulation with 10 target regions of 1% each
./coalsim -N 10000 -c 50 -r 1.0 -T 10,0.01,0,0 -m 0.01 -V output.vcf

# Verbose mode shows progress and pool statistics
./coalsim -N 10000 -c 50 -r 1.0 -T 10,0.01,0,0 -m 0.01 -v -V output.vcf
```

## Date

January 2026
