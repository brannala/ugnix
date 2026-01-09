# Per-Region MRCA Tracking Optimization

## Problem

`TestMRCAForTargetRegions` takes 52.7% of simulation time in sparse mode. It's called after every coalescence event and scans all chromosomes/segments to check if each target region has reached MRCA.

Current complexity: O(coalescence_events × chromosomes × segments × regions)

## Proposed Solution

Track active segment count per target region incrementally, rather than scanning everything after each coalescence.

### Data Structure

```c
typedef struct {
    int* active_count;      /* active_count[i] = segments overlapping region i */
    int* reached_mrca;      /* reached_mrca[i] = 1 if region i has MRCA */
    int regions_remaining;  /* count of regions not yet at MRCA */
} region_mrca_state;
```

### Algorithm

1. **Initialization**: For each sample's initial chromosome, count segments overlapping each target region.

2. **On coalescence**: When two chromosomes merge:
   - For merged segments in target regions: decrement count (two become one)
   - Check if any region's count dropped to 1 → mark as MRCA

3. **On recombination**: When chromosome splits:
   - Left/right children may have different target region overlaps
   - Update counts accordingly (usually no change, just redistribution)

4. **MRCA check**: Simply check `regions_remaining == 0`
   - O(1) instead of O(chromosomes × segments × regions)

### Key Functions to Modify

1. `coalescence()` in coalescent.c (~line 560)
   - After `mergeChr()`, update region counts
   - Decrement for coalesced segments

2. `recombination()` in coalescent.c (~line 340)
   - Track which regions each child chromosome overlaps
   - Usually no count change (segments redistribute)

3. `TestMRCAForTargetRegions()` in coalescent.c (~line 1675)
   - Replace full scan with: `return state->regions_remaining == 0;`

4. New helper functions:
   - `init_region_mrca_state(chrsample*, target_region_set*)`
   - `update_region_counts_coalescence(chromosome* merged, chromosome* chr1, chromosome* chr2)`
   - `update_region_counts_recombination(chromosome* left, chromosome* right, chromosome* parent)`
   - `free_region_mrca_state()`

### Complexity After Optimization

- Per coalescence: O(segments_in_merged_chromosome × regions)
- Per recombination: O(segments_in_parent × regions)
- MRCA check: O(1)

### Estimated Speedup

- TestMRCAForTargetRegions: 52.7% → ~5-10% of time
- Overall: ~40-50% faster

### Edge Cases

1. Chromosome pruning: When a chromosome is pruned (no target overlap), no count changes needed
2. Segment merging in `combineIdentAdjAncSegs`: May need to update if segments spanning region boundaries are merged
3. Single region case: Simpler tracking, can optimize further

### Testing

1. Verify final MRCA detection matches current implementation
2. Compare mutation counts and positions (should be identical)
3. Profile to confirm speedup

## Date

January 2026
