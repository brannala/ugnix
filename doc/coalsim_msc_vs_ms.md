# coalsim_msc vs Hudson's ms: Comparison and Benchmarks

## Overview

This document compares our `coalsim_msc` implementation with Hudson's `ms` program, the classic coalescent simulator. Both implement the coalescent with recombination and migration, but use different parameterizations.

## Parameter Comparison

| Parameter | Hudson's ms | coalsim_msc | Conversion |
|-----------|-------------|-------------|------------|
| Time units | 4N₀ generations | Generations | ms_time × 4N₀ |
| Population size | N₀ (reference) | N (per population) | Same meaning |
| Coalescence rate | n(n-1)/(4N₀) | n(n-1)/(4N) | Equivalent (both diploid) |
| Migration rate | M = 4N₀m | m (per-lineage/gen) | m = M / (4N₀) |
| Recombination | ρ = 4N₀rL | r (per-chrom/gen) | r = ρ / (4N₀) |
| Mutation | θ = 4N₀μ | μ (per-chrom/gen) | μ = θ / (4N₀) |

### Migration Direction Convention

Both use **forward-time** migration conventions:

| Program | Notation | Meaning |
|---------|----------|---------|
| ms | `-m i j M` | M_ij = flow INTO pop i FROM pop j |
| coalsim_msc | `A->B:m` | Flow FROM A TO B |

Example equivalence:
```
ms: -m 1 2 0.01          # Forward flow: 2 → 1
msc: B->A:0.000025       # Forward flow: B → A (with N=100)
```

## Statistical Validation

TMRCA distributions were compared using 1000-5000 replicates per test. Both t-tests and Kolmogorov-Smirnov tests confirm equivalent distributions (p > 0.05).

| Test | ms Mean | msc Mean | Difference | p-value (t-test) |
|------|---------|----------|------------|------------------|
| Single population | 362.77 gen | 355.29 gen | 2.1% | 0.08 |
| Isolated populations | 478.80 gen | 467.73 gen | 2.3% | >0.05 |
| Symmetric migration | 896.09 gen | 914.14 gen | 2.0% | 0.20 |
| Asymmetric migration | 709.27 gen | 732.56 gen | 3.3% | >0.05 |
| Different pop sizes | 1198.12 gen | 1193.52 gen | 0.4% | >0.05 |

All differences are within sampling variance.

## Performance Comparison

### Batch Mode

coalsim_msc supports batch mode for running multiple replicates efficiently:

```bash
./coalsim_msc -r 1000 -L -q control.txt
```

Options:
- `-r N` : Run N replicates
- `-L` : Output TMRCA per replicate (like ms -L)
- `-q` : Quiet mode (suppress banner)

### Runtime: Small Chromosomes (low recombination)

Parameters: N=100, n=10, no recombination

| Replicates | ms | coalsim_msc | Ratio |
|------------|-----|-------------|-------|
| 1,000 | 0.004s | 0.005s | 1.3x |
| 10,000 | 0.033s | 0.041s | 1.25x |
| 50,000 | 0.154s | 0.202s | 1.31x |

**Result:** ms is ~1.3x faster for simple simulations.

### Runtime: Large Chromosomes (high recombination)

Parameters: N=100, n=20, varying ρ

| ρ (rho) | Chromosome Size* | ms | coalsim_msc | Winner |
|---------|------------------|-----|-------------|--------|
| 10 | 2.5 Mb | 0.001s | 0.002s | ms (1.6x) |
| 100 | 25 Mb | 0.005s | 0.005s | tie |
| 500 | 125 Mb | 0.119s | 0.114s | msc (1.04x) |
| 1000 | 250 Mb | 0.941s | 0.635s | **msc (1.5x)** |
| 2000 | 500 Mb | 5.354s | 2.885s | **msc (1.9x)** |

*Assuming human-like recombination rate (r = 10⁻⁸/bp) with N=100

**Result:** coalsim_msc is faster for large chromosomes with high recombination.

### Runtime: Sample Size Scaling at High Recombination

Parameters: N=100, ρ=1000, 2 replicates

| Samples | ms | coalsim_msc | Speedup |
|---------|-----|-------------|---------|
| 10 | 0.378s | 0.208s | **1.8x faster** |
| 20 | 0.644s | 0.387s | **1.7x faster** |
| 50 | 0.789s | 0.696s | **1.1x faster** |
| 100 | 1.026s | 0.963s | **1.1x faster** |

### Runtime: With Migration

Parameters: N=100, n=5+5, M=1.0 symmetric migration

| Replicates | ms | coalsim_msc | Ratio |
|------------|-----|-------------|-------|
| 1,000 | 0.006s | 0.007s | 1.2x |
| 10,000 | 0.070s | 0.064s | **0.91x** |

**Result:** coalsim_msc is slightly faster with migration.

## Chromosome Size Translation

The relationship between ρ and physical chromosome size:

```
ρ = 4 × N × r × L

Where:
  N = effective population size
  r = per-base recombination rate (~10⁻⁸ for humans)
  L = chromosome length in bp
```

### For N=100 (test parameters)

| ρ | Chromosome Size |
|---|-----------------|
| 10 | 2.5 Mb |
| 100 | 25 Mb |
| 500 | 125 Mb |
| 1000 | 250 Mb |
| 2000 | 500 Mb |

### For N=10,000 (human-like)

| ρ | Chromosome Size |
|---|-----------------|
| 100 | 250 kb |
| 1000 | 2.5 Mb |
| 10,000 | 25 Mb |
| 40,000 | 100 Mb |

## MSC Overhead Analysis

Compared to the original `coalsim`, `coalsim_msc` has additional overhead:

### Original coalsim - O(1) rate calculation:
```c
totRate = (noChrom*(noChrom-1)/2.0)*(1.0/(2.0*popSize)) + (recRate+mutRate)*ancLength;
```

### coalsim_msc - O(n × k) rate calculation:
```c
msc_count_lineages(sample, tree);  // O(n_samples × n_nodes) per event!
for (int i = 0; i < tree->n_nodes; i++) { ... }  // O(n_nodes)
```

### Overhead sources:
1. **`msc_count_lineages`**: O(n_samples × n_nodes) linear search called every event
2. **Per-population rates**: O(n_nodes) iteration even for single population
3. **Rate structure**: malloc/free each event
4. **Population ID tracking**: Each chromosome carries population_id

### Single-Population Fast Path (Implemented)

When `n_tips == 1` and no migration is configured, coalsim_msc automatically uses an optimized fast path (`msc_simulate_single_pop_fast`) that:
- Uses O(1) rate calculation like the original coalsim
- Skips population tracking and lineage counting
- Avoids per-event structure allocation

### Fast Path Performance (coalsim_msc / coalsim):

| ρ (rho) | coalsim | msc_fast | Ratio | vs ms |
|---------|---------|----------|-------|-------|
| 10 | 0.0045s | 0.0030s | **0.66x** | 1.39x |
| 50 | 0.0204s | 0.0184s | **0.90x** | 1.33x |
| 100 | 0.0508s | 0.0524s | 1.03x | 1.31x |
| 200 | 0.1596s | 0.1628s | 1.02x | 1.04x |
| 500 | 0.8073s | 0.9170s | 1.14x | 1.01x |

*Parameters: N=100, n=10, 50 replicates each*

The fast path performs within 3-14% of the original coalsim at medium-to-high recombination, and is actually faster at low recombination rates.

### Multi-Population Overhead (without fast path):

| Recombination | 10 samples | 20 samples |
|---------------|------------|------------|
| ρ=100 | 0.66x | 0.89x |
| ρ=500 | 1.71x | 1.85x |

At high recombination rates with multiple populations, the per-event overhead from `msc_count_lineages` becomes significant.

### Remaining optimization opportunities:
1. **Incremental lineage counts**: Update counts only after coal/migration, not every event
2. **Direct ID lookup**: Use array instead of linear search for node-by-ID
3. **Rate caching**: Reuse rate structure, only update changed values

## Summary

1. **Correctness:** Both programs produce statistically equivalent TMRCA distributions.

2. **Performance trade-offs:**
   - ms is faster (~1.3x) for small, simple simulations
   - coalsim_msc is faster (1.5-2x) for large chromosomes with high recombination
   - coalsim_msc is slightly faster with migration

3. **Use cases:**
   - Quick, small-scale simulations: ms may be faster
   - Genome-scale simulations with recombination: coalsim_msc is faster
   - Species tree / MSC simulations: coalsim_msc (ms doesn't support this)
   - Migration between ancestral populations: coalsim_msc only

4. **Batch mode:** Essential for performance - eliminates process startup overhead.

## Running the Comparison

```bash
# Full comparison script
python3 scripts/compare_ms_msc.py

# Quick runtime test
./coalsim_msc -r 1000 -L -q control.txt
./msdir/ms 10 1000 -L -t 1.0
```

## Files

- `scripts/compare_ms_msc.py` - Statistical comparison script
- `doc/msc_migration_notes.md` - Detailed migration parameterization
- `msdir/` - Hudson's ms source and binary
