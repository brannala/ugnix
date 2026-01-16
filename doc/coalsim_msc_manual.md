# coalsim_msc User Manual

Multispecies Coalescent Simulator with Migration

## Synopsis

```
coalsim_msc [options] <control_file>
```

## Description

`coalsim_msc` simulates gene genealogies under the multispecies coalescent model with optional migration between populations. The program reads a control file specifying the species tree, sample sizes, and simulation parameters, then generates coalescent histories consistent with the model.

## Options

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Print detailed progress information to stderr |
| `-h, --help` | Display help message and exit |

## Control File Format

The control file consists of key-value pairs, one per line. Comments begin with `#`.

### Required Parameters

#### species_tree

The species tree in extended Newick format with theta annotations.

**Format**: `(children)NAME:branch_length#theta`

- `children`: Comma-separated subtrees or tip names
- `NAME`: Optional name for internal nodes (required for ancestral migration)
- `branch_length`: Time to parent node in generations
- `theta`: Population parameter 4Nu for this branch

**Example**:
```
species_tree: ((A:1000#0.001,B:1000#0.002)AB:500#0.003,C:1500#0.001)ROOT#0.002;
```

This defines:
- Species A and B diverged 1000 generations ago (theta=0.001 and 0.002)
- Their ancestor AB existed for 500 generations (theta=0.003)
- Species C diverged 1500 generations ago (theta=0.001)
- Root population has theta=0.002

#### samples

Number of haploid samples from each species.

**Format**: `species:count,species:count,...`

**Example**:
```
samples: A:10,B:10,C:5
```

### Optional Parameters

#### recombination_rate

Recombination rate in cM (crossovers per generation per unit genetic distance). Default: 0.0

```
recombination_rate: 1.0
```

#### mutation_rate

Mutation rate for generating sequence variation. Default: 0.0

```
mutation_rate: 0.5
```

#### migration

Migration rates between populations using the M = 4Nm parameterization.

**Format**: `source->dest:M,source->dest:M,...`

- `source`: Population where individuals originate (forward in time)
- `dest`: Population where individuals arrive (forward in time)
- `M`: Migration rate (4Nm)

Migration only occurs between populations that exist at the same time. When populations merge at a divergence event, migration between them ceases.

**Example**:
```
migration: A->B:0.5,B->A:0.3
```

For ancestral migration, use internal node names:
```
migration: A->B:0.5,B->A:0.3,AB->C:0.1,C->AB:0.1
```

#### seed

Random number seed for reproducibility. Default: current time

```
seed: 12345
```

#### output_gene_trees

Output gene trees in Newick format. Values: `true` or `false`. Default: false

```
output_gene_trees: true
```

#### output_vcf

Generate VCF output file. Values: `true` or `false`. Default: false

```
output_vcf: true
```

#### vcf_file

Output filename for VCF (required if output_vcf is true).

```
vcf_file: output.vcf
```

#### verbose

Enable verbose output. Values: `true` or `false`. Default: false

```
verbose: true
```

## Examples

### Basic Two-Species Simulation

```
# Two species, no migration
species_tree: (A:500#0.01,B:500#0.01)#0.01;
samples: A:5,B:5
recombination_rate: 0.0
mutation_rate: 0.0
seed: 12345
```

### Two Species with Symmetric Migration

```
# Island model with equal migration
species_tree: (A:1000#0.01,B:1000#0.01)ROOT#0.01;
samples: A:10,B:10
recombination_rate: 0.0
mutation_rate: 1.0
migration: A->B:1.0,B->A:1.0
seed: 12345
output_vcf: true
vcf_file: symmetric_migration.vcf
```

### Three Species with Ancestral Migration

```
# ((A,B),C) tree with migration at multiple levels
species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;
samples: A:10,B:10,C:10
recombination_rate: 1.0
mutation_rate: 0.5
migration: A->B:0.5,B->A:0.3,AB->C:0.15,C->AB:0.1
seed: 42
output_vcf: true
vcf_file: three_species.vcf
```

Timeline for this example (going backwards in time):
- t=0 to t=300: Species A, B, C exist; A↔B migration active
- t=300 to t=500: Ancestral AB and C exist; AB↔C migration active
- t>500: Only ROOT exists; no migration

### Asymmetric Migration

```
# More migration from A to B than B to A
species_tree: (A:500#0.01,B:500#0.01)ROOT#0.01;
samples: A:10,B:10
migration: A->B:2.0,B->A:0.5
seed: 12345
```

### High Migration (Nearly Panmictic)

```
# Very high migration approximates a single population
species_tree: (A:1000#0.01,B:1000#0.01)ROOT#0.01;
samples: A:10,B:10
migration: A->B:100.0,B->A:100.0
seed: 12345
```

## Output

### Standard Output

The program prints a summary upon completion:

```
MSC Simulation Complete:
  Coalescences: 19
  Migrations: 156
  Recombinations: 12
  Mutations: 45
  Final time: 1523.67 generations
```

### Verbose Output (stderr)

With `-v` or `verbose: true`, detailed event information is printed:

```
Starting MSC simulation:
  Total samples: 20
  Species: 2
  Recombination rate: 1.0000
  Mutation rate: 0.5000
  Migration bands: 2
  Seed: 12345
...
Time 15.3: Migration - lineage from A to B
Time 23.7: Migration - lineage from B to A
Time 500.0: Divergence - merging populations 1 and 2 into 0
```

### VCF Output

When `output_vcf: true`, variants are written in standard VCF format.

## Parameter Guidelines

### Choosing Theta

Theta (θ = 4Nu) controls the coalescence rate within populations:
- Higher theta → slower coalescence → deeper gene trees
- Typical values: 0.001 to 0.1 for most species
- For humans: θ ≈ 0.001

### Choosing Migration Rates

M = 4Nm is the population-scaled migration rate:
- M = 0: Complete isolation
- M = 0.1-1: Low to moderate gene flow
- M = 1-10: Substantial gene flow
- M > 10: Populations behave nearly as one

Rule of thumb: M > 1 means more than one migrant per generation on average.

### Branch Lengths

Branch lengths are in generations:
- Should reflect actual divergence times
- Tree should be ultrametric (all tips at same time = present)

## Interpreting Results

### Migration Events

A high migration count indicates substantial gene flow. With symmetric migration M=1 between two populations with n samples each, expect roughly:

```
E[migrations] ≈ n × M × E[time_to_divergence]
```

### Effect on Gene Trees

- **No migration**: Gene tree topology constrained by species tree
- **Low migration**: Occasional discordance between gene and species trees
- **High migration**: Gene tree may differ substantially from species tree

### Coalescence Times

Migration generally reduces coalescence times because lineages can find common ancestors in either population. Very high migration makes the system behave like a single population with combined theta.

## Troubleshooting

### "Species 'X' not found in tree"

The species name in `samples` or `migration` doesn't match any node in the species tree. Check spelling and that internal node names are specified correctly.

### "Migration must come after species_tree"

The control file must specify `species_tree` before `migration` so population names can be validated.

### Simulation runs very long

This can happen with:
- Very small theta (slow coalescence)
- Many samples
- High recombination rate (creates many fragments)

Try reducing sample sizes or increasing theta for testing.

### No migration events occurring

Check that:
1. Migration populations are spelled correctly
2. Both populations exist at the simulation time
3. Migration rates are non-zero

## References

The multispecies coalescent model:
- Rannala B, Yang Z (2003) Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics 164:1645-1656

Isolation with migration:
- Hey J, Nielsen R (2004) Multilocus methods for estimating population sizes, migration rates and divergence time, with applications to the divergence of Drosophila pseudoobscura and D. persimilis. Genetics 167:747-760

## See Also

- `coalsim` - Single-population coalescent simulator
- Species tree estimation methods (ASTRAL, *BEAST, BPP)
