# uGnix User Manual

**Evolutionary Genomics Simulation and Analysis Toolkit**

Version 1.0

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Simulation Programs](#simulation-programs)
   - [coalsim - Coalescent Simulator](#coalsim)
   - [coalsim_msc - Multispecies Coalescent](#coalsim_msc)
   - [pedsim - Pedigree Simulator](#pedsim)
   - [pedtrans - Pedigree Transmission](#pedtrans)
   - [pedsim_multipop - Multi-Population Pedigree](#pedsim_multipop)
5. [Integrated Pipelines](#integrated-pipelines)
   - [pedsim_seq - Pedigree Sequence Pipeline](#pedsim_seq)
   - [pedsim_vcf - Pedigree VCF Pipeline](#pedsim_vcf)
   - [pedsim_vcf_multipop - Multi-Population VCF Pipeline](#pedsim_vcf_multipop)
6. [Assembly Tools](#assembly-tools)
   - [seqassemble - Sequence Assembly](#seqassemble)
   - [vcfassemble - VCF Assembly](#vcfassemble)
7. [Analysis Tools](#analysis-tools)
   - [kinship - Relatedness Calculator](#kinship)
   - [het - Heterozygosity Analysis](#het)
   - [gsum - Genotype Summary](#gsum)
   - [hwe-dis - Hardy-Weinberg Test](#hwe-dis)
   - [sample - VCF Subsampling](#sample)
8. [File Formats](#file-formats)
9. [Workflows and Pipelines](#workflows-and-pipelines)
10. [Theoretical Background](#theoretical-background)
11. [Troubleshooting](#troubleshooting)

---

## Introduction

uGnix is a comprehensive toolkit for simulating and analyzing evolutionary genomic data. It provides:

- **Coalescent simulation** for generating gene genealogies and sequence data
- **Pedigree simulation** for modeling inheritance in finite populations
- **Multi-population models** with migration
- **Multispecies coalescent** for phylogenomic studies
- **Population genetics analysis** tools

The programs are designed to work together in pipelines, enabling complex simulation scenarios from simple building blocks.

---

## Installation

### Requirements

- GCC compiler
- GNU Scientific Library (GSL)
- GLib 2.0

### Building

```bash
git clone https://github.com/ugnix/ugnix.git
cd ugnix
make all
```

### Testing

```bash
make tests
./runtests
```

---

## Quick Start

### Simulate a simple coalescent dataset

```bash
# Generate 20 sequences with recombination and mutation
./coalsim -c 20 -N 10000 -r 1.0 -m 0.5 -V output.vcf
```

### Simulate sequences through a pedigree

```bash
# All-in-one pipeline
./pedsim_vcf -n 10 -N 1000 -g 5 -r 1.0 -m 0.5 -o output.vcf
```

### Multispecies coalescent with migration

```bash
# Create control file
cat > msc.ctl << EOF
species_tree: (A:1000#0.01,B:1000#0.01)ROOT#0.01;
samples: A:10,B:10
migration: A->B:1.0,B->A:1.0
mutation_rate: 0.5
output_vcf: true
vcf_file: species.vcf
EOF

./coalsim_msc msc.ctl
```

---

## Simulation Programs

<a name="coalsim"></a>
### coalsim - Coalescent Simulator

Simulates gene genealogies and sequences under the standard coalescent model for a single population.

#### Synopsis

```
coalsim [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-c NUM` | Sample size (haploid chromosomes) | Required |
| `-N NUM` | Effective population size | 1000 |
| `-r RATE` | Recombination rate (cM) | 0.05 |
| `-m RATE` | Mutation rate | 0.5 |
| `-s SEED` | Random seed | time-based |
| `-o FILE` | FASTA output file | - |
| `-V FILE` | VCF output file | - |
| `-u UNITS` | Sequence units (b, Kb, Mb) | Mb |
| `-M MODEL` | Substitution model (JC69, HKY) | JC69 |
| `-k KAPPA` | Ti/Tv ratio for HKY | 2.0 |
| `-p FREQS` | Base frequencies (piA,piC,piG,piT) | equal |
| `-T PARAMS` | Sparse target regions (n,len,start,gap) | - |
| `-R FILE` | Load target regions from file | - |
| `-a MODE` | MRCA output (r/i/s) | - |
| `-g MODE` | Gene tree output (s/f) | - |
| `-v` | Verbose output | off |
| `-q` | Quiet mode (no progress bar) | off |

#### Examples

**Basic simulation with VCF output:**
```bash
coalsim -c 100 -N 10000 -r 1.0 -m 0.5 -V variants.vcf
```

**HKY model with custom parameters:**
```bash
coalsim -c 50 -M HKY -k 4.0 -p 0.3,0.2,0.2,0.3 -o sequences.fa
```

**Sparse simulation (target regions only):**
```bash
coalsim -c 100 -T 1000,100,0,1000 -V sparse.vcf
```

---

<a name="coalsim_msc"></a>
### coalsim_msc - Multispecies Coalescent with Migration

Simulates gene genealogies under the multispecies coalescent model, with optional migration between populations.

#### Synopsis

```
coalsim_msc [options] <control_file>
```

#### Options

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Verbose output |
| `-h, --help` | Show help |

#### Control File Format

```
# Comments start with #
species_tree: <newick_with_theta>
samples: <species:count,...>
recombination_rate: <rate>
mutation_rate: <rate>
migration: <source->dest:M,...>
seed: <number>
output_gene_trees: true|false
output_vcf: true|false
vcf_file: <filename>
```

#### Species Tree Format

Extended Newick with theta (4Nu) annotations:

```
(children)NAME:branch_length#theta
```

- `NAME`: Optional name for internal nodes (required for ancestral migration)
- `branch_length`: Divergence time in generations
- `theta`: Population parameter 4Nu

#### Migration Format

```
source->dest:M
```

- `source`: Population where individuals originate (forward in time)
- `dest`: Population where individuals arrive
- `M`: Migration rate (4Nm)

Migration only occurs between contemporary populations.

#### Examples

**Two species with migration:**
```
species_tree: (A:1000#0.01,B:1000#0.01)ROOT#0.01;
samples: A:20,B:20
migration: A->B:1.0,B->A:0.5
mutation_rate: 0.5
output_vcf: true
vcf_file: two_species.vcf
seed: 12345
```

**Three species with ancestral migration:**
```
species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;
samples: A:10,B:10,C:10
migration: A->B:0.5,B->A:0.5,AB->C:0.1,C->AB:0.1
mutation_rate: 0.5
output_vcf: true
vcf_file: three_species.vcf
```

---

<a name="pedsim"></a>
### pedsim - Backwards Pedigree Simulator

Generates random pedigrees backwards in time from a sample of diploid individuals in a Wright-Fisher population.

#### Synopsis

```
pedsim [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-n NUM` | Sample size (individuals at present) | 10 |
| `-N NUM` | Population size (must be even) | 1000 |
| `-k NUM` | Generations to simulate backwards | 5 |
| `-s SEED` | Random seed | time-based |
| `-o FILE` | Output pedigree file | stdout |
| `-d FILE` | Output DOT file for visualization | - |
| `-S` | Print statistics to stderr | off |
| `-E` | Print expected values | off |
| `-v` | Verbose output | off |

#### Output Format

```
individualID fatherID motherID
```

Founders have fatherID = motherID = 0.

#### Example

```bash
# Generate pedigree with visualization
pedsim -n 20 -N 500 -k 10 -o pedigree.txt -d pedigree.dot -S

# Convert to image
dot -Tpng pedigree.dot -o pedigree.png
```

---

<a name="pedtrans"></a>
### pedtrans - Pedigree Chromosome Transmission

Simulates chromosome transmission through a pedigree, tracking founder origin for each genomic segment.

#### Synopsis

```
pedtrans [options] <pedigree_file>
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-c NUM` | Number of chromosomes | 1 |
| `-r RATE` | Recombination rate (crossovers/meiosis) | 1.0 |
| `-s SEED` | Random seed | time-based |
| `-o FILE` | Output file | stdout |
| `-v` | Verbose output | off |

#### Output Format

Chromosome segments with founder origin and positions normalized to [0, 1].

#### Example

```bash
pedtrans -c 22 -r 1.5 pedigree.txt -o segments.txt
```

---

<a name="pedsim_multipop"></a>
### pedsim_multipop - Multi-Population Pedigree Simulator

Extends pedsim to multiple populations with migration.

#### Synopsis

```
pedsim_multipop [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output pedigree file | Required |
| `-J NUM` | Number of populations | 1 |
| `-N LIST` | Population sizes (comma-separated) | 1000 |
| `-n LIST` | Sample sizes (comma-separated) | 10 |
| `-k NUM` | Generations | 5 |
| `-m STR` | Migration matrix ("m00,m01;m10,m11") | - |
| `-M NUM` | Symmetric island model rate | - |
| `-s SEED` | Random seed | time-based |
| `-d` | Output DOT format | off |
| `-v` | Verbose output | off |

#### Examples

**Two populations with island model:**
```bash
pedsim_multipop -J 2 -N 500,500 -n 10,10 -M 0.01 -k 10 -o pedigree.txt
```

**Custom migration matrix:**
```bash
pedsim_multipop -J 2 -N 500,500 -n 10,10 -m "0,0.02;0.01,0" -k 10 -o pedigree.txt
```

---

## Integrated Pipelines

<a name="pedsim_seq"></a>
### pedsim_seq - Pedigree Sequence Pipeline

Combines pedsim, coalsim, pedtrans, and seqassemble into a single tool.

#### Synopsis

```
pedsim_seq [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-n NUM` | Sample size | 10 |
| `-N NUM` | Effective population size | 1000 |
| `-k NUM` | Generations | 5 |
| `-r RATE` | Coalescent recombination rate | 0.1 |
| `-m RATE` | Mutation rate | 0.5 |
| `-R RATE` | Pedtrans recombination rate | 1.0 |
| `-o FILE` | Output file | stdout |
| `-f FORMAT` | Output format (fasta/vcf) | fasta |
| `-S` | Samples only (exclude founders) | off |
| `-s SEED` | Random seed | time-based |
| `-v` | Verbose output | off |

#### Example

```bash
pedsim_seq -n 50 -N 5000 -k 10 -r 1.0 -m 0.5 -f vcf -o output.vcf
```

---

<a name="pedsim_vcf"></a>
### pedsim_vcf - Pedigree VCF Pipeline

Streamlined pipeline for VCF output.

#### Synopsis

```
pedsim_vcf [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-g NUM` | Generations | 5 |
| `-n NUM` | Sample size | 10 |
| `-P NUM` | Pedigree population size | 1000 |
| `-N NUM` | Coalescent effective size | 10000 |
| `-r NUM` | Recombination rate (cM) | 1.0 |
| `-m NUM` | Mutation rate | 0.5 |
| `-s NUM` | Random seed | time-based |
| `-o FILE` | Output VCF file | Required |
| `-S` | Samples only | off |
| `-k` | Keep temporary files | off |
| `-v` | Verbose output | off |

#### Example

```bash
pedsim_vcf -g 10 -n 100 -N 10000 -r 1.5 -m 0.5 -o population.vcf
```

---

<a name="pedsim_vcf_multipop"></a>
### pedsim_vcf_multipop - Multi-Population VCF Pipeline

Complete multi-population simulation to VCF.

#### Synopsis

```
pedsim_vcf_multipop [options]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output VCF file | Required |
| `-J NUM` | Number of populations | 1 |
| `-N LIST` | Population sizes | 1000 |
| `-n LIST` | Sample sizes | 10 |
| `-g NUM` | Generations | 5 |
| `-c NUM` | Number of chromosomes (max 64) | 1 |
| `-m STR` | Migration matrix | - |
| `-M NUM` | Island model migration rate | - |
| `-E NUM` | Coalescent effective size | 10000 |
| `-r NUM` | Recombination rate (cM) | 1.0 |
| `-u NUM` | Mutation rate | 0.5 |
| `-s NUM` | Random seed | time-based |
| `-S` | Samples only | off |
| `-k` | Keep temporary files | off |
| `-v` | Verbose output | off |

#### Example

```bash
pedsim_vcf_multipop -J 3 -N 1000,500,500 -n 20,10,10 \
    -M 0.01 -g 10 -c 22 -o multipop.vcf
```

---

## Assembly Tools

<a name="seqassemble"></a>
### seqassemble - Sequence Assembly

Assembles sample sequences from pedigree segment data and founder sequences.

#### Synopsis

```
seqassemble [options] <segments_file> <founders_fasta>
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output file | stdout |
| `-f FORMAT` | Output format (fasta/vcf) | fasta |
| `-s` | Samples only | off |
| `-m` | Memory-efficient streaming mode | off |
| `-v` | Verbose output | off |

#### Example

```bash
seqassemble segments.txt founders.fa -f vcf -o assembled.vcf
```

---

<a name="vcfassemble"></a>
### vcfassemble - VCF Assembly

Assembles sample VCF from founder VCF and pedigree segments.

#### Synopsis

```
vcfassemble -v <founder.vcf> -p <segments.txt> -o <output.vcf> [-s]
```

#### Options

| Option | Description |
|--------|-------------|
| `-v FILE` | Founder VCF file (Required) |
| `-p FILE` | Pedtrans segment file (Required) |
| `-o FILE` | Output VCF file (Required) |
| `-s` | Samples only |

---

## Analysis Tools

<a name="kinship"></a>
### kinship - Relatedness Calculator

Calculates kinship, inbreeding, and relatedness coefficients from pedigrees.

#### Synopsis

```
kinship [options] <pedigree_file>
```

#### Options

| Option | Description |
|--------|-------------|
| `-i` | Calculate inbreeding coefficients |
| `-k` | Calculate pairwise kinship coefficients |
| `-r` | Calculate pairwise relatedness coefficients |
| `-v` | Verbose output |

#### Example

```bash
kinship -k -i pedigree.txt > coefficients.txt
```

---

<a name="het"></a>
### het - Heterozygosity Analysis

Calculates average heterozygosity for individuals and loci.

#### Synopsis

```
het [options] <genotype_file>
```

#### Options

| Option | Description |
|--------|-------------|
| `-i` | Per-individual heterozygosity |
| `-l` | Per-locus heterozygosity |

#### Output

Heterozygosity values with 95% confidence intervals. Missing data marked as "M".

---

<a name="gsum"></a>
### gsum - Genotype Summary

Provides summary statistics for genotype data files.

#### Synopsis

```
gsum [options] <genotype_file>
```

#### Options

| Option | Description |
|--------|-------------|
| `-i` | Summarize individuals |
| `-p` | Summarize populations |
| `-l` | Summarize loci |

---

<a name="hwe-dis"></a>
### hwe-dis - Hardy-Weinberg Equilibrium Test

Tests for departure from Hardy-Weinberg equilibrium.

#### Synopsis

```
hwe-dis <genotype_file>
```

#### Output

Disequilibrium coefficient (D) with 95% CI. Significant departures marked with asterisk (*).

---

<a name="sample"></a>
### sample - VCF Subsampling

Subsample individuals and/or markers from VCF files.

#### Synopsis

**Single-file mode:**
```
sample [options] <input.vcf>
```

**Multi-file mode (marker intersection):**
```
sample -O <suffix> -n <num> -m <method> <input1.vcf> <input2.vcf> [input3.vcf]
```

#### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output file | stdout |
| `-I FILE` | File with individual IDs to keep | - |
| `-n NUM` | Number of markers to select | - |
| `-m METHOD` | Selection method (u=uniform, r=random) | random |
| `-R REGION` | Region filter (CHROM or CHROM:START-END) | - |
| `-s SEED` | Random seed | time-based |
| `-O SUFFIX` | Output suffix (multi-file mode) | - |
| `-v` | Verbose output | off |

#### Examples

**Subsample individuals:**
```bash
sample -I keep_ids.txt input.vcf -o subset.vcf
```

**Random marker subsample:**
```bash
sample -n 1000 -m r input.vcf -o subset.vcf
```

**Uniformly spaced markers:**
```bash
sample -n 1000 -m u input.vcf -o subset.vcf
```

**Region extraction:**
```bash
sample -R chr1:1000000-2000000 input.vcf -o region.vcf
```

---

## File Formats

### Pedigree Format

Used by pedsim, pedtrans, kinship.

```
individualID fatherID motherID
```

- One individual per line
- Founders have fatherID = motherID = 0
- IDs are positive integers

**Example:**
```
1 0 0
2 0 0
3 1 2
4 1 2
5 3 4
```

### Genotype Format (BA3/Immanc)

Used by het, gsum, hwe-dis.

Population-structured genotype data with allele codes.

### VCF Format

Standard Variant Call Format for variant data. Used as output by coalsim, seqassemble, vcfassemble, and pipeline programs.

### FASTA Format

Standard biological sequence format. Used by coalsim and seqassemble.

```
>sequence_id
ACGTACGTACGT...
```

### Control File Format (coalsim_msc)

Key-value pairs with `#` comments:

```
# This is a comment
key: value
```

---

## Workflows and Pipelines

### Manual Pipeline: Pedigree to Sequences

```bash
# 1. Generate pedigree
pedsim -n 20 -N 1000 -k 5 -o pedigree.txt

# 2. Count founders
FOUNDERS=$(grep " 0 0$" pedigree.txt | wc -l)

# 3. Generate founder sequences (2 chromosomes per diploid)
coalsim -c $((2 * FOUNDERS)) -N 10000 -r 1.0 -m 0.5 -o founders.fa

# 4. Simulate transmission
pedtrans -r 1.0 pedigree.txt -o segments.txt

# 5. Assemble sample sequences
seqassemble segments.txt founders.fa -o samples.fa
```

### One-Command Pipeline

```bash
# Equivalent to manual pipeline above
pedsim_vcf -n 20 -P 1000 -N 10000 -g 5 -r 1.0 -m 0.5 -o samples.vcf
```

### Multi-Population Study

```bash
# Simulate 3 populations with island model migration
pedsim_vcf_multipop \
    -J 3 \
    -N 1000,1000,1000 \
    -n 50,50,50 \
    -M 0.01 \
    -g 20 \
    -c 22 \
    -E 10000 \
    -r 1.0 \
    -u 0.5 \
    -o three_pops.vcf
```

### Species Divergence Study

```bash
cat > species.ctl << 'EOF'
# Human-Chimp-Gorilla simulation
species_tree: ((Human:6000000#0.0001,Chimp:6000000#0.0001)HC:4000000#0.0001,Gorilla:10000000#0.0001)ROOT#0.0001;
samples: Human:50,Chimp:50,Gorilla:20
recombination_rate: 1.0
mutation_rate: 0.5
migration: Human->Chimp:0.1,Chimp->Human:0.1
seed: 42
output_vcf: true
vcf_file: primates.vcf
EOF

coalsim_msc species.ctl
```

---

## Theoretical Background

### Coalescent Model

The coalescent describes the genealogical history of a sample of genes, tracing lineages backwards in time until they coalesce to a common ancestor. Key parameters:

- **θ = 4Nu**: Population-scaled mutation rate
- **ρ = 4Nr**: Population-scaled recombination rate
- **Coalescence rate**: n(n-1)/(4N) for n lineages

### Wright-Fisher Model

The discrete-generation model underlying pedigree simulation:

- Fixed population size N (N/2 males, N/2 females)
- Random mating each generation
- Non-overlapping generations

### Multispecies Coalescent

Extends the coalescent to multiple species/populations:

- Lineages can only coalesce within the same population
- Populations merge at divergence times (going backwards)
- Migration allows lineage exchange between contemporary populations

### Migration (Isolation-with-Migration)

- **M = 4Nm**: Population-scaled migration rate
- Migration rate per lineage = M/2
- Higher M → more gene flow → less differentiation

---

## Troubleshooting

### Common Errors

**"Population size must be even"**
- pedsim requires even N for equal males/females
- Solution: Use an even number for -N

**"Species 'X' not found in tree"**
- Name mismatch in migration or samples specification
- Check spelling matches species tree exactly

**"Migration must come after species_tree"**
- Control file order matters for coalsim_msc
- Put species_tree line before migration line

**Simulation runs very long**
- Very small theta causes slow coalescence
- Many samples with recombination creates many fragments
- Try: increase theta, reduce samples, reduce recombination

**Memory errors with large simulations**
- Use streaming mode (-m) in seqassemble
- Reduce number of chromosomes or samples

### Performance Tips

1. **Use integrated pipelines** (pedsim_vcf) instead of manual steps
2. **Sparse simulation** (-T option in coalsim) for target regions only
3. **Streaming mode** (-m in seqassemble) for large datasets
4. **Appropriate theta values**: too small causes very deep trees

### Getting Help

- Use `-h` or `--help` with any program
- Check stderr output with `-v` for debugging
- Report issues at: https://github.com/ugnix/ugnix/issues

---

## Citation

If you use uGnix in your research, please cite:

[Citation information to be added]

---

## License

[License information to be added]

---

## Authors

[Author information to be added]
