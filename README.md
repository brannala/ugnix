# uGnix

**Evolutionary Genomics Simulation and Analysis Toolkit**

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]()

uGnix is a comprehensive C toolkit for simulating and analyzing evolutionary genomic data. It provides coalescent simulation, pedigree-based inheritance modeling, multi-population dynamics with migration, and population genetics analysis tools.

---

## Features

- **Coalescent Simulation** - Generate gene genealogies and sequence data under the standard coalescent
- **Multispecies Coalescent** - Simulate across species trees with migration (Isolation-with-Migration model)
- **Pedigree Simulation** - Model inheritance in finite Wright-Fisher populations
- **Multi-Population Models** - Support for multiple populations with custom migration matrices
- **Integrated Pipelines** - Single-command workflows from pedigree to VCF
- **Analysis Tools** - Heterozygosity, kinship, Hardy-Weinberg tests, and more

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Programs](#programs)
  - [Simulation Programs](#simulation-programs)
  - [Integrated Pipelines](#integrated-pipelines)
  - [Assembly Tools](#assembly-tools)
  - [Analysis Tools](#analysis-tools)
- [Detailed Usage](#detailed-usage)
  - [coalsim](#coalsim---coalescent-simulator)
  - [coalsim_msc](#coalsim_msc---multispecies-coalescent-with-migration)
  - [pedsim](#pedsim---backwards-pedigree-simulator)
  - [pedtrans](#pedtrans---pedigree-chromosome-transmission)
  - [pedsim_multipop](#pedsim_multipop---multi-population-pedigree-simulator)
  - [pedsim_seq](#pedsim_seq---pedigree-sequence-pipeline)
  - [pedsim_vcf](#pedsim_vcf---pedigree-vcf-pipeline)
  - [pedsim_vcf_multipop](#pedsim_vcf_multipop---multi-population-vcf-pipeline)
  - [seqassemble](#seqassemble---sequence-assembly)
  - [vcfassemble](#vcfassemble---vcf-assembly)
  - [kinship](#kinship---relatedness-calculator)
  - [het](#het---heterozygosity-analysis)
  - [gsum](#gsum---genotype-summary)
  - [hwe-dis](#hwe-dis---hardy-weinberg-equilibrium-test)
  - [sample](#sample---vcf-subsampling)
- [File Formats](#file-formats)
- [Workflows](#workflows)
- [Theoretical Background](#theoretical-background)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

---

## Installation

### Requirements

- GCC compiler
- GNU Scientific Library (GSL)
- GLib 2.0

**Ubuntu/Debian:**
```bash
sudo apt-get install build-essential libgsl-dev libglib2.0-dev
```

**macOS (Homebrew):**
```bash
brew install gsl glib
```

### Building

```bash
git clone https://github.com/ugnix/ugnix.git
cd ugnix
make all
```

### Testing

```bash
make test_msc test_coalescent test_ugnix
./test_msc && ./test_coalescent && ./test_ugnix
```

### Install (optional)

```bash
sudo cp coalsim coalsim_msc pedsim pedtrans pedsim_vcf /usr/local/bin/
```

---

## Quick Start

### Simulate a coalescent dataset

```bash
# Generate 20 sequences with recombination and mutation
./coalsim -c 20 -N 10000 -r 1.0 -m 0.5 -V output.vcf
```

### Simulate sequences through a pedigree

```bash
# All-in-one pipeline: pedigree → coalescent → VCF
./pedsim_vcf -n 10 -N 1000 -g 5 -r 1.0 -m 0.5 -o output.vcf
```

### Multispecies coalescent with migration

```bash
cat > species.ctl << 'EOF'
species_tree: (A:1000#0.01,B:1000#0.01)ROOT#0.01;
samples: A:10,B:10
migration: A->B:1.0,B->A:1.0
mutation_rate: 0.5
output_vcf: true
vcf_file: species.vcf
EOF

./coalsim_msc species.ctl
```

### Multi-population simulation

```bash
# 3 populations with island-model migration
./pedsim_vcf_multipop -J 3 -N 1000,1000,1000 -n 20,20,20 \
    -M 0.01 -g 10 -o multipop.vcf
```

---

## Programs

### Simulation Programs

| Program | Description |
|---------|-------------|
| `coalsim` | Single-population coalescent simulator with recombination and mutation |
| `coalsim_msc` | Multispecies coalescent with migration between populations |
| `pedsim` | Backwards-in-time pedigree simulator (Wright-Fisher model) |
| `pedtrans` | Chromosome transmission simulator tracking founder segments |
| `pedsim_multipop` | Multi-population pedigree simulator with migration |

### Integrated Pipelines

| Program | Description |
|---------|-------------|
| `pedsim_seq` | Complete pedigree-to-sequence pipeline (FASTA/VCF output) |
| `pedsim_vcf` | Streamlined pedigree-to-VCF pipeline |
| `pedsim_vcf_multipop` | Multi-population pedigree-to-VCF pipeline |

### Assembly Tools

| Program | Description |
|---------|-------------|
| `seqassemble` | Assemble sequences from founder FASTA and pedigree segments |
| `vcfassemble` | Assemble VCF from founder VCF and pedigree segments |

### Analysis Tools

| Program | Description |
|---------|-------------|
| `kinship` | Calculate kinship, inbreeding, and relatedness coefficients |
| `het` | Heterozygosity analysis for individuals and loci |
| `gsum` | Summary statistics for genotype data |
| `hwe-dis` | Hardy-Weinberg equilibrium disequilibrium test |
| `sample` | Subsample individuals and markers from VCF files |

---

## Detailed Usage

### coalsim - Coalescent Simulator

Simulates gene genealogies and sequences under the standard coalescent model.

```
coalsim [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-c NUM` | Sample size (haploid chromosomes) | **Required** |
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
| `-a MODE` | MRCA output (r=regions, i=intervals, s=summary) | - |
| `-g MODE` | Gene tree output (s=screen, f=file) | - |
| `-v` | Verbose output | off |
| `-q` | Quiet mode | off |

**Examples:**

```bash
# Basic simulation
coalsim -c 100 -N 10000 -r 1.0 -m 0.5 -V variants.vcf

# HKY model with custom parameters
coalsim -c 50 -M HKY -k 4.0 -p 0.3,0.2,0.2,0.3 -o sequences.fa

# Sparse simulation (target regions only)
coalsim -c 100 -T 1000,100,0,1000 -V sparse.vcf
```

---

### coalsim_msc - Multispecies Coalescent with Migration

Simulates gene genealogies under the multispecies coalescent with optional migration.

```
coalsim_msc [options] <control_file>
```

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Verbose output |
| `-h, --help` | Show help |

**Control File Format:**

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

**Species Tree Format:**

Extended Newick with theta (4Nu) annotations:
```
(children)NAME:branch_length#theta
```

- `NAME`: Optional name for internal nodes (required for ancestral migration)
- `branch_length`: Divergence time in generations
- `theta`: Population parameter 4Nu

**Migration Format:**

```
source->dest:M
```
Where M = 4Nm (population-scaled migration rate). Migration only occurs between contemporary populations.

**Examples:**

Two species with symmetric migration:
```
species_tree: (A:1000#0.01,B:1000#0.01)ROOT#0.01;
samples: A:20,B:20
migration: A->B:1.0,B->A:1.0
mutation_rate: 0.5
output_vcf: true
vcf_file: two_species.vcf
seed: 12345
```

Three species with ancestral migration:
```
species_tree: ((A:300#100,B:300#100)AB:200#100,C:500#100)ROOT#100;
samples: A:10,B:10,C:10
migration: A->B:0.5,B->A:0.5,AB->C:0.1,C->AB:0.1
mutation_rate: 0.5
output_vcf: true
vcf_file: three_species.vcf
```

---

### pedsim - Backwards Pedigree Simulator

Generates random pedigrees backwards in time from a sample in a Wright-Fisher population.

```
pedsim [options]
```

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

**Output Format:**
```
individualID fatherID motherID
```
Founders have fatherID = motherID = 0.

**Example:**
```bash
# Generate pedigree with visualization
pedsim -n 20 -N 500 -k 10 -o pedigree.txt -d pedigree.dot -S

# Convert DOT to image
dot -Tpng pedigree.dot -o pedigree.png
```

---

### pedtrans - Pedigree Chromosome Transmission

Simulates chromosome transmission through a pedigree with recombination.

```
pedtrans [options] <pedigree_file>
```

| Option | Description | Default |
|--------|-------------|---------|
| `-c NUM` | Number of chromosomes | 1 |
| `-r RATE` | Recombination rate (crossovers/meiosis) | 1.0 |
| `-s SEED` | Random seed | time-based |
| `-o FILE` | Output file | stdout |
| `-v` | Verbose output | off |

**Example:**
```bash
pedtrans -c 22 -r 1.5 pedigree.txt -o segments.txt
```

---

### pedsim_multipop - Multi-Population Pedigree Simulator

Extends pedsim to multiple populations with migration.

```
pedsim_multipop [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output pedigree file | **Required** |
| `-J NUM` | Number of populations | 1 |
| `-N LIST` | Population sizes (comma-separated) | 1000 |
| `-n LIST` | Sample sizes (comma-separated) | 10 |
| `-k NUM` | Generations | 5 |
| `-m STR` | Migration matrix ("m00,m01;m10,m11") | - |
| `-M NUM` | Symmetric island model rate | - |
| `-s SEED` | Random seed | time-based |
| `-d` | Output DOT format | off |
| `-v` | Verbose output | off |

**Examples:**
```bash
# Two populations with island model
pedsim_multipop -J 2 -N 500,500 -n 10,10 -M 0.01 -k 10 -o pedigree.txt

# Custom asymmetric migration matrix
pedsim_multipop -J 2 -N 500,500 -n 10,10 -m "0,0.02;0.01,0" -k 10 -o pedigree.txt
```

---

### pedsim_seq - Pedigree Sequence Pipeline

Combines pedsim, coalsim, pedtrans, and seqassemble into one tool.

```
pedsim_seq [options]
```

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

**Example:**
```bash
pedsim_seq -n 50 -N 5000 -k 10 -r 1.0 -m 0.5 -f vcf -o output.vcf
```

---

### pedsim_vcf - Pedigree VCF Pipeline

Streamlined pipeline for VCF output.

```
pedsim_vcf [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-g NUM` | Generations | 5 |
| `-n NUM` | Sample size | 10 |
| `-P NUM` | Pedigree population size | 1000 |
| `-N NUM` | Coalescent effective size | 10000 |
| `-r NUM` | Recombination rate (cM) | 1.0 |
| `-m NUM` | Mutation rate | 0.5 |
| `-s NUM` | Random seed | time-based |
| `-o FILE` | Output VCF file | **Required** |
| `-S` | Samples only | off |
| `-k` | Keep temporary files | off |
| `-v` | Verbose output | off |

**Example:**
```bash
pedsim_vcf -g 10 -n 100 -N 10000 -r 1.5 -m 0.5 -o population.vcf
```

---

### pedsim_vcf_multipop - Multi-Population VCF Pipeline

Complete multi-population simulation to VCF.

```
pedsim_vcf_multipop [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output VCF file | **Required** |
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

**Example:**
```bash
pedsim_vcf_multipop -J 3 -N 1000,500,500 -n 20,10,10 \
    -M 0.01 -g 10 -c 22 -o multipop.vcf
```

---

### seqassemble - Sequence Assembly

Assembles sample sequences from pedigree segment data and founder sequences.

```
seqassemble [options] <segments_file> <founders_fasta>
```

| Option | Description | Default |
|--------|-------------|---------|
| `-o FILE` | Output file | stdout |
| `-f FORMAT` | Output format (fasta/vcf) | fasta |
| `-s` | Samples only | off |
| `-m` | Memory-efficient streaming mode | off |
| `-v` | Verbose output | off |

**Example:**
```bash
seqassemble segments.txt founders.fa -f vcf -o assembled.vcf
```

---

### vcfassemble - VCF Assembly

Assembles sample VCF from founder VCF and pedigree segments.

```
vcfassemble -v <founder.vcf> -p <segments.txt> -o <output.vcf> [-s]
```

| Option | Description |
|--------|-------------|
| `-v FILE` | Founder VCF file |
| `-p FILE` | Pedtrans segment file |
| `-o FILE` | Output VCF file |
| `-s` | Samples only |

---

### kinship - Relatedness Calculator

Calculates kinship, inbreeding, and relatedness coefficients from pedigrees.

```
kinship [options] <pedigree_file>
```

| Option | Description |
|--------|-------------|
| `-i` | Calculate inbreeding coefficients |
| `-k` | Calculate pairwise kinship coefficients |
| `-r` | Calculate pairwise relatedness coefficients |
| `-v` | Verbose output |

**Example:**
```bash
kinship -k -i pedigree.txt > coefficients.txt
```

---

### het - Heterozygosity Analysis

Calculates average heterozygosity for individuals and loci.

```
het [options] <genotype_file>
```

| Option | Description |
|--------|-------------|
| `-i` | Per-individual heterozygosity |
| `-l` | Per-locus heterozygosity |

Output includes 95% confidence intervals. Missing data marked as "M".

---

### gsum - Genotype Summary

Provides summary statistics for genotype data files.

```
gsum [options] <genotype_file>
```

| Option | Description |
|--------|-------------|
| `-i` | Summarize individuals |
| `-p` | Summarize populations |
| `-l` | Summarize loci |

---

### hwe-dis - Hardy-Weinberg Equilibrium Test

Tests for departure from Hardy-Weinberg equilibrium using additive disequilibrium.

```
hwe-dis <genotype_file>
```

Output shows disequilibrium coefficient (D) with 95% CI. Significant departures marked with `*`.

---

### sample - VCF Subsampling

Subsample individuals and/or markers from VCF files.

**Single-file mode:**
```
sample [options] <input.vcf>
```

**Multi-file mode (marker intersection):**
```
sample -O <suffix> -n <num> -m <method> <input1.vcf> <input2.vcf>
```

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

**Examples:**
```bash
# Subsample individuals
sample -I keep_ids.txt input.vcf -o subset.vcf

# Random 1000 markers
sample -n 1000 -m r input.vcf -o subset.vcf

# Uniformly spaced markers
sample -n 1000 -m u input.vcf -o subset.vcf

# Extract region
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
- Founders have `fatherID = motherID = 0`
- IDs are positive integers

**Example:**
```
1 0 0
2 0 0
3 1 2
4 1 2
5 3 4
```

### VCF Format

Standard Variant Call Format. Used as output by coalsim, pipeline programs, and assembly tools.

### FASTA Format

Standard biological sequence format.
```
>sequence_id
ACGTACGTACGT...
```

### Control File Format (coalsim_msc)

Key-value pairs with `#` comments:
```
# Comment
species_tree: ((A:1000#0.01,B:1000#0.01)AB:500#0.01,C:1500#0.01)ROOT#0.01;
samples: A:10,B:10,C:10
migration: A->B:0.5,B->A:0.5
```

---

## Workflows

### Manual Pipeline: Pedigree to Sequences

```bash
# 1. Generate pedigree
pedsim -n 20 -N 1000 -k 5 -o pedigree.txt

# 2. Count founders
FOUNDERS=$(grep " 0 0$" pedigree.txt | wc -l)

# 3. Generate founder sequences
coalsim -c $((2 * FOUNDERS)) -N 10000 -r 1.0 -m 0.5 -o founders.fa

# 4. Simulate transmission
pedtrans -r 1.0 pedigree.txt -o segments.txt

# 5. Assemble sequences
seqassemble segments.txt founders.fa -o samples.fa
```

### One-Command Equivalent

```bash
pedsim_vcf -n 20 -P 1000 -N 10000 -g 5 -r 1.0 -m 0.5 -o samples.vcf
```

### Multi-Population Study

```bash
pedsim_vcf_multipop \
    -J 3 \
    -N 1000,1000,1000 \
    -n 50,50,50 \
    -M 0.01 \
    -g 20 \
    -c 22 \
    -E 10000 \
    -o three_pops.vcf
```

### Species Divergence with Gene Flow

```bash
cat > primates.ctl << 'EOF'
species_tree: ((Human:6000000#0.0001,Chimp:6000000#0.0001)HC:4000000#0.0001,Gorilla:10000000#0.0001)ROOT#0.0001;
samples: Human:50,Chimp:50,Gorilla:20
recombination_rate: 1.0
mutation_rate: 0.5
migration: Human->Chimp:0.1,Chimp->Human:0.1
output_vcf: true
vcf_file: primates.vcf
EOF

coalsim_msc primates.ctl
```

---

## Theoretical Background

### Coalescent Model

The coalescent describes genealogical history backwards in time:

- **θ = 4Nu**: Population-scaled mutation rate
- **ρ = 4Nr**: Population-scaled recombination rate
- **Coalescence rate**: n(n-1)/(4N) for n lineages

### Wright-Fisher Model

Discrete-generation model for pedigree simulation:

- Fixed population size N (N/2 males, N/2 females)
- Random mating each generation
- Non-overlapping generations

### Multispecies Coalescent

Extension for multiple species/populations:

- Lineages coalesce only within the same population
- Populations merge at divergence times (backwards)
- Migration enables lineage exchange between contemporary populations

### Migration Parameters

- **M = 4Nm**: Population-scaled migration rate
- M > 1 means >1 migrant per generation on average
- Higher M → more gene flow → less differentiation

---

## Troubleshooting

### Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Population size must be even" | pedsim requires even N | Use even number for `-N` |
| "Species 'X' not found" | Name mismatch in control file | Check spelling matches tree |
| "Migration must come after species_tree" | Wrong order in control file | Put `species_tree` before `migration` |

### Performance Issues

**Simulation runs very long:**
- Small theta causes slow coalescence → increase theta
- Many samples with recombination → reduce samples or recombination

**Memory errors:**
- Use streaming mode (`-m`) in seqassemble
- Reduce chromosomes or sample size

### Getting Help

```bash
# Any program
<program> -h
<program> --help

# Verbose debugging
<program> -v ...
```

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

Report bugs at: https://github.com/ugnix/ugnix/issues

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Citation

If you use uGnix in published research, please cite:

```
[Citation to be added]
```

---

## Acknowledgments

uGnix builds on foundational work in coalescent theory and population genetics simulation.

**Key References:**

- Kingman JFC (1982) The coalescent. *Stochastic Processes and their Applications*
- Hudson RR (1983) Properties of a neutral allele model with intragenic recombination. *Theoretical Population Biology*
- Rannala B, Yang Z (2003) Bayes estimation of species divergence times. *Genetics*
- Hey J, Nielsen R (2004) Multilocus methods for estimating population sizes, migration rates and divergence time. *Genetics*
