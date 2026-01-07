#!/usr/bin/env python3
"""
Test pedtrans by verifying expected kinship coefficients between relatives.

Kinship coefficient = E[(# IBD allele pairs) / 4] = expected proportion of alleles IBD

Expected values (no inbreeding):
- Parent-child: 0.25
- Full siblings: 0.25
- Half-siblings: 0.125
- Grandparent-grandchild: 0.125
- Aunt/Uncle - Niece/Nephew: 0.125
- First cousins: 0.0625
"""

import subprocess
import tempfile
import os
import re
from collections import defaultdict
import statistics

def create_sibling_pedigree():
    """Create a simple pedigree with 2 founders and 2 siblings."""
    return """# Two founders, two siblings
Father 0 0
Mother 0 0
Sib1 Father Mother
Sib2 Father Mother
"""

def create_extended_pedigree():
    """
    Create a 4-generation pedigree for testing various relationships.

    Gen 0 (founders): F1-F8 (8 founders = 4 couples)
    Gen 1: 4 individuals (children of founder couples)
    Gen 2: 4 individuals (2 pairs of siblings from 2 couples)
    Gen 3: 4 individuals (children, including cousins)

    Relationships to test:
    - Siblings (same parents)
    - Aunt/Uncle to Niece/Nephew
    - First cousins
    """
    return """# Extended pedigree for relationship testing
# Generation 0 - Founders (4 couples)
F1 0 0
F2 0 0
F3 0 0
F4 0 0
F5 0 0
F6 0 0
F7 0 0
F8 0 0
# Generation 1 - Children of founders
G1_A F1 F2
G1_B F3 F4
G1_C F5 F6
G1_D F7 F8
# Generation 2 - Two pairs of siblings
G2_A1 G1_A G1_B
G2_A2 G1_A G1_B
G2_B1 G1_C G1_D
G2_B2 G1_C G1_D
# Generation 3 - Children (cousins to each other within groups)
G3_A1 G2_A1 G2_B1
G3_A2 G2_A1 G2_B1
G3_B1 G2_A2 G2_B2
G3_B2 G2_A2 G2_B2
"""

def run_pedtrans(ped_content, rec_rate=1.0, seed=None):
    """Run pedtrans and return the output."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.ped', delete=False) as f:
        f.write(ped_content)
        ped_file = f.name

    try:
        cmd = ['./pedtrans', '-r', str(rec_rate)]
        if seed is not None:
            cmd.extend(['-s', str(seed)])
        cmd.append(ped_file)

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            return None
        return result.stdout
    finally:
        os.unlink(ped_file)

def parse_output(output):
    """
    Parse pedtrans output into a dictionary of individuals -> chromosomes -> segments.

    Returns: {indiv_name: {'Paternal': [(start, end, founder, homolog), ...],
                           'Maternal': [(start, end, founder, homolog), ...]}}
    """
    individuals = {}
    current_indiv = None
    current_chrom = None

    for line in output.split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        if line.startswith('Individual:'):
            current_indiv = line.split(':')[1].strip()
            individuals[current_indiv] = {'Paternal': [], 'Maternal': []}
        elif line.startswith('Paternal:'):
            current_chrom = 'Paternal'
        elif line.startswith('Maternal:'):
            current_chrom = 'Maternal'
        elif current_indiv and current_chrom:
            # Parse segment: "0.000000 0.456789 FounderName:pat"
            parts = line.split()
            if len(parts) >= 3:
                start = float(parts[0])
                end = float(parts[1])
                origin = parts[2]
                # Parse origin "FounderName:pat" or "FounderName:mat"
                if ':' in origin:
                    founder, homolog = origin.rsplit(':', 1)
                    homolog = 0 if homolog == 'pat' else 1
                    individuals[current_indiv][current_chrom].append(
                        (start, end, founder, homolog))

    return individuals

def calculate_ibd_sharing(indiv1_data, indiv2_data):
    """
    Calculate the kinship coefficient between two individuals.

    The kinship coefficient is the expected proportion of alleles shared IBD,
    which equals E[(# of IBD allele pairs) / 4] at each locus, averaged over genome.

    For diploid individuals, there are 4 possible allele pairings to check.
    """
    # Get all breakpoints
    breakpoints = set([0.0, 1.0])
    for chrom in ['Paternal', 'Maternal']:
        for start, end, _, _ in indiv1_data[chrom]:
            breakpoints.add(start)
            breakpoints.add(end)
        for start, end, _, _ in indiv2_data[chrom]:
            breakpoints.add(start)
            breakpoints.add(end)

    breakpoints = sorted(breakpoints)

    total_kinship = 0.0

    for i in range(len(breakpoints) - 1):
        seg_start = breakpoints[i]
        seg_end = breakpoints[i + 1]
        seg_mid = (seg_start + seg_end) / 2
        seg_len = seg_end - seg_start

        # Find which founder chromosome each individual has at this position
        def get_origin_at_pos(indiv_data, chrom, pos):
            for start, end, founder, homolog in indiv_data[chrom]:
                if start <= pos < end:
                    return (founder, homolog)
            return None

        # Get origins for both individuals at this position
        i1_pat = get_origin_at_pos(indiv1_data, 'Paternal', seg_mid)
        i1_mat = get_origin_at_pos(indiv1_data, 'Maternal', seg_mid)
        i2_pat = get_origin_at_pos(indiv2_data, 'Paternal', seg_mid)
        i2_mat = get_origin_at_pos(indiv2_data, 'Maternal', seg_mid)

        # Count IBD pairs out of 4 possible pairings
        # (i1_pat, i2_pat), (i1_pat, i2_mat), (i1_mat, i2_pat), (i1_mat, i2_mat)
        ibd_count = 0
        if i1_pat and i2_pat and i1_pat == i2_pat:
            ibd_count += 1
        if i1_pat and i2_mat and i1_pat == i2_mat:
            ibd_count += 1
        if i1_mat and i2_pat and i1_mat == i2_pat:
            ibd_count += 1
        if i1_mat and i2_mat and i1_mat == i2_mat:
            ibd_count += 1

        # Kinship contribution from this segment = (ibd_count / 4) * seg_len
        total_kinship += (ibd_count / 4.0) * seg_len

    return total_kinship

def run_simulation(ped_content, n_reps=1000, rec_rate=1.0, pairs_to_test=None):
    """
    Run multiple simulations and calculate average IBD sharing.

    pairs_to_test: list of (name1, name2, expected_sharing, description) tuples
    """
    if pairs_to_test is None:
        pairs_to_test = []

    # Store results for each pair
    results = {desc: [] for _, _, _, desc in pairs_to_test}

    for rep in range(n_reps):
        if (rep + 1) % 100 == 0:
            print(f"  Replicate {rep + 1}/{n_reps}...")

        output = run_pedtrans(ped_content, rec_rate=rec_rate, seed=rep + 1)
        if output is None:
            continue

        individuals = parse_output(output)

        for name1, name2, expected, desc in pairs_to_test:
            if name1 in individuals and name2 in individuals:
                sharing = calculate_ibd_sharing(individuals[name1], individuals[name2])
                results[desc].append(sharing)

    return results

def main():
    print("=" * 60)
    print("Testing pedtrans IBD sharing")
    print("=" * 60)

    # Test 1: Sibling sharing (expected 50%)
    print("\n" + "-" * 60)
    print("Test 1: Sibling kinship coefficient (expected: 0.25)")
    print("-" * 60)

    ped = create_sibling_pedigree()
    pairs = [('Sib1', 'Sib2', 0.25, 'Siblings')]  # Kinship coefficient for siblings = 0.25

    results = run_simulation(ped, n_reps=1000, rec_rate=1.0, pairs_to_test=pairs)

    for desc, values in results.items():
        if values:
            mean_sharing = statistics.mean(values)
            std_sharing = statistics.stdev(values) if len(values) > 1 else 0
            print(f"\n{desc}:")
            print(f"  Mean IBD sharing: {mean_sharing:.4f} ({mean_sharing*100:.2f}%)")
            print(f"  Std deviation:    {std_sharing:.4f}")
            print(f"  Expected:         0.2500 (25.00%)")
            print(f"  N replicates:     {len(values)}")

    # Test 2: Extended pedigree with multiple relationship types
    print("\n" + "-" * 60)
    print("Test 2: Extended pedigree - multiple relationships")
    print("-" * 60)

    ped = create_extended_pedigree()
    pairs = [
        ('G2_A1', 'G2_A2', 0.25, 'Siblings (G2)'),
        ('G3_A1', 'G3_A2', 0.25, 'Siblings (G3)'),
        ('G2_A1', 'G3_A1', 0.25, 'Parent-Child'),
        ('G2_A2', 'G3_A1', 0.125, 'Aunt/Uncle-Niece/Nephew'),
    ]

    results = run_simulation(ped, n_reps=1000, rec_rate=1.0, pairs_to_test=pairs)

    print("\nResults:")
    print(f"{'Relationship':<30} {'Mean':>8} {'Std':>8} {'Expected':>10}")
    print("-" * 60)

    for desc, values in results.items():
        if values:
            mean_sharing = statistics.mean(values)
            std_sharing = statistics.stdev(values) if len(values) > 1 else 0
            expected = [e for n1, n2, e, d in pairs if d == desc][0]
            status = "OK" if abs(mean_sharing - expected) < 0.02 else "CHECK"
            print(f"{desc:<30} {mean_sharing:>8.4f} {std_sharing:>8.4f} {expected:>10.4f}  {status}")

    print("\n" + "=" * 60)
    print("Testing complete")
    print("=" * 60)

if __name__ == '__main__':
    main()
