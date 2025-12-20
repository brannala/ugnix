#!/usr/bin/env python3
"""
Test that coalescent gene tree topologies occur with correct frequencies.

For n=4 samples under the Kingman coalescent:
- 15 labeled rooted topologies total
- 3 "balanced" topologies: ((a,b),(c,d)) - each has P = 1/9
- 12 "caterpillar" topologies: (((a,b),c),d) - each has P = 1/18

Overall: P(balanced) = 1/3, P(caterpillar) = 2/3
"""

import subprocess
import re
from collections import Counter
import math

def extract_topology(newick):
    """
    Extract topology from Newick string (remove branch lengths).
    Returns canonical form for comparison.
    """
    # Remove branch lengths
    topo = re.sub(r':\d+\.?\d*', '', newick)
    # Remove semicolon
    topo = topo.rstrip(';').strip()
    return topo

def canonicalize_topology(topo):
    """
    Convert topology to canonical form by sorting subtrees.
    This allows comparing topologies regardless of left/right ordering.
    """
    def parse_and_sort(s):
        s = s.strip()
        if not s.startswith('('):
            return s  # Leaf node

        # Find matching parentheses and split
        level = 0
        parts = []
        start = 1
        for i, c in enumerate(s):
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1
            elif c == ',' and level == 1:
                parts.append(s[start:i])
                start = i + 1
        parts.append(s[start:-1])  # Last part before closing paren

        # Recursively canonicalize and sort
        sorted_parts = sorted(parse_and_sort(p) for p in parts)
        return '(' + ','.join(sorted_parts) + ')'

    return parse_and_sort(topo)

def get_topology_class(topo):
    """
    Determine if topology is 'balanced' or 'caterpillar'.
    Balanced: ((a,b),(c,d)) - both children of root are cherries
    Caterpillar: (((a,b),c),d) - one child of root is a single tip
    """
    # Count depth of each tip
    depths = {}
    current_depth = 0
    current_tip = ''

    for c in topo:
        if c == '(':
            current_depth += 1
        elif c == ')':
            if current_tip:
                depths[current_tip] = current_depth
                current_tip = ''
            current_depth -= 1
        elif c == ',':
            if current_tip:
                depths[current_tip] = current_depth
                current_tip = ''
        elif c.isdigit():
            current_tip += c

    if current_tip:
        depths[current_tip] = current_depth

    # Balanced: all tips at same depth (2 for 4 taxa)
    # Caterpillar: tips at different depths
    depth_values = list(depths.values())
    if len(set(depth_values)) == 1:
        return 'balanced'
    else:
        return 'caterpillar'

def run_simulation(n_samples, seed):
    """Run coalsim and extract tree topology."""
    cmd = ['./coalsim', '-c', str(n_samples), '-N', '1000', '-r', '0', '-m', '0',
           '-s', str(seed), '-g', 's']
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Extract tree from output
    for line in result.stderr.split('\n'):
        if line.startswith('['):
            match = re.match(r'\[[\d.]+,[\d.]+\]\s+(.+)', line)
            if match:
                return match.group(1)
    return None

def main():
    print("=" * 70)
    print("COALESCENT TOPOLOGY FREQUENCY TEST")
    print("=" * 70)
    print()
    print("Theory for n=4 samples:")
    print("  - 3 balanced topologies ((a,b),(c,d)): P = 1/9 each = 0.1111")
    print("  - 12 caterpillar topologies (((a,b),c),d): P = 1/18 each = 0.0556")
    print("  - P(balanced class) = 1/3 = 0.3333")
    print("  - P(caterpillar class) = 2/3 = 0.6667")
    print()

    n_sims = 5000
    n_samples = 4

    print(f"Running {n_sims} simulations with n={n_samples} samples...")
    print()

    topologies = []
    canonical_topos = []
    classes = []

    for seed in range(1, n_sims + 1):
        tree = run_simulation(n_samples, seed)
        if tree:
            topo = extract_topology(tree)
            canon = canonicalize_topology(topo)
            topo_class = get_topology_class(topo)

            topologies.append(topo)
            canonical_topos.append(canon)
            classes.append(topo_class)

        if seed % 1000 == 0:
            print(f"  Completed {seed}/{n_sims} simulations...")

    print()

    # Count topology classes
    class_counts = Counter(classes)
    print("=" * 70)
    print("TOPOLOGY CLASS FREQUENCIES")
    print("=" * 70)
    print()
    print(f"{'Class':<15} {'Observed':<10} {'Expected':<10} {'Obs Freq':<12} {'Exp Freq':<12}")
    print("-" * 60)

    n_total = len(classes)
    balanced_obs = class_counts.get('balanced', 0)
    caterpillar_obs = class_counts.get('caterpillar', 0)

    balanced_exp = n_total / 3
    caterpillar_exp = 2 * n_total / 3

    print(f"{'Balanced':<15} {balanced_obs:<10} {balanced_exp:<10.1f} {balanced_obs/n_total:<12.4f} {1/3:<12.4f}")
    print(f"{'Caterpillar':<15} {caterpillar_obs:<10} {caterpillar_exp:<10.1f} {caterpillar_obs/n_total:<12.4f} {2/3:<12.4f}")

    # Chi-square test for class frequencies
    chi2_class = ((balanced_obs - balanced_exp)**2 / balanced_exp +
                  (caterpillar_obs - caterpillar_exp)**2 / caterpillar_exp)
    # Critical value for df=1, alpha=0.05 is 3.841
    class_pass = chi2_class < 3.841

    print()
    print(f"Chi-square statistic: {chi2_class:.4f} (critical value at p=0.05: 3.841)")
    print(f"Result: {'PASS' if class_pass else 'FAIL'}")
    print()

    # Count individual canonical topologies
    topo_counts = Counter(canonical_topos)

    print("=" * 70)
    print("INDIVIDUAL TOPOLOGY FREQUENCIES")
    print("=" * 70)
    print()

    # Separate balanced and caterpillar topologies
    balanced_topos = {}
    caterpillar_topos = {}

    for topo, count in topo_counts.items():
        if get_topology_class(topo) == 'balanced':
            balanced_topos[topo] = count
        else:
            caterpillar_topos[topo] = count

    print("Balanced topologies (expected freq = 1/9 = 0.1111 each):")
    print("-" * 60)
    exp_balanced = n_total / 9
    chi2_balanced = 0
    for topo, count in sorted(balanced_topos.items()):
        freq = count / n_total
        chi2_balanced += (count - exp_balanced)**2 / exp_balanced
        status = "✓" if abs(freq - 1/9) < 0.02 else "?"
        print(f"  {status} {topo}: {count:5d} ({freq:.4f})")

    print()
    print(f"  Chi-square for balanced uniformity: {chi2_balanced:.4f}")
    # df = 3-1 = 2, critical value at p=0.05 is 5.991
    print(f"  (critical value at p=0.05, df=2: 5.991) {'PASS' if chi2_balanced < 5.991 else 'FAIL'}")

    print()
    print("Caterpillar topologies (expected freq = 1/18 = 0.0556 each):")
    print("-" * 60)
    exp_caterpillar = n_total / 18
    chi2_caterpillar = 0
    for topo, count in sorted(caterpillar_topos.items()):
        freq = count / n_total
        chi2_caterpillar += (count - exp_caterpillar)**2 / exp_caterpillar
        status = "✓" if abs(freq - 1/18) < 0.015 else "?"
        print(f"  {status} {topo}: {count:5d} ({freq:.4f})")

    print()
    print(f"  Chi-square for caterpillar uniformity: {chi2_caterpillar:.4f}")
    # df = 12-1 = 11, critical value at p=0.05 is 19.675
    print(f"  (critical value at p=0.05, df=11: 19.675) {'PASS' if chi2_caterpillar < 19.675 else 'FAIL'}")

    # Overall test
    print()
    print("=" * 70)
    print("OVERALL UNIFORMITY TEST (all 15 topologies)")
    print("=" * 70)

    chi2_all = 0
    for topo, count in topo_counts.items():
        if get_topology_class(topo) == 'balanced':
            expected = n_total / 9
        else:
            expected = n_total / 18
        chi2_all += (count - expected)**2 / expected

    # df = 15-1 = 14, critical value at p=0.05 is 23.685
    print(f"Chi-square statistic: {chi2_all:.4f}")
    print(f"Critical value at p=0.05, df=14: 23.685")
    overall_pass = chi2_all < 23.685
    print(f"Result: {'PASS' if overall_pass else 'FAIL'}")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    all_pass = class_pass and overall_pass
    print(f"  Class frequencies (balanced vs caterpillar): {'PASS' if class_pass else 'FAIL'}")
    print(f"  Overall topology uniformity: {'PASS' if overall_pass else 'FAIL'}")
    print(f"  Final result: {'ALL TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")
    print("=" * 70)

if __name__ == "__main__":
    main()
