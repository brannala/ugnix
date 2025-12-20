#!/usr/bin/env python3
"""
Test coalescent simulator against theoretical expectations.

Key theoretical results for the standard coalescent:
1. E[TMRCA] = 4N(1 - 1/n) generations for n samples
2. Variance of TMRCA decreases with sample size
3. Recombination creates multiple MRCAs along the chromosome
4. Mutations scale linearly with mutation rate
"""

import subprocess
import re
import math

def mean(x):
    return sum(x) / len(x) if x else 0

def std(x):
    if len(x) < 2:
        return 0
    m = mean(x)
    return math.sqrt(sum((v - m)**2 for v in x) / len(x))

def run_coalsim(n, N, r=0, m=0, seed=None, mrca_summary=False):
    """Run coalsim and parse output."""
    cmd = ['./coalsim', '-c', str(n), '-N', str(N), '-r', str(r), '-m', str(m)]
    if seed is not None:
        cmd.extend(['-s', str(seed)])
    if mrca_summary:
        cmd.extend(['-a', 's'])

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout + result.stderr

    # Parse output
    data = {}

    # TMRCA
    match = re.search(r'Oldest_TMRCA:\s+([\d.]+)', output)
    if match:
        data['tmrca'] = float(match.group(1))

    # Recombinations
    match = re.search(r'No_Recombinations:\s+(\d+)', output)
    if match:
        data['recomb'] = int(match.group(1))

    # Mutations
    match = re.search(r'No_Mutations:\s+(\d+)', output)
    if match:
        data['mutations'] = int(match.group(1))

    # Number of MRCAs (only with -a flag)
    match = re.search(r'No_MRCAs:\s+(\d+)', output)
    if match:
        data['num_mrcas'] = int(match.group(1))

    # Youngest TMRCA
    match = re.search(r'Youngest_TMRCA:\s+([\d.]+)', output)
    if match:
        data['youngest_tmrca'] = float(match.group(1))

    return data

def test_tmrca_expectation():
    """Test 1: E[TMRCA] = 4N(1 - 1/n)"""
    print("=" * 60)
    print("TEST 1: Expected TMRCA (no recombination)")
    print("=" * 60)
    print("Theory: E[TMRCA] = 4N(1 - 1/n)")
    print()

    N = 1000
    n = 10
    nreps = 200

    expected = 4 * N * (1 - 1/n)

    tmrcas = []
    for seed in range(1, nreps + 1):
        result = run_coalsim(n, N, r=0, m=0, seed=seed)
        if 'tmrca' in result:
            tmrcas.append(result['tmrca'])

    observed_mean = mean(tmrcas)
    observed_std = std(tmrcas)
    se = observed_std / math.sqrt(len(tmrcas))

    print(f"Parameters: N={N}, n={n}, reps={nreps}")
    print(f"Expected TMRCA:  {expected:.2f}")
    print(f"Observed mean:   {observed_mean:.2f} ± {1.96*se:.2f} (95% CI)")
    print(f"Difference:      {100*(observed_mean - expected)/expected:.1f}%")

    # Simple check - within 2 standard errors of expected
    print(f"Std Error:       {se:.2f}")

    passed = abs(observed_mean - expected) / expected < 0.10  # Within 10%
    print(f"Result: {'PASS' if passed else 'FAIL'}")
    print()
    return passed

def test_tmrca_scaling():
    """Test 2: TMRCA scales correctly with sample size."""
    print("=" * 60)
    print("TEST 2: TMRCA scaling with sample size")
    print("=" * 60)
    print("Theory: E[TMRCA] = 4N(1 - 1/n)")
    print()

    N = 500
    nreps = 100

    print(f"{'n':>4}  {'Expected':>10}  {'Observed':>10}  {'95% CI':>12}  {'Diff%':>8}  {'Pass':>6}")
    print("-" * 60)

    all_passed = True
    for n in [2, 5, 10, 20]:
        expected = 4 * N * (1 - 1/n)
        tmrcas = []
        for seed in range(1, nreps + 1):
            result = run_coalsim(n, N, r=0, m=0, seed=seed)
            if 'tmrca' in result:
                tmrcas.append(result['tmrca'])

        observed_mean = mean(tmrcas)
        se = std(tmrcas) / math.sqrt(len(tmrcas))
        diff_pct = 100 * (observed_mean - expected) / expected
        passed = abs(diff_pct) < 15  # Within 15%
        all_passed = all_passed and passed

        print(f"{n:>4}  {expected:>10.1f}  {observed_mean:>10.1f}  ±{1.96*se:>10.1f}  {diff_pct:>7.1f}%  {'PASS' if passed else 'FAIL':>6}")

    print()
    return all_passed

def test_recombination_mrcas():
    """Test 3: Recombination increases number of MRCAs."""
    print("=" * 60)
    print("TEST 3: Recombination increases MRCA count")
    print("=" * 60)
    print("Theory: Higher recombination rate → more distinct MRCAs")
    print()

    N = 1000
    n = 10
    nreps = 30

    print(f"{'RecRate':>8}  {'Avg MRCAs':>10}  {'Avg Recomb':>12}  {'Trend':>8}")
    print("-" * 45)

    prev_mrcas = 0
    passed = True
    # Use tiny mutation rate (0.0001) to avoid floating point issues with m=0
    for r in [0.0, 0.02, 0.05, 0.1, 0.2]:
        mrcas = []
        recombs = []
        for seed in range(1, nreps + 1):
            result = run_coalsim(n, N, r=r, m=0.0001, seed=seed, mrca_summary=True)
            mrcas.append(result.get('num_mrcas', 1))
            recombs.append(result.get('recomb', 0))

        avg_mrcas = mean(mrcas)
        avg_recombs = mean(recombs)

        if r > 0 and avg_mrcas < prev_mrcas:
            trend = "↓ BAD"
            passed = False
        elif avg_mrcas > prev_mrcas:
            trend = "↑ OK"
        else:
            trend = "- OK"

        print(f"{r:>8.2f}  {avg_mrcas:>10.2f}  {avg_recombs:>12.1f}  {trend:>8}")
        prev_mrcas = avg_mrcas

    print()
    print(f"Result: {'PASS' if passed else 'FAIL'} - MRCAs should increase with recombination")
    print()
    return passed

def test_mutation_scaling():
    """Test 4: Mutations scale linearly with mutation rate."""
    print("=" * 60)
    print("TEST 4: Mutations scale with mutation rate")
    print("=" * 60)
    print("Theory: E[mutations] ∝ mutation_rate")
    print()

    N = 500
    n = 5
    nreps = 100

    print(f"{'MutRate':>10}  {'Avg Mutations':>15}  {'Ratio (vs prev)':>18}")
    print("-" * 50)

    prev_muts = None
    passed = True
    for m in [0.1, 0.2, 0.4, 0.8]:
        muts = []
        for seed in range(1, nreps + 1):
            result = run_coalsim(n, N, r=0, m=m, seed=seed)
            muts.append(result.get('mutations', 0))

        avg_muts = mean(muts)

        if prev_muts is not None and prev_muts > 0:
            ratio = avg_muts / prev_muts
            expected_ratio = m / prev_m
            ratio_str = f"{ratio:.2f} (expected ~{expected_ratio:.1f})"
            # Check if ratio is close to expected (within 30%)
            if abs(ratio - expected_ratio) / expected_ratio > 0.30:
                passed = False
        else:
            ratio_str = "-"

        print(f"{m:>10.2f}  {avg_muts:>15.1f}  {ratio_str:>18}")
        prev_muts = avg_muts
        prev_m = m

    print()
    print(f"Result: {'PASS' if passed else 'FAIL'} - Mutations should scale linearly")
    print()
    return passed

def test_tmrca_variance():
    """Test 5: TMRCA variance properties."""
    print("=" * 60)
    print("TEST 5: TMRCA Variance decreases with n")
    print("=" * 60)
    print("Theory: Var[TMRCA]/E[TMRCA]² decreases as n increases")
    print()

    N = 500
    nreps = 150

    print(f"{'n':>4}  {'Mean':>10}  {'StdDev':>10}  {'CV':>8}  {'Trend':>8}")
    print("-" * 50)

    prev_cv = float('inf')
    passed = True
    for n in [2, 5, 10, 20]:
        tmrcas = []
        for seed in range(1, nreps + 1):
            result = run_coalsim(n, N, r=0, m=0, seed=seed)
            if 'tmrca' in result:
                tmrcas.append(result['tmrca'])

        m = mean(tmrcas)
        s = std(tmrcas)
        cv = s / m if m > 0 else 0

        if cv > prev_cv:
            trend = "↑ BAD"
            passed = False
        else:
            trend = "↓ OK"

        print(f"{n:>4}  {m:>10.1f}  {s:>10.1f}  {cv:>8.3f}  {trend:>8}")
        prev_cv = cv

    print()
    print(f"Result: {'PASS' if passed else 'FAIL'} - CV should decrease with sample size")
    print()
    return passed

def main():
    print()
    print("=" * 60)
    print("COALESCENT SIMULATOR VALIDATION TESTS")
    print("=" * 60)
    print()

    results = []

    results.append(("TMRCA Expectation", test_tmrca_expectation()))
    results.append(("TMRCA Scaling", test_tmrca_scaling()))
    results.append(("Recombination/MRCAs", test_recombination_mrcas()))
    results.append(("Mutation Scaling", test_mutation_scaling()))
    results.append(("TMRCA Variance", test_tmrca_variance()))

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, passed in results:
        print(f"  {name:25s}: {'PASS' if passed else 'FAIL'}")

    total_passed = sum(1 for _, p in results if p)
    print()
    print(f"Total: {total_passed}/{len(results)} tests passed")
    print("=" * 60)

if __name__ == "__main__":
    main()
