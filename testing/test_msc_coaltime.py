#!/usr/bin/env python3
"""
Test MSC coalescence time distribution against theoretical prediction.

For 2 lineages with coalescence rate 1/θ, conditional on coalescing before
divergence time τ, the density is:

f(t | t < τ) = (1/θ) * exp(-t/θ) / (1 - exp(-τ/θ))  for 0 < t < τ

This is a truncated exponential distribution.

Test approach:
- Use 2 samples from species A only (single species effectively)
- With divergence time τ very large, coalescence always occurs before τ
- Compare empirical distribution to Exp(1/θ)

Then test the conditional case:
- Use 2 samples from A, with moderate τ
- The proportion that coalesce before τ should match 1 - exp(-τ/θ)
- Conditional on pre-τ coalescence, times should follow truncated exponential
"""

import subprocess
import numpy as np
import re
import sys
from scipy import stats
import tempfile
import os

def run_msc_simulation_single_species(theta, seed):
    """
    Run MSC simulation with just 2 samples from one species.
    With very large divergence time, this is effectively single-population coalescent.
    Returns the coalescence time.
    """
    # Use very large divergence time so coalescence always happens before it
    tau = 1e9  # Effectively infinite

    control = f"""
species_tree: A:0#{theta};
samples: A:2
recombination_rate: 0.0
mutation_rate: 0.0
seed: {seed}
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.ctl', delete=False) as f:
        f.write(control)
        ctl_file = f.name

    try:
        result = subprocess.run(
            ['./coalsim_msc', ctl_file],
            capture_output=True, text=True, cwd='/home/bruce/repos/ugnix'
        )
        output = result.stdout + result.stderr

        # Parse the final time from output
        match = re.search(r'Final time:\s+([\d.]+)', output)
        if match:
            return float(match.group(1))
        return None
    finally:
        os.unlink(ctl_file)


def run_msc_two_species(theta_A, theta_anc, tau, seed):
    """
    Run MSC with 2 samples from A, 1 from B.
    Returns (final_time, divergence_time).

    If final_time ≈ tau, the A samples coalesced before divergence.
    If final_time > tau, the A samples coalesced in the ancestral population.
    """
    control = f"""
species_tree: (A:{tau}#{theta_A},B:{tau}#{theta_A})#{theta_anc};
samples: A:2,B:1
recombination_rate: 0.0
mutation_rate: 0.0
seed: {seed}
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.ctl', delete=False) as f:
        f.write(control)
        ctl_file = f.name

    try:
        result = subprocess.run(
            ['./coalsim_msc', ctl_file],
            capture_output=True, text=True, cwd='/home/bruce/repos/ugnix'
        )
        output = result.stdout + result.stderr

        match = re.search(r'Final time:\s+([\d.]+)', output)
        if match:
            return float(match.group(1))
        return None
    finally:
        os.unlink(ctl_file)

def theoretical_cdf(t, theta, tau):
    """CDF of truncated exponential: P(T <= t | T < tau)"""
    if t <= 0:
        return 0.0
    if t >= tau:
        return 1.0
    # P(T <= t | T < tau) = P(T <= t) / P(T < tau) = (1 - exp(-t/θ)) / (1 - exp(-τ/θ))
    return (1 - np.exp(-t/theta)) / (1 - np.exp(-tau/theta))

def theoretical_quantile(p, theta, tau):
    """Quantile function of truncated exponential"""
    # F(t) = (1 - exp(-t/θ)) / (1 - exp(-τ/θ)) = p
    # 1 - exp(-t/θ) = p * (1 - exp(-τ/θ))
    # exp(-t/θ) = 1 - p * (1 - exp(-τ/θ))
    # t = -θ * ln(1 - p * (1 - exp(-τ/θ)))
    return -theta * np.log(1 - p * (1 - np.exp(-tau/theta)))

def test_single_species_exponential(theta, n_sims):
    """
    Test 1: Single species with 2 samples.
    Coalescence time should follow Exp(1/θ) with E[T] = θ.
    """
    print("=" * 60)
    print("TEST 1: Single-species coalescence ~ Exp(1/θ)")
    print("=" * 60)
    print(f"  θ = {theta}")
    print(f"  E[T] = θ = {theta}")
    print(f"  Running {n_sims} simulations...")

    times = []
    for i in range(n_sims):
        seed = 10000 + i
        t = run_msc_simulation_single_species(theta, seed)
        if t is not None:
            times.append(t)
        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    times = np.array(times)
    print(f"\nResults ({len(times)} successful simulations):")
    print(f"  Empirical mean:    {np.mean(times):.2f}")
    print(f"  Theoretical mean:  {theta:.2f}")
    print(f"  Empirical std:     {np.std(times):.2f}")
    print(f"  Theoretical std:   {theta:.2f}")  # std of Exp(1/θ) is also θ

    # K-S test against exponential distribution
    # scipy.stats.expon has scale parameter = mean = θ
    ks_stat, p_value = stats.kstest(times, 'expon', args=(0, theta))
    print(f"\nKolmogorov-Smirnov test vs Exp(1/θ):")
    print(f"  K-S statistic: {ks_stat:.4f}")
    print(f"  p-value: {p_value:.4f}")
    if p_value > 0.05:
        print(f"  Result: PASS (p > 0.05)")
    else:
        print(f"  Result: FAIL (p <= 0.05)")

    return times, p_value > 0.05


def test_pre_divergence_proportion(theta, tau, n_sims):
    """
    Test 2: Proportion of coalescences occurring before divergence time τ.
    Should match P(T < τ) = 1 - exp(-τ/θ).

    We use 2 samples from A and 1 from B. Track when final coalescence occurs.
    If A samples coalesce before τ, they form one lineage, then merge with B at τ,
    then the final coalescence happens in the ancestral population.

    Actually, the FINAL time tells us when the last coalescence happened.
    With θ_anc very small, ancestral coalescence is fast, so:
    - If A coalesces before τ: final time ≈ τ + small
    - If A doesn't coalesce before τ: final time > τ (could be much larger)

    Better approach: Use only species A with 2 samples, single tip tree.
    Compare final time to τ directly.
    """
    print("\n" + "=" * 60)
    print("TEST 2: Pre-divergence coalescence proportion")
    print("=" * 60)
    print(f"  θ = {theta}")
    print(f"  τ = {tau}")
    print(f"  P(T < τ) = 1 - exp(-τ/θ) = {1 - np.exp(-tau/theta):.4f}")
    print(f"  Running {n_sims} simulations...")

    times = []
    for i in range(n_sims):
        seed = 20000 + i
        t = run_msc_simulation_single_species(theta, seed)
        if t is not None:
            times.append(t)
        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    times = np.array(times)
    n_before_tau = np.sum(times < tau)
    empirical_prop = n_before_tau / len(times)
    theoretical_prop = 1 - np.exp(-tau/theta)

    print(f"\nResults ({len(times)} simulations):")
    print(f"  Coalescences before τ: {n_before_tau}/{len(times)}")
    print(f"  Empirical proportion:   {empirical_prop:.4f}")
    print(f"  Theoretical proportion: {theoretical_prop:.4f}")

    # Binomial test
    # H0: true proportion = theoretical_prop
    result = stats.binomtest(n_before_tau, len(times), theoretical_prop)
    print(f"\nBinomial test:")
    print(f"  p-value: {result.pvalue:.4f}")
    if result.pvalue > 0.05:
        print(f"  Result: PASS (p > 0.05)")
    else:
        print(f"  Result: FAIL (p <= 0.05)")

    return empirical_prop, theoretical_prop, result.pvalue > 0.05


def test_truncated_exponential(theta, tau, n_sims):
    """
    Test 3: Conditional distribution of coalescence times given T < τ.
    Should follow truncated exponential.
    """
    print("\n" + "=" * 60)
    print("TEST 3: Truncated exponential distribution (T | T < τ)")
    print("=" * 60)
    print(f"  θ = {theta}")
    print(f"  τ = {tau}")

    # Expected mean of truncated exponential
    exp_term = np.exp(-tau/theta)
    theoretical_mean = theta * (1 - (tau/theta + 1)*exp_term) / (1 - exp_term)
    print(f"  E[T | T < τ] = {theoretical_mean:.2f}")
    print(f"  Running {n_sims} simulations...")

    times = []
    for i in range(n_sims):
        seed = 30000 + i
        t = run_msc_simulation_single_species(theta, seed)
        if t is not None and t < tau:
            times.append(t)
        if (i + 1) % 200 == 0:
            print(f"    Completed {i + 1}/{n_sims} (collected {len(times)} pre-τ)")

    if len(times) < 30:
        print(f"\nToo few samples ({len(times)}). Increase n_sims or τ/θ ratio.")
        return None, False

    times = np.array(times)
    print(f"\nResults ({len(times)} pre-divergence coalescences):")
    print(f"  Empirical mean:    {np.mean(times):.2f}")
    print(f"  Theoretical mean:  {theoretical_mean:.2f}")

    # K-S test against truncated exponential
    sorted_times = np.sort(times)
    n = len(sorted_times)
    empirical_cdf = np.arange(1, n + 1) / n
    theoretical_cdf_values = np.array([theoretical_cdf(t, theta, tau) for t in sorted_times])

    ks_stat = np.max(np.abs(empirical_cdf - theoretical_cdf_values))
    # Approximate p-value using scipy's kstest internals
    # For custom CDF, we use the K-S statistic directly
    critical_value_05 = 1.36 / np.sqrt(n)
    critical_value_01 = 1.63 / np.sqrt(n)

    print(f"\nKolmogorov-Smirnov test vs truncated exponential:")
    print(f"  K-S statistic: {ks_stat:.4f}")
    print(f"  Critical value (α=0.05): {critical_value_05:.4f}")
    print(f"  Critical value (α=0.01): {critical_value_01:.4f}")

    passed = ks_stat < critical_value_05
    if passed:
        print(f"  Result: PASS (K-S stat < critical value)")
    else:
        print(f"  Result: FAIL (K-S stat >= critical value)")

    # Q-Q comparison
    print(f"\nQ-Q comparison (percentiles):")
    print("  %ile\tTheoretical\tEmpirical")
    for p in [10, 25, 50, 75, 90]:
        emp_q = np.percentile(times, p)
        theo_q = theoretical_quantile(p/100, theta, tau)
        print(f"  {p}\t{theo_q:.2f}\t\t{emp_q:.2f}")

    return times, passed


def test_population_isolation(theta, tau, n_sims):
    """
    Test 4: Population isolation - TMRCA must be >= divergence time.

    For a two-species tree (A,B):τ with samples from both species,
    the final coalescence (between A and B lineages) can only occur
    at time >= τ.
    """
    print("\n" + "=" * 60)
    print("TEST 4: Population isolation (TMRCA >= τ)")
    print("=" * 60)
    print(f"  θ = {theta}")
    print(f"  τ = {tau}")
    print(f"  Prediction: All final times must be >= {tau}")
    print(f"  Running {n_sims} simulations...")

    times = []
    violations = 0

    for i in range(n_sims):
        seed = 40000 + i
        t = run_msc_two_species(theta, theta, tau, seed)
        if t is not None:
            times.append(t)
            if t < tau - 0.01:  # Small tolerance for floating point
                violations += 1
        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    times = np.array(times)
    print(f"\nResults ({len(times)} simulations):")
    print(f"  Min TMRCA: {np.min(times):.2f}")
    print(f"  Max TMRCA: {np.max(times):.2f}")
    print(f"  Mean TMRCA: {np.mean(times):.2f}")
    print(f"  Divergence time τ: {tau:.2f}")
    print(f"  Violations (TMRCA < τ): {violations}")

    passed = violations == 0
    if passed:
        print(f"  Result: PASS (all TMRCA >= τ)")
    else:
        print(f"  Result: FAIL ({violations} violations)")

    return times, passed


def main():
    print("MSC Coalescence Time Distribution Tests")
    print("=" * 60)
    print()

    n_sims = 500
    results = []

    # Test 1: Basic exponential distribution
    theta1 = 100.0
    _, passed1 = test_single_species_exponential(theta1, n_sims)
    results.append(("Exp(1/θ) distribution", passed1))

    # Test 2: Pre-divergence proportion
    theta2 = 100.0
    tau2 = 150.0  # τ/θ = 1.5, so P(T<τ) ≈ 0.78
    _, _, passed2 = test_pre_divergence_proportion(theta2, tau2, n_sims)
    results.append(("Pre-divergence proportion", passed2))

    # Test 3: Truncated exponential
    theta3 = 100.0
    tau3 = 200.0  # τ/θ = 2, so P(T<τ) ≈ 0.86
    _, passed3 = test_truncated_exponential(theta3, tau3, n_sims * 2)  # More samples needed
    results.append(("Truncated exponential", passed3))

    # Test 4: Population isolation
    theta4 = 50.0  # Small theta so within-species coalescence is fast
    tau4 = 500.0   # Large divergence time
    _, passed4 = test_population_isolation(theta4, tau4, n_sims)
    results.append(("Population isolation", passed4))

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    all_passed = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("All tests PASSED!")
        return 0
    else:
        print("Some tests FAILED - check implementation.")
        return 1

if __name__ == '__main__':
    main()
