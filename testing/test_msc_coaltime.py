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

def run_msc_simulation_single_species(theta, seed, n_samples=2):
    """
    Run MSC simulation with n samples from one species.
    With very large divergence time, this is effectively single-population coalescent.
    Returns the coalescence time (TMRCA).
    """
    control = f"""
species_tree: A:0#{theta};
samples: A:{n_samples}
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


def run_msc_three_species(theta, tau1, tau2, seed, samples_A=1, samples_B=1, samples_C=1):
    """
    Run MSC with three species: ((A,B):τ1, C):τ2
    τ1 = divergence time of A and B
    τ2 = divergence time of AB ancestor and C (must be > τ1)
    Returns final time (TMRCA).
    """
    control = f"""
species_tree: ((A:{tau1}#{theta},B:{tau1}#{theta}):{tau2-tau1}#{theta},C:{tau2}#{theta})#{theta};
samples: A:{samples_A},B:{samples_B},C:{samples_C}
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


def test_multi_sample_tmrca(theta, n_samples, n_sims):
    """
    Test 5: Multi-sample TMRCA.

    For n samples with coalescence rate n(n-1)/(2θ), the expected TMRCA is:
    E[TMRCA] = 2θ * (1 - 1/n)

    This is the sum of expected waiting times at each stage:
    E[T_k] = 2θ / (k(k-1)) for k = n, n-1, ..., 2
    Sum = 2θ * Σ(1/(k(k-1))) = 2θ * (1 - 1/n)
    """
    print("\n" + "=" * 60)
    print(f"TEST 5: Multi-sample TMRCA (n={n_samples})")
    print("=" * 60)
    print(f"  θ = {theta}")
    print(f"  n = {n_samples} samples")
    theoretical_mean = 2 * theta * (1 - 1/n_samples)
    print(f"  E[TMRCA] = 2θ(1 - 1/n) = {theoretical_mean:.2f}")
    print(f"  Running {n_sims} simulations...")

    times = []
    for i in range(n_sims):
        seed = 50000 + i
        t = run_msc_simulation_single_species(theta, seed, n_samples)
        if t is not None:
            times.append(t)
        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    times = np.array(times)
    empirical_mean = np.mean(times)

    print(f"\nResults ({len(times)} simulations):")
    print(f"  Empirical mean:    {empirical_mean:.2f}")
    print(f"  Theoretical mean:  {theoretical_mean:.2f}")
    print(f"  Ratio (emp/theo):  {empirical_mean/theoretical_mean:.4f}")

    # Use t-test to check if empirical mean matches theoretical
    # H0: true mean = theoretical_mean
    t_stat, p_value = stats.ttest_1samp(times, theoretical_mean)
    print(f"\nt-test (H0: mean = {theoretical_mean:.2f}):")
    print(f"  t-statistic: {t_stat:.4f}")
    print(f"  p-value: {p_value:.4f}")

    passed = p_value > 0.05
    if passed:
        print(f"  Result: PASS (p > 0.05)")
    else:
        print(f"  Result: FAIL (p <= 0.05)")

    return times, passed


def test_theta_scaling(theta1, theta2, n_sims):
    """
    Test 6: Theta scaling.

    Coalescence times should scale linearly with θ.
    If we double θ, E[T] should double.
    """
    print("\n" + "=" * 60)
    print("TEST 6: Theta scaling")
    print("=" * 60)
    print(f"  θ1 = {theta1}")
    print(f"  θ2 = {theta2}")
    ratio = theta2 / theta1
    print(f"  θ2/θ1 = {ratio:.2f}")
    print(f"  Prediction: E[T2]/E[T1] ≈ {ratio:.2f}")
    print(f"  Running {n_sims} simulations for each θ...")

    times1 = []
    times2 = []
    for i in range(n_sims):
        seed1 = 60000 + i
        seed2 = 60000 + i  # Same seeds for paired comparison
        t1 = run_msc_simulation_single_species(theta1, seed1)
        t2 = run_msc_simulation_single_species(theta2, seed2)
        if t1 is not None:
            times1.append(t1)
        if t2 is not None:
            times2.append(t2)
        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    times1 = np.array(times1)
    times2 = np.array(times2)
    mean1 = np.mean(times1)
    mean2 = np.mean(times2)
    empirical_ratio = mean2 / mean1

    print(f"\nResults:")
    print(f"  E[T1] = {mean1:.2f} (θ1 = {theta1})")
    print(f"  E[T2] = {mean2:.2f} (θ2 = {theta2})")
    print(f"  Empirical ratio E[T2]/E[T1] = {empirical_ratio:.4f}")
    print(f"  Theoretical ratio = {ratio:.4f}")

    # Test if the ratio is approximately correct
    # We normalize times by theta and check if distributions match
    normalized1 = times1 / theta1
    normalized2 = times2 / theta2

    # Both should be Exp(1) after normalization
    ks_stat, p_value = stats.ks_2samp(normalized1, normalized2)
    print(f"\nK-S test (normalized times should have same distribution):")
    print(f"  K-S statistic: {ks_stat:.4f}")
    print(f"  p-value: {p_value:.4f}")

    passed = p_value > 0.05
    if passed:
        print(f"  Result: PASS (p > 0.05)")
    else:
        print(f"  Result: FAIL (p <= 0.05)")

    return passed


def test_three_species_divergence(theta, tau1, tau2, n_sims):
    """
    Test 7: Three-species nested divergence.

    For species tree ((A,B):τ1, C):τ2 with 1 sample each:
    - TMRCA must be >= τ2 (all three must coalesce in root population)
    - The gene tree could match species tree or be discordant

    This tests that nested divergence events are handled correctly.
    """
    print("\n" + "=" * 60)
    print("TEST 7: Three-species nested divergence")
    print("=" * 60)
    print(f"  θ = {theta}")
    print(f"  τ1 (A-B divergence) = {tau1}")
    print(f"  τ2 (AB-C divergence) = {tau2}")
    print(f"  Prediction: All TMRCA >= {tau2}")
    print(f"  Running {n_sims} simulations...")

    times = []
    violations = 0

    for i in range(n_sims):
        seed = 70000 + i
        t = run_msc_three_species(theta, tau1, tau2, seed)
        if t is not None:
            times.append(t)
            if t < tau2 - 0.01:
                violations += 1
        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    times = np.array(times)
    print(f"\nResults ({len(times)} simulations):")
    print(f"  Min TMRCA: {np.min(times):.2f}")
    print(f"  Max TMRCA: {np.max(times):.2f}")
    print(f"  Mean TMRCA: {np.mean(times):.2f}")
    print(f"  Root divergence τ2: {tau2:.2f}")
    print(f"  Violations (TMRCA < τ2): {violations}")

    passed = violations == 0
    if passed:
        print(f"  Result: PASS (all TMRCA >= τ2)")
    else:
        print(f"  Result: FAIL ({violations} violations)")

    return times, passed


def test_ancestral_waiting_time(theta_tip, theta_anc, tau, n_sims):
    """
    Test 8: Ancestral population waiting time.

    For 2 samples from different species (1 from A, 1 from B) with divergence at τ:
    - After divergence, both lineages enter ancestral population
    - Time to coalescence in ancestral pop follows Exp(1/θ_anc)
    - So TMRCA - τ ~ Exp(1/θ_anc) with E[TMRCA - τ] = θ_anc
    """
    print("\n" + "=" * 60)
    print("TEST 8: Ancestral population waiting time")
    print("=" * 60)
    print(f"  θ_tip = {theta_tip}")
    print(f"  θ_anc = {theta_anc}")
    print(f"  τ = {tau}")
    print(f"  E[TMRCA - τ] = θ_anc = {theta_anc}")
    print(f"  Running {n_sims} simulations...")

    waiting_times = []
    for i in range(n_sims):
        seed = 80000 + i
        # Use 1 sample from each species - they can only coalesce after τ
        control = f"""
species_tree: (A:{tau}#{theta_tip},B:{tau}#{theta_tip})#{theta_anc};
samples: A:1,B:1
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
                tmrca = float(match.group(1))
                waiting_times.append(tmrca - tau)
        finally:
            os.unlink(ctl_file)

        if (i + 1) % 100 == 0:
            print(f"    Completed {i + 1}/{n_sims}")

    waiting_times = np.array(waiting_times)
    print(f"\nResults ({len(waiting_times)} simulations):")
    print(f"  Empirical mean (TMRCA - τ):  {np.mean(waiting_times):.2f}")
    print(f"  Theoretical mean (θ_anc):    {theta_anc:.2f}")
    print(f"  Empirical std:               {np.std(waiting_times):.2f}")
    print(f"  Theoretical std:             {theta_anc:.2f}")

    # K-S test against exponential with mean θ_anc
    ks_stat, p_value = stats.kstest(waiting_times, 'expon', args=(0, theta_anc))
    print(f"\nKolmogorov-Smirnov test vs Exp(1/θ_anc):")
    print(f"  K-S statistic: {ks_stat:.4f}")
    print(f"  p-value: {p_value:.4f}")

    passed = p_value > 0.05
    if passed:
        print(f"  Result: PASS (p > 0.05)")
    else:
        print(f"  Result: FAIL (p <= 0.05)")

    return waiting_times, passed


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

    # Test 5: Multi-sample TMRCA
    theta5 = 50.0
    n_samples5 = 10
    _, passed5 = test_multi_sample_tmrca(theta5, n_samples5, n_sims)
    results.append(("Multi-sample TMRCA", passed5))

    # Test 6: Theta scaling
    theta6a = 50.0
    theta6b = 150.0  # 3x theta
    passed6 = test_theta_scaling(theta6a, theta6b, n_sims)
    results.append(("Theta scaling", passed6))

    # Test 7: Three-species nested divergence
    theta7 = 100.0
    tau7a = 200.0   # A-B divergence
    tau7b = 500.0   # AB-C divergence
    _, passed7 = test_three_species_divergence(theta7, tau7a, tau7b, n_sims)
    results.append(("Three-species divergence", passed7))

    # Test 8: Ancestral waiting time
    theta8_tip = 100.0  # tip population theta (irrelevant for 1 sample each)
    theta8_anc = 75.0   # ancestral population theta
    tau8 = 300.0
    _, passed8 = test_ancestral_waiting_time(theta8_tip, theta8_anc, tau8, n_sims)
    results.append(("Ancestral waiting time", passed8))

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
