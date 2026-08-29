"""
run_simulation.py

Runs the full simulation-only proof of concept for
"Spectral Silences in Raman Spectroscopy of Tooth Enamel" (Manuscript ID 2979716)
and prints every numeric result reported in the manuscript's Sections 5-6 and
Tables 1-2, using the exact fixed parameters in Table 1 (seed = 42).

Run with:  python run_simulation.py
Outputs numeric results to stdout and writes table2_detection_benchmark.csv
and silence_scores.csv to the working directory.
"""

import os
import sys
import numpy as np
import csv

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from spectral_silence import (
    RANDOM_SEED, SPECTRAL_MIN, SPECTRAL_MAX, STEP, K_VALUES, AMPLITUDES,
    CANDIDATE_WINDOWS, MASK_INTERVALS,
    generate_simulated_spectrum, process_spectrum, monte_carlo_benchmark,
)

def main():
    nu = np.arange(SPECTRAL_MIN, SPECTRAL_MAX + STEP, STEP)
    rng = np.random.default_rng(RANDOM_SEED)

    # 1. Generate the base synthetic spectrum (Figure 2)
    I, baseline, M, O = generate_simulated_spectrum(nu, rng=rng)

    # 2. Silence-score processing across candidate windows and k grid (Figure 3)
    proc = process_spectrum(I, nu)
    print(f"Robust noise sigma (first-difference MAD): {proc['sigma']:.4f}\n")

    print("Silence scores S*(W):")
    mean_by_window_k2 = {}
    with open('silence_scores.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['window_low', 'window_high', 'k', 'S_star', 'Q90'])
        for (W, k), val in proc['results'].items():
            writer.writerow([W[0], W[1], k, round(val['S'], 4), round(val['q90'], 4)])
            if k == 2.0:
                mean_by_window_k2[W] = val['S']
            print(f"  Window {W}, k={k}: S*(W)={val['S']:.3f}, Q90={val['q90']:.3f}")

    best_window = max(mean_by_window_k2, key=mean_by_window_k2.get)
    print(f"\nBest window at k=2.0: {best_window}, S*={mean_by_window_k2[best_window]:.3f}")

    # 3. Monte Carlo spike-in benchmark (Table 2, Figures 4-5)
    print(f"\nRunning Monte Carlo benchmark: {len(AMPLITUDES)} amplitudes x "
          f"{200} trials (plus {200} null trials) ...")
    bench = monte_carlo_benchmark(baseline, M, O, nu, seed=RANDOM_SEED)

    print(f"\nZero-amplitude false-positive rates:")
    print(f"  Matched-filter test:      {bench['mf_fpr0']:.3f}")
    print(f"  Conventional local SNR:   {bench['snr_fpr0']:.3f}")

    print(f"\nPooled AUC:")
    print(f"  Matched-filter test:      {bench['auc_mf']:.3f}")
    print(f"  Conventional local SNR:   {bench['auc_snr']:.3f}")

    print("\nDetection probabilities by amplitude:")
    with open('table2_detection_benchmark.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['injected_amplitude', 'method', 'detection_probability',
                          'ci_95_halfwidth', 'false_positive_rate', 'auc'])
        # amplitude 0 rows
        writer.writerow([0.000, 'Matched-filter test', round(bench['mf_fpr0'], 3),
                          '', round(bench['mf_fpr0'], 3), round(bench['auc_mf'], 3)])
        writer.writerow([0.000, 'Conventional local SNR', round(bench['snr_fpr0'], 3),
                          '', round(bench['snr_fpr0'], 3), round(bench['auc_snr'], 3)])
        for A in AMPLITUDES:
            d = bench['detection_summary'][A]
            print(f"  A={A}: matched-filter={d['mf_detection_prob']:.3f}  "
                  f"local-SNR={d['snr_detection_prob']:.3f}")
            writer.writerow([A, 'Matched-filter test', round(d['mf_detection_prob'], 3),
                              round(d['mf_ci_halfwidth'], 3), '', round(bench['auc_mf'], 3)])
            writer.writerow([A, 'Conventional local SNR', round(d['snr_detection_prob'], 3),
                              round(d['snr_ci_halfwidth'], 3), '', round(bench['auc_snr'], 3)])

    print("\nWrote silence_scores.csv and table2_detection_benchmark.csv")
    print("\nDone. Compare the printed values above against Table 2 and Section 5-6 "
          "of the manuscript.")

if __name__ == "__main__":
    main()
