"""
run_simulation.py

Runs the complete simulation-only proof of concept for:

"Spectral Silences in Raman Spectroscopy of Tooth Enamel"

Manuscript ID 2979716
Journal of Raman Spectroscopy

The script:
1. Generates the synthetic Raman-like spectrum.
2. Computes the silence scores for the candidate windows.
3. Runs the Monte-Carlo weak-feature detection benchmark.
4. Writes the numerical results to the repository's results/ directory.

Run from the repository root:

    python scripts/run_simulation.py

Outputs:

    results/silence_scores.csv
    results/table2_detection_benchmark.csv

No experimental enamel spectra are used.
All spectra and benchmark trials are synthetic.
"""

import csv
import os
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Repository paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

SRC_DIR = os.path.join(PROJECT_ROOT, "src")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")

# Allow this script to import the reusable analysis module from src/.
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)


# ---------------------------------------------------------------------------
# Import the reproducible analysis module
# ---------------------------------------------------------------------------

from spectral_silence import (
    RANDOM_SEED,
    SPECTRAL_MIN,
    SPECTRAL_MAX,
    STEP,
    K_VALUES,
    AMPLITUDES,
    CANDIDATE_WINDOWS,
    MASK_INTERVALS,
    N_TRIALS,
    generate_simulated_spectrum,
    process_spectrum,
    monte_carlo_benchmark,
)


# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------

SILENCE_SCORES_PATH = os.path.join(
    RESULTS_DIR,
    "silence_scores.csv",
)

BENCHMARK_PATH = os.path.join(
    RESULTS_DIR,
    "table2_detection_benchmark.csv",
)


# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------

def main():
    """Run the complete reproducible simulation workflow."""

    # Ensure that the repository's results directory exists.
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 72)
    print("SPECTRAL SILENCE RAMAN ENAMEL")
    print("Reproducible simulation")
    print("=" * 72)

    print("\nSimulation configuration")
    print("-" * 72)
    print(f"Random seed:       {RANDOM_SEED}")
    print(f"Spectral range:    {SPECTRAL_MIN}–{SPECTRAL_MAX} cm^-1")
    print(f"Spectral step:     {STEP} cm^-1")
    print(f"Candidate windows: {len(CANDIDATE_WINDOWS)}")
    print(f"k values:          {K_VALUES}")
    print(f"Amplitudes:        {AMPLITUDES}")
    print(f"Trials/amplitude:  {N_TRIALS}")
    print(f"Masks:             {MASK_INTERVALS}")

    # -----------------------------------------------------------------------
    # 1. Generate the base synthetic spectrum
    # -----------------------------------------------------------------------

    nu = np.arange(
        SPECTRAL_MIN,
        SPECTRAL_MAX + STEP,
        STEP,
    )

    rng = np.random.default_rng(RANDOM_SEED)

    I, baseline, M, O = generate_simulated_spectrum(
        nu,
        rng=rng,
    )

    print("\nBase synthetic spectrum generated.")
    print(f"Number of spectral points: {len(nu)}")

    # -----------------------------------------------------------------------
    # 2. Silence-score analysis
    # -----------------------------------------------------------------------

    print("\n" + "=" * 72)
    print("SILENCE-SCORE ANALYSIS")
    print("=" * 72)

    proc = process_spectrum(
        I,
        nu,
        mask_intervals=MASK_INTERVALS,
        candidate_windows=CANDIDATE_WINDOWS,
        k_values=K_VALUES,
    )

    print(
        f"\nRobust noise sigma "
        f"(first-difference MAD): {proc['sigma']:.4f}"
    )

    print("\nSilence scores S*(W):")

    mean_by_window_k2 = {}

    with open(
        SILENCE_SCORES_PATH,
        "w",
        newline="",
    ) as f:

        writer = csv.writer(f)

        writer.writerow(
            [
                "window_low",
                "window_high",
                "k",
                "S_star",
                "Q90",
            ]
        )

        for (window, k), values in proc["results"].items():

            S_value = values["S"]
            q90_value = values["q90"]

            writer.writerow(
                [
                    window[0],
                    window[1],
                    k,
                    round(S_value, 4),
                    round(q90_value, 4),
                ]
            )

            if k == 2.0:
                mean_by_window_k2[window] = S_value

            print(
                f"  Window {window}, "
                f"k={k}: "
                f"S*(W)={S_value:.3f}, "
                f"Q90={q90_value:.3f}"
            )

    best_window = max(
        mean_by_window_k2,
        key=mean_by_window_k2.get,
    )

    best_score = mean_by_window_k2[best_window]

    print(
        f"\nBest window at k=2.0: "
        f"{best_window}, "
        f"S*={best_score:.3f}"
    )

    print(
        f"\nWrote silence scores to:\n"
        f"  {SILENCE_SCORES_PATH}"
    )

    # -----------------------------------------------------------------------
    # 3. Monte-Carlo weak-feature detection benchmark
    # -----------------------------------------------------------------------

    print("\n" + "=" * 72)
    print("MONTE-CARLO DETECTION BENCHMARK")
    print("=" * 72)

    print(
        f"\nRunning Monte-Carlo benchmark: "
        f"{len(AMPLITUDES)} amplitudes × "
        f"{N_TRIALS} trials "
        f"(plus {N_TRIALS} null trials) ..."
    )

    bench = monte_carlo_benchmark(
        baseline,
        M,
        O,
        nu,
        seed=RANDOM_SEED,
    )

    # -----------------------------------------------------------------------
    # 4. Report false-positive rates
    # -----------------------------------------------------------------------

    print("\nZero-amplitude false-positive rates:")

    print(
        f"  Matched-filter test:    "
        f"{bench['mf_fpr0']:.3f}"
    )

    print(
        f"  Conventional local SNR: "
        f"{bench['snr_fpr0']:.3f}"
    )

    # -----------------------------------------------------------------------
    # 5. Report pooled AUC
    # -----------------------------------------------------------------------

    print("\nPooled AUC:")

    print(
        f"  Matched-filter test:    "
        f"{bench['auc_mf']:.3f}"
    )

    print(
        f"  Conventional local SNR: "
        f"{bench['auc_snr']:.3f}"
    )

    # -----------------------------------------------------------------------
    # 6. Write Table 2 benchmark results
    # -----------------------------------------------------------------------

    print("\nDetection probabilities by amplitude:")

    with open(
        BENCHMARK_PATH,
        "w",
        newline="",
    ) as f:

        writer = csv.writer(f)

        writer.writerow(
            [
                "injected_amplitude",
                "method",
                "detection_probability",
                "ci_95_halfwidth",
                "false_positive_rate",
                "auc",
            ]
        )

        # ---------------------------------------------------------------
        # Zero-amplitude rows
        # ---------------------------------------------------------------

        writer.writerow(
            [
                0.000,
                "Matched-filter test",
                round(bench["mf_fpr0"], 3),
                "",
                round(bench["mf_fpr0"], 3),
                round(bench["auc_mf"], 3),
            ]
        )

        writer.writerow(
            [
                0.000,
                "Conventional local SNR",
                round(bench["snr_fpr0"], 3),
                "",
                round(bench["snr_fpr0"], 3),
                round(bench["auc_snr"], 3),
            ]
        )

        # ---------------------------------------------------------------
        # Signal-present rows
        # ---------------------------------------------------------------

        for amplitude in AMPLITUDES:

            summary = bench["detection_summary"][amplitude]

            mf_probability = summary["mf_detection_prob"]
            mf_ci = summary["mf_ci_halfwidth"]

            snr_probability = summary["snr_detection_prob"]
            snr_ci = summary["snr_ci_halfwidth"]

            print(
                f"  A={amplitude}: "
                f"matched-filter={mf_probability:.3f}  "
                f"local-SNR={snr_probability:.3f}"
            )

            writer.writerow(
                [
                    amplitude,
                    "Matched-filter test",
                    round(mf_probability, 3),
                    round(mf_ci, 3),
                    "",
                    round(bench["auc_mf"], 3),
                ]
            )

            writer.writerow(
                [
                    amplitude,
                    "Conventional local SNR",
                    round(snr_probability, 3),
                    round(snr_ci, 3),
                    "",
                    round(bench["auc_snr"], 3),
                ]
            )

    print(
        f"\nWrote benchmark results to:\n"
        f"  {BENCHMARK_PATH}"
    )

    # -----------------------------------------------------------------------
    # 7. Final summary
    # -----------------------------------------------------------------------

    print("\n" + "=" * 72)
    print("SIMULATION COMPLETE")
    print("=" * 72)

    print("\nGenerated files:")

    print(
        f"  {SILENCE_SCORES_PATH}"
    )

    print(
        f"  {BENCHMARK_PATH}"
    )

    print("\nKey benchmark values:")

    print(
        f"  Matched-filter AUC:    "
        f"{bench['auc_mf']:.3f}"
    )

    print(
        f"  Local-SNR AUC:         "
        f"{bench['auc_snr']:.3f}"
    )

    print(
        f"  Matched-filter FPR:    "
        f"{bench['mf_fpr0']:.3f}"
    )

    print(
        f"  Local-SNR FPR:         "
        f"{bench['snr_fpr0']:.3f}"
    )

    print("\nDone.")


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
