Reproducibility Guide

Spectral Silences in Raman Spectroscopy of Tooth Enamel

This repository contains the code, simulation outputs, and figures supporting the methodological study:

“Spectral Silences in Raman Spectroscopy of Tooth Enamel: A reproducible framework for low-signal window identification, baseline anchoring, noise estimation and weak-feature detectability”

The repository is intended to make the computational proof-of-concept transparent and reproducible.

Important scope statement

This repository does not contain experimental enamel Raman spectra.

All numerical results presented in the computational analysis are generated from a fully reproducible synthetic simulation. The simulation is intended to demonstrate the mathematical and computational behavior of the proposed framework.

No universal Raman “silent window” for tooth enamel is claimed on the basis of these simulations.

Experimental validation using real enamel spectra remains a separate requirement.

Repository structure

src/
    spectral_silence.py

scripts/
    run_simulation.py
    make_figures.py

results/
    silence_scores.csv
    table2_detection_benchmark.csv

figures/
    Figure2_simulated_spectrum.png
    Figure3_silence_scores.png
    Figure4_detection_benchmark.png
    Figure5_roc_comparison.png

docs/
    reproducibility.md

Simulation configuration

The synthetic simulation uses the following fixed parameters:

Parameter| Specification
Spectral grid| 200–3200 cm⁻¹ at 1 cm⁻¹
Masks| 930–990 and 1040–1100 cm⁻¹
Baseline| Asymmetric least squares (ALS)
ALS λ| 1 × 10⁵
ALS p| 0.05
ALS iterations| 6
Gaussian noise σ| 1
Isolated-spike rate| 0.001 per sample
Random seed| 42
Candidate windows| 1800–1900, 1900–2000, 2000–2100, 2100–2200, 2200–2300, 2300–2400 cm⁻¹
Weak feature| Gaussian
Weak-feature centre| 2000 cm⁻¹
Weak-feature FWHM| 8 cm⁻¹
Trials| 200 Monte-Carlo trials per amplitude

The weak-feature amplitudes evaluated in the simulation are 0.5, 1.0, 1.5 and 2.0 arbitrary units.

Analytical workflow

The computational workflow consists of:

1. Constructing the synthetic Raman-like spectrum.
2. Applying the predefined spectral masks.
3. Estimating the baseline using fixed ALS parameters.
4. Subtracting the estimated baseline to obtain the residual spectrum.
5. Estimating the noise level from first differences using a robust MAD estimator.
6. Evaluating candidate spectral windows using the silence score.
7. Injecting controlled weak Gaussian features.
8. Comparing matched-filter detection with conventional local SNR.
9. Generating numerical result tables and figures.

Silence score

For a candidate window W, the residual spectrum is defined as:

[
J(\nu)=I(\nu)-\hat{B}(\nu)
]

The robust noise estimate is obtained from first differences:

[
\sigma =
\frac{\operatorname{MAD}(\Delta J)}
{0.6745\sqrt{2}}
]

The windowed silence score is:

[
S^*(W)

1-
\min
\left{
1,
\frac{Q_{90}(|J|_W)}
{k\sigma}
\right}
]

where Q_{90} is the 90th percentile of the absolute residual amplitude within the candidate window and k is the specified tolerance parameter.

The simulation evaluates:

k = 1.5
k = 2.0
k = 2.5

A value approaching 1 indicates greater analytical quietness relative to the specified noise model. A value approaching 0 indicates that residual structure reaches or exceeds the specified tolerance.

Running the simulation

The main simulation is implemented in:

scripts/run_simulation.py

The reusable analysis functions are contained in:

src/spectral_silence.py

The figure-generation workflow is contained in:

scripts/make_figures.py

The exact execution commands should be checked against the repository's current script imports and dependencies before being treated as the final command-line protocol.

Simulation outputs

The repository preserves the numerical outputs in:

results/silence_scores.csv
results/table2_detection_benchmark.csv

The corresponding visual outputs are:

figures/Figure2_simulated_spectrum.png
figures/Figure3_silence_scores.png
figures/Figure4_detection_benchmark.png
figures/Figure5_roc_comparison.png

Reported computational benchmark

The simulation reported a pooled AUC of:

- 0.925 for the Gaussian matched-filter test.
- 0.800 for conventional local SNR.

Both methods used thresholds calibrated from independent zero-amplitude simulation trials at the 95th percentile, targeting a nominal 5% false-positive rate. The empirical zero-amplitude false-positive rate was 0.050 for both methods.

These values describe the synthetic simulation and should not be interpreted as experimental performance on real enamel spectra.

Interpretation

The purpose of the repository is to allow another researcher to inspect the implementation, reproduce the synthetic calculations, examine the generated outputs, and evaluate the assumptions of the proposed framework.

The simulation demonstrates computational behavior under a controlled synthetic model. It does not establish that any particular spectral interval is silent in real tooth enamel.

Future experimental validation

A future experimental validation should use real enamel spectra from multiple specimens with technical replicates, acquisition metadata, explicit artifact handling, and a predefined processing pipeline.

The experimental study should compare the proposed silence score with alternative measures such as RMS/Q95 and conventional SNR and should evaluate controlled weak-feature injections, sensitivity, specificity, false-positive rate, detection bias and ROC AUC.

Baseline parameters, mask boundaries, tolerance k, window width, spectral resolution, calibration uncertainty and alternative noise estimators should also be stress-tested.

Multi-instrument repeatability is particularly important because spectral silence is defined relative to a measurement and preprocessing pipeline.

Reproducibility principle

The central reproducibility principle of this repository is:

«The computational result must be traceable from the documented parameters and source code to the numerical outputs and figures.»

No experimental conclusion should be inferred from the simulation beyond what the synthetic model directly demonstrates.
