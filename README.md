Spectral Silences in Raman Spectroscopy of Tooth Enamel

Reproducible simulation and analysis framework for low-signal window identification, baseline anchoring, noise estimation, and weak-feature detectability

Author: Adrian N. Avramescu
Manuscript: Spectral Silences in Raman Spectroscopy of Tooth Enamel
Journal: Journal of Raman Spectroscopy
Manuscript ID: 2979716

---

Overview

This repository contains the computational materials supporting the methodological study:

«Spectral Silences in Raman Spectroscopy of Tooth Enamel: A reproducible framework for low-signal window identification, baseline anchoring, noise estimation and weak-feature detectability»

The project defines spectral silence as an operational, measurement-specific property of a Raman spectrum and its documented preprocessing pipeline.

The framework combines:

- explicit spectral masking;
- fixed asymmetric least-squares baseline estimation;
- robust noise estimation from first differences;
- a windowed silence score S^*(W);
- controlled weak-feature injection;
- Gaussian matched-filter detection;
- comparison with conventional local SNR.

---

Important scope statement

This repository contains no experimental enamel Raman spectra.

All numerical results in the computational proof-of-concept are generated from an explicitly defined synthetic simulation.

The simulation demonstrates the behavior of the proposed analytical framework under controlled conditions. It does not establish a universal silent spectral interval for tooth enamel, nor does it demonstrate experimental detection performance on real enamel.

Experimental validation remains a separate requirement.

---

Repository contents

spectral-silence-raman-enamel/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── src/
│   └── spectral_silence.py
│
├── scripts/
│   ├── run_simulation.py
│   └── make_figures.py
│
├── results/
│   ├── silence_scores.csv
│   └── table2_detection_benchmark.csv
│
├── figures/
│   ├── Figure2_simulated_spectrum.png
│   ├── Figure3_silence_scores.png
│   ├── Figure4_detection_benchmark.png
│   └── Figure5_roc_comparison.png
│
└── docs/
    └── reproducibility.md

---

Note on committed outputs

The CSV files in results/ and the PNG files in figures/ are committed to this repository intentionally, not as build artifacts. They are the exact outputs produced by scripts/run_simulation.py and scripts/make_figures.py at the time of manuscript submission, included so a reader can inspect the reported results without first installing dependencies or running any code. Re-running the scripts (see Reproducibility below) regenerates these same files from scratch.

---

Simulation

The synthetic spectrum uses a spectral grid from 200–3200 cm⁻¹ at 1 cm⁻¹ resolution.

The simulation includes:

- representative mineral-like Gaussian bands;
- smooth baseline curvature;
- Gaussian stochastic noise;
- isolated spikes;
- predefined spectral masks;
- controlled weak-feature injections.

Fixed parameters

Parameter| Value
Spectral range| 200–3200 cm⁻¹
Spectral spacing| 1 cm⁻¹
Masks| 930–990 and 1040–1100 cm⁻¹
Baseline| ALS
ALS λ| 1\times10^5
ALS p| 0.05
ALS iterations| 6
Gaussian noise σ| 1
Spike rate| 0.001/sample
Random seed| 42
Weak-feature centre| 2000 cm⁻¹
Weak-feature FWHM| 8 cm⁻¹
Trials| 200 per amplitude

Candidate windows are:

1800–1900 cm⁻¹
1900–2000 cm⁻¹
2000–2100 cm⁻¹
2100–2200 cm⁻¹
2200–2300 cm⁻¹
2300–2400 cm⁻¹

The tolerance parameter k is evaluated at:

k = 1.5
k = 2.0
k = 2.5

---

Silence score

After baseline estimation, the residual spectrum is:

[
J(\nu)=I(\nu)-\hat{B}(\nu)
]

Noise is estimated from first differences using:

[
\sigma =
\frac{\operatorname{MAD}(\Delta J)}
{0.6745\sqrt{2}}
]

For a candidate window W:

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

where Q_{90} is the 90th percentile of absolute residual amplitude within the window.

Higher S^*(W) indicates greater analytical quietness relative to the specified noise model.

The score is not a claim of molecular absence.

---

Detection benchmark

The simulation compares:

1. a Gaussian matched-filter test;
2. conventional local SNR.

Decision thresholds were calibrated using independent zero-amplitude simulation trials at the 95th percentile, targeting a nominal 5% false-positive rate.

The reported simulation results were:

Method| AUC| Empirical zero-amplitude false-positive rate
Gaussian matched filter| 0.925| 0.050
Conventional local SNR| 0.800| 0.050

Matched-filter detection probabilities at injected amplitudes of 0.5, 1.0, 1.5 and 2.0 arbitrary units were 0.315, 0.685, 0.945 and 0.995, respectively.

For conventional local SNR, the corresponding detection probabilities were 0.105, 0.240, 0.500 and 0.690.

These are simulation-only results and should not be interpreted as measurements of real enamel.

---

Figures

The repository contains the four simulation figures corresponding to the manuscript:

- "Figure2_simulated_spectrum.png" — synthetic Raman-like spectrum;
- "Figure3_silence_scores.png" — simulation-only silence scores;
- "Figure4_detection_benchmark.png" — weak-feature detection benchmark;
- "Figure5_roc_comparison.png" — ROC comparison of the two detection methods.

No figure contains experimental enamel data.

---

Reproducibility

The computational implementation is organized as follows:

Core analysis

src/spectral_silence.py

Contains the reusable analysis functions.

Simulation

scripts/run_simulation.py

Runs the synthetic simulation and produces numerical outputs.

Figure generation

scripts/make_figures.py

Generates the simulation figures from the computational workflow.

Numerical outputs

results/

Contains the CSV outputs used for the reported computational results.

Reproducibility documentation

docs/reproducibility.md

Contains the detailed simulation configuration, analytical workflow, scope statement, and future validation requirements.

---

Experimental validation

The next stage of the research is experimental validation using real enamel Raman spectra.

A suitable validation study should include:

- multiple enamel specimens;
- technical replicates;
- raw or minimally processed spectra;
- complete acquisition metadata;
- explicit artifact handling;
- predefined masks and baseline parameters;
- evaluation of S^*(W);
- comparison with alternative noise and signal measures;
- controlled weak-feature spike-in experiments;
- sensitivity and specificity;
- false-positive rate;
- detection bias;
- ROC AUC;
- stress testing of baseline parameters, k, window width, spectral resolution, and noise estimation;
- multi-instrument repeatability.

The purpose of this repository is therefore to make the computational component reproducible while clearly identifying the experimental questions that remain open.

---

Scientific scope

The central methodological proposition is that spectral silence should be treated as a reproducible analytical object rather than a visual impression or a claim of molecular absence.

The repository therefore separates:

Demonstrated computationally

- implementation of the proposed silence metric;
- behavior of the metric under a controlled synthetic spectrum;
- controlled weak-feature injection;
- matched-filter versus local-SNR benchmarking.

Not demonstrated by this repository

- a universal silent interval in real tooth enamel;
- experimental enamel spectra;
- clinical performance;
- in vivo performance;
- safety;
- universal superiority of the proposed statistic.

---

Citation

If you use this repository, please cite the associated manuscript:

Avramescu, A. N. Spectral Silences in Raman Spectroscopy of Tooth Enamel: A reproducible framework for low-signal window identification, baseline anchoring, noise estimation and weak-feature detectability.

Journal of Raman Spectroscopy. Manuscript ID 2979716.

---

License

The source code in this repository is distributed under the license specified in "LICENSE".

Please note that the manuscript and associated scholarly text may be subject to different copyright or licensing terms from the source code.

---

Contact

Adrian N. Avramescu
Independent Researcher
Corresponding author: avram3scu@gmail.com
