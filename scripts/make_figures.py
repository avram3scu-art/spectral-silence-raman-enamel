import os
import sys
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# Allow the script to import the reusable analysis module from src/
sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "src")
    )
)

from spectral_silence import *

rng = np.random.default_rng(RANDOM_SEED)
nu = np.arange(SPECTRAL_MIN, SPECTRAL_MAX + STEP, STEP)
I, baseline, M, O = generate_simulated_spectrum(nu, rng=rng)
mask = make_global_mask(nu, MASK_INTERVALS)
b_hat = baseline_als(I)
residual = I - b_hat

# ---------- Figure 2: example simulated spectrum ----------
fig, ax = plt.subplots(figsize=(2670/300, 1302/300), dpi=300)
ax.plot(nu, I, color='#1f77b4', lw=0.5, label='Simulated spectrum')
ax.plot(nu, b_hat, color='#ff7f0e', lw=1.2, label='ALS baseline')
for a, b in MASK_INTERVALS:
    ax.axvspan(a, b, color='gray', alpha=0.25)
for W in CANDIDATE_WINDOWS:
    ax.axvspan(W[0], W[1], color='steelblue', alpha=0.10)
ax.set_xlabel('Wavenumber/cm$^{-1}$')
ax.set_ylabel('Raman Intensity/Arbitr. Units')
ax.set_title('Simulation-only proof of concept')
ax.legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig('Figure2_simulated_spectrum.png', dpi=300, facecolor='white')
plt.close()

# ---------- Figure 3: silence scores ----------
proc = process_spectrum(I, nu)
windows = CANDIDATE_WINDOWS
k_vals = K_VALUES
width = 0.25
x = np.arange(len(windows))
fig, ax = plt.subplots(figsize=(2670/300, 1393/300), dpi=300)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
for i, k in enumerate(k_vals):
    vals = [proc['results'][(W, k)]['S'] for W in windows]
    ax.bar(x + (i-1)*width, vals, width=width, label=f'k={k}', color=colors[i])
ax.set_xticks(x)
ax.set_xticklabels([f'{a}-{b}' for a, b in windows], rotation=20)
ax.set_ylim(0, 1.0)
ax.set_xlabel('Candidate window/cm$^{-1}$')
ax.set_ylabel('S*(W)')
ax.set_title('Silence scores across candidate windows')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('Figure3_silence_scores.png', dpi=300, facecolor='white')
plt.close()

# ---------- Run Monte Carlo benchmark once, reuse for Fig 4 & 5 ----------
bench = monte_carlo_benchmark(baseline, M, O, nu, seed=RANDOM_SEED)

# ---------- Figure 4: detection probability vs amplitude ----------
amps = AMPLITUDES
mf_probs = [bench['detection_summary'][A]['mf_detection_prob'] for A in amps]
mf_ci = [bench['detection_summary'][A]['mf_ci_halfwidth'] for A in amps]
snr_probs = [bench['detection_summary'][A]['snr_detection_prob'] for A in amps]
snr_ci = [bench['detection_summary'][A]['snr_ci_halfwidth'] for A in amps]

fig, ax = plt.subplots(figsize=(2220/300, 1406/300), dpi=300)
ax.errorbar(amps, mf_probs, yerr=mf_ci, marker='o', color='#1f77b4', label='Matched-filter test', capsize=3)
ax.errorbar(amps, snr_probs, yerr=snr_ci, marker='o', color='#ff7f0e', label='Conventional local SNR', capsize=3)
ax.set_ylim(0, 1.05)
ax.set_xlabel('Injected amplitude/Arbitr. Units')
ax.set_ylabel('Detection probability')
ax.set_title('Weak-feature detection benchmark')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('Figure4_detection_benchmark.png', dpi=300, facecolor='white')
plt.close()

# ---------- Figure 5: ROC comparison ----------
fpr_mf, tpr_mf = bench['roc_mf']
fpr_snr, tpr_snr = bench['roc_snr']
fig, ax = plt.subplots(figsize=(1920/300, 1615/300), dpi=300)
ax.plot(fpr_mf, tpr_mf, color='#1f77b4', lw=1.5,
        label=f"Matched-filter test (AUC={bench['auc_mf']:.3f})")
ax.plot(fpr_snr, tpr_snr, color='#ff7f0e', lw=1.5, linestyle='--',
        label=f"Conventional local SNR (AUC={bench['auc_snr']:.3f})")
ax.plot([0, 1], [0, 1], color='gray', linestyle=':', lw=1, label='Random (AUC=0.500)')
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('Pooled ROC comparison (all amplitudes)')
ax.legend(fontsize=8, loc='lower right')
plt.tight_layout()
plt.savefig('Figure5_roc_comparison.png', dpi=300, facecolor='white')
plt.close()

print("mf_fpr0", bench['mf_fpr0'], "snr_fpr0", bench['snr_fpr0'])
print("auc_mf", bench['auc_mf'], "auc_snr", bench['auc_snr'])
print("mf_probs", mf_probs, "ci", mf_ci)
print("snr_probs", snr_probs, "ci", snr_ci)
print("saved all 4 figures")
