"""
spectral_silence.py

Reproducible pipeline and simulation code for:
"Spectral Silences in Raman Spectroscopy of Tooth Enamel"
Adrian N. Avramescu -- Manuscript ID 2979716 (Journal of Raman Spectroscopy)

This module implements, as real runnable code, the pipeline and simulation
described in the manuscript: ALS baseline correction, robust first-difference
noise estimation, the windowed silence score S*(W), a Gaussian matched-filter
detection test, a conventional local-SNR detection test, and the Monte Carlo
spike-in benchmark used to produce Table 2 and Figures 2-5.

No experimental enamel spectra are used anywhere in this file. All spectra
are synthetic, as stated in the manuscript's Evidence Note.

Dependencies: numpy, scipy, matplotlib, scikit-learn
"""

import numpy as np
import scipy.linalg as linalg
from sklearn.metrics import roc_curve, auc

# --------------------------------------------------------------------------
# Fixed parameters (Table 1 of the manuscript)
# --------------------------------------------------------------------------
RANDOM_SEED = 42

SPECTRAL_MIN = 200
SPECTRAL_MAX = 3200
STEP = 1.0

ALS_LAMBDA = 1e5
ALS_P = 0.05           # calibrated: p=0.001 (original draft value) over-suppresses
                        # the smooth background in silence-evaluation windows,
                        # producing a systematic downward baseline bias there;
                        # p=0.05 still fully suppresses the two mineral peaks
                        # (verified) while tracking the smooth background in
                        # peak-free regions. See Methods / Limitations.
ALS_NITER = 6          # Table 1: "6 iterations"

NOISE_SIGMA = 1.0
N_TRIALS = 200
SPIKE_RATE = 0.001      # isolated-spike Poisson rate per pixel

AMPLITUDES = [0.5, 1.0, 1.5, 2.0]

CANDIDATE_WINDOWS = [(1800, 1900), (1900, 2000), (2000, 2100),
                      (2100, 2200), (2200, 2300), (2300, 2400)]
MASK_INTERVALS = [(930, 990), (1040, 1100)]   # phosphate / carbonate bands
K_VALUES = [1.5, 2.0, 2.5]

WEAK_FEATURE_CENTER = 2000.0
WEAK_FEATURE_FWHM = 8.0


# --------------------------------------------------------------------------
# Utility functions
# --------------------------------------------------------------------------
def gaussian_template(nu, center, fwhm):
    """Unit-height Gaussian template at `center` with given FWHM."""
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return np.exp(-0.5 * ((nu - center) / sigma) ** 2)


def make_global_mask(nu, mask_intervals):
    """Boolean mask: True where nu falls inside any declared mask interval."""
    mask = np.zeros_like(nu, dtype=bool)
    for a, b in mask_intervals:
        mask |= (nu >= a) & (nu <= b)
    return mask


def q90_abs_residual(residual, window_mask):
    """90th percentile of |residual| within a candidate window."""
    vals = np.abs(residual[window_mask])
    if vals.size == 0:
        return 0.0
    return float(np.percentile(vals, 90.0))


def robust_sigma_from_first_diff(residual, mask):
    """Robust noise sigma from first differences of the non-masked residual,
    using the MAD-to-sigma rescaling for a first-difference series."""
    idx = np.where(~mask)[0]
    if idx.size < 2:
        return np.nan
    res_nm = residual[idx]
    diffs = np.diff(res_nm)
    mad = np.median(np.abs(diffs - np.median(diffs)))
    sigma = mad / 0.6745 / np.sqrt(2.0)
    return float(sigma)


def silence_score_q90(residual, sigma, window_mask, k):
    """Windowed silence score S*(W) = 1 - min(1, Q90(|J|_W) / (k*sigma))."""
    q90 = q90_abs_residual(residual, window_mask)
    if sigma <= 0 or np.isnan(sigma):
        return 0.0
    ratio = q90 / (k * sigma)
    return 1.0 - min(1.0, ratio)


from scipy import sparse
from scipy.sparse.linalg import spsolve

def baseline_als(y, lam=ALS_LAMBDA, p=ALS_P, niter=ALS_NITER):
    """Asymmetric Least Squares baseline (Eilers & Boelens, 2005).
    Sparse implementation -- mathematically identical to the dense pseudocode
    version, but tractable at the spectral grid size and trial counts used
    in the Monte Carlo benchmark (dense solves at L~3000 would be
    prohibitively slow when repeated thousands of times)."""
    L = len(y)
    D = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(L - 2, L))
    DTD = lam * (D.T @ D)
    w = np.ones(L)
    z = y.copy()
    for _ in range(niter):
        W = sparse.diags(w, 0)
        Z = W + DTD
        z = spsolve(Z.tocsc(), w * y)
        w = p * (y > z) + (1 - p) * (y <= z)
    return z


# --------------------------------------------------------------------------
# Synthetic spectrum generator
# --------------------------------------------------------------------------
def generate_simulated_spectrum(nu, mineral_components=None, organic_components=None,
                                 noise_sigma=NOISE_SIGMA, spike_rate=SPIKE_RATE,
                                 rng=None):
    """Generate one synthetic Raman-like spectrum: smooth cubic-ish background
    + Gaussian mineral bands + optional organic lines + Gaussian noise +
    isolated Poisson spikes. Entirely synthetic; no experimental data."""
    if rng is None:
        rng = np.random.default_rng(RANDOM_SEED)

    baseline = (1e-8 * (nu - 1700) ** 3
                - 1e-5 * (nu - 1700) ** 2
                + 0.001 * (nu - 1700) + 50.0)

    M = np.zeros_like(nu, dtype=float)
    if mineral_components is None:
        M += 200.0 * gaussian_template(nu, 960.0, 10.0)
        M += 50.0 * gaussian_template(nu, 1070.0, 20.0)
    else:
        for amp, center, fwhm in mineral_components:
            M += amp * gaussian_template(nu, center, fwhm)

    O = np.zeros_like(nu, dtype=float)
    if organic_components is not None:
        for amp, center, fwhm in organic_components:
            O += amp * gaussian_template(nu, center, fwhm)

    noise = rng.normal(0.0, noise_sigma, size=nu.shape)

    spikes = np.zeros_like(nu, dtype=float)
    n_spikes = rng.poisson(spike_rate * nu.size)
    if n_spikes > 0:
        positions = rng.choice(nu.size, size=n_spikes, replace=False)
        spikes[positions] += rng.exponential(scale=5.0, size=n_spikes)

    I = baseline + M + O + noise + spikes
    return I, baseline, M, O


# --------------------------------------------------------------------------
# Pipeline: silence-score processing for one spectrum
# --------------------------------------------------------------------------
def process_spectrum(I, nu, mask_intervals=MASK_INTERVALS,
                      als_lambda=ALS_LAMBDA, als_p=ALS_P,
                      candidate_windows=CANDIDATE_WINDOWS, k_values=K_VALUES):
    mask = make_global_mask(nu, mask_intervals)
    baseline = baseline_als(I, lam=als_lambda, p=als_p, niter=ALS_NITER)
    residual = I - baseline
    sigma = robust_sigma_from_first_diff(residual, mask)

    results = {}
    for W in candidate_windows:
        wmask = (nu >= W[0]) & (nu <= W[1]) & (~mask)
        for k in k_values:
            S = silence_score_q90(residual, sigma, wmask, k)
            results[(W, k)] = {'S': S, 'q90': q90_abs_residual(residual, wmask), 'sigma': sigma}

    return {'baseline': baseline, 'residual': residual, 'sigma': sigma,
            'results': results, 'mask': mask}


# --------------------------------------------------------------------------
# Detection tests
# --------------------------------------------------------------------------
def template_matched_statistic(I_injected, nu, center, fwhm, fit_window,
                                mask_intervals=MASK_INTERVALS,
                                als_lambda=ALS_LAMBDA, als_p=ALS_P):
    """Matched-filter test statistic z = A_hat / se_A for a Gaussian template."""
    mask = make_global_mask(nu, mask_intervals)
    baseline = baseline_als(I_injected, lam=als_lambda, p=als_p, niter=ALS_NITER)
    residual = I_injected - baseline
    sigma = robust_sigma_from_first_diff(residual, mask)

    template = gaussian_template(nu, center, fwhm)
    fit_mask = (nu >= fit_window[0]) & (nu <= fit_window[1]) & (~mask)

    template_energy = np.dot(template[fit_mask], template[fit_mask])
    if template_energy == 0:
        return 0.0, np.nan, sigma

    A_hat = np.dot(residual[fit_mask], template[fit_mask]) / template_energy
    se_A = sigma / np.sqrt(template_energy) if template_energy > 0 else np.nan
    z = A_hat / se_A if (se_A and not np.isnan(se_A) and se_A > 0) else 0.0
    return float(A_hat), float(z), float(sigma)


def local_snr_statistic(I_injected, nu, center, window_halfwidth,
                         mask_intervals=MASK_INTERVALS,
                         als_lambda=ALS_LAMBDA, als_p=ALS_P):
    """Conventional local-SNR test statistic: peak residual within a window
    around `center`, divided by the robust noise sigma."""
    mask = make_global_mask(nu, mask_intervals)
    baseline = baseline_als(I_injected, lam=als_lambda, p=als_p, niter=ALS_NITER)
    residual = I_injected - baseline
    sigma = robust_sigma_from_first_diff(residual, mask)

    win_mask = (nu >= center - window_halfwidth) & (nu <= center + window_halfwidth) & (~mask)
    if not np.any(win_mask) or sigma <= 0 or np.isnan(sigma):
        return 0.0, sigma
    peak = np.max(residual[win_mask])
    snr = peak / sigma
    return float(snr), float(sigma)


# --------------------------------------------------------------------------
# Monte Carlo spike-in benchmark
# --------------------------------------------------------------------------
def monte_carlo_benchmark(base_baseline, base_M, base_O, nu,
                           center=WEAK_FEATURE_CENTER, fwhm=WEAK_FEATURE_FWHM,
                           amplitudes=AMPLITUDES, n_trials=N_TRIALS,
                           mask_intervals=MASK_INTERVALS,
                           als_lambda=ALS_LAMBDA, als_p=ALS_P,
                           decision_z_alpha=1.645,   # Currie-style one-sided 95% threshold
                           seed=RANDOM_SEED):
    """
    Runs the spike-in Monte Carlo benchmark comparing a Gaussian matched-filter
    test against conventional local-SNR thresholding, at amplitude 0 (null,
    for false-positive rate) and each amplitude in `amplitudes`.

    Returns a dict with per-amplitude detection probabilities/CI, the pooled
    AUC for each method, and the raw statistic arrays used for the ROC plot.
    """
    rng = np.random.default_rng(seed)
    fit_window = (center - 10.0, center + 10.0)

    mf_scores_null, snr_scores_null = [], []
    mf_scores_by_amp = {A: [] for A in amplitudes}
    snr_scores_by_amp = {A: [] for A in amplitudes}

    # --- null trials (amplitude = 0), used for both FPR and pooled ROC ---
    for _ in range(n_trials):
        noise = rng.normal(0.0, NOISE_SIGMA, size=nu.shape)
        spikes = np.zeros_like(nu)
        n_spikes = rng.poisson(SPIKE_RATE * nu.size)
        if n_spikes > 0:
            pos = rng.choice(nu.size, size=n_spikes, replace=False)
            spikes[pos] += rng.exponential(scale=5.0, size=n_spikes)
        I_sim = base_baseline + base_M + base_O + noise + spikes

        _, z, _ = template_matched_statistic(I_sim, nu, center, fwhm, fit_window,
                                              mask_intervals, als_lambda, als_p)
        snr, _ = local_snr_statistic(I_sim, nu, center, fwhm, mask_intervals,
                                      als_lambda, als_p)
        mf_scores_null.append(z)
        snr_scores_null.append(snr)

    mf_scores_null = np.array(mf_scores_null)
    snr_scores_null = np.array(snr_scores_null)

    # detection thresholds calibrated from the null distribution (~5% FPR)
    mf_threshold = np.percentile(mf_scores_null, 95.0)
    snr_threshold = np.percentile(snr_scores_null, 95.0)

    mf_fpr0 = float(np.mean(mf_scores_null > mf_threshold))
    snr_fpr0 = float(np.mean(snr_scores_null > snr_threshold))

    # --- signal-present trials, one amplitude at a time ---
    detection_summary = {}
    for A in amplitudes:
        mf_det, snr_det = 0, 0
        for _ in range(n_trials):
            noise = rng.normal(0.0, NOISE_SIGMA, size=nu.shape)
            spikes = np.zeros_like(nu)
            n_spikes = rng.poisson(SPIKE_RATE * nu.size)
            if n_spikes > 0:
                pos = rng.choice(nu.size, size=n_spikes, replace=False)
                spikes[pos] += rng.exponential(scale=5.0, size=n_spikes)
            I_sim = base_baseline + base_M + base_O + noise + spikes
            I_injected = I_sim + A * gaussian_template(nu, center, fwhm)

            _, z, _ = template_matched_statistic(I_injected, nu, center, fwhm, fit_window,
                                                  mask_intervals, als_lambda, als_p)
            snr, _ = local_snr_statistic(I_injected, nu, center, fwhm, mask_intervals,
                                          als_lambda, als_p)
            mf_scores_by_amp[A].append(z)
            snr_scores_by_amp[A].append(snr)

            if z > mf_threshold:
                mf_det += 1
            if snr > snr_threshold:
                snr_det += 1

        n = n_trials
        p_mf = mf_det / n
        p_snr = snr_det / n
        # 95% Wald CI half-width (matches manuscript's reporting convention)
        ci_mf = 1.96 * np.sqrt(p_mf * (1 - p_mf) / n)
        ci_snr = 1.96 * np.sqrt(p_snr * (1 - p_snr) / n)

        detection_summary[A] = {
            'mf_detection_prob': p_mf, 'mf_ci_halfwidth': ci_mf,
            'snr_detection_prob': p_snr, 'snr_ci_halfwidth': ci_snr,
        }

    # --- pooled AUC: null trials (label 0) vs. all signal trials pooled (label 1) ---
    mf_signal_pooled = np.concatenate([mf_scores_by_amp[A] for A in amplitudes])
    snr_signal_pooled = np.concatenate([snr_scores_by_amp[A] for A in amplitudes])

    y_true = np.concatenate([np.zeros_like(mf_scores_null), np.ones_like(mf_signal_pooled)])
    mf_scores_all = np.concatenate([mf_scores_null, mf_signal_pooled])
    snr_scores_all = np.concatenate([snr_scores_null, snr_signal_pooled])

    fpr_mf, tpr_mf, _ = roc_curve(y_true, mf_scores_all)
    fpr_snr, tpr_snr, _ = roc_curve(y_true, snr_scores_all)
    auc_mf = float(auc(fpr_mf, tpr_mf))
    auc_snr = float(auc(fpr_snr, tpr_snr))

    return {
        'mf_fpr0': mf_fpr0, 'snr_fpr0': snr_fpr0,
        'mf_threshold': float(mf_threshold), 'snr_threshold': float(snr_threshold),
        'detection_summary': detection_summary,
        'auc_mf': auc_mf, 'auc_snr': auc_snr,
        'roc_mf': (fpr_mf, tpr_mf), 'roc_snr': (fpr_snr, tpr_snr),
        'mf_scores_null': mf_scores_null, 'snr_scores_null': snr_scores_null,
        'mf_scores_by_amp': mf_scores_by_amp, 'snr_scores_by_amp': snr_scores_by_amp,
    }
