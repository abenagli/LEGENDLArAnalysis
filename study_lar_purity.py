#!/usr/bin/env python3
"""LAr purity from the lifetime of the slow (triplet) scintillation component.

The raw CAEN wavedump files (<run>/0_wave{0,1}.dat) are read directly: the
reconstructed analyzed*.root files keep only a fraction of the waveforms.

Per event (sum of the two channels): baseline from the pre-trigger, t50 of the
leading edge, charge and fprompt in windows relative to the t50. Events that
are saturated (front-end clipping), on a sloped baseline, or with a second
pulse in the tail are rejected. Electron-like events (low fprompt, high
charge: muons, gammas, 39Ar) are aligned at the t50, normalized to their
charge and averaged, per run and per dataset.

The average waveform is fitted from TAIL_FIT_START after the t50 with
    A_T exp(-t/tau_T) + A_S exp(-t/TAU_SPE) + c
where the second term is the tail of the single-photon response of the
SiPM + front-end (tau ~ 0.6 us) and c a constant offset. Variants of the fit
(start, TAU_SPE, single exponential at late times) give the systematic error.
tau_T is converted to an equivalent N2 (O2) contamination with
    1/tau_T = 1/tau_T0 + k_Q [X].
"""

import argparse
import calendar
import csv
import datetime
import importlib.util
import math
import re
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


SCRIPT_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = SCRIPT_DIR / "config_purity.py"

# Raw data format (CAEN wavedump, binary with header).
HEADER_WORDS = 6
DT_NS = 4.0
ADC_MAX = 16383
CHANNEL_FILES = ("0_wave0.dat", "0_wave1.dat")
CHANNELS = ("ch0", "ch1", "sum")

# Event processing, in samples.
BASELINE_RANGE = (100, 900)
T50_SEARCH = (900, 1700)
PRE_T50 = 10                  # windows start 40 ns before the t50
FAST_SAMPLES = 150            # 600 ns, as Read::fast
ALIGN_SAMPLE = 1000           # t50 of the aligned waveforms
REBIN = 10                    # 40 ns bins for the fit
CHUNK_EVENTS = 2000

# The front-end clips at ~4400 ADC above the baseline (not at the 14-bit
# full scale): events reaching this level in either channel are rejected.
SATURATION_ABOVE_BASELINE = 4300.0
# Events on the tail of a previous pulse: baseline slope in the pre-trigger.
BASELINE_SLOPE_MAX_ADC = 2.0
# Pile-up: a fast rise (over RISE_SAMPLES) in the tail larger than
# PILEUP_RISE_FRACTION of the rise of the main pulse, from 2 us after the t50.
RISE_SAMPLES = 5
PILEUP_START = 500
PILEUP_RISE_FRACTION = 0.05
PILEUP_RISE_MIN_ADC = 30.0

# Tail fit (times in ns from the t50).
TAIL_FIT_START = 1500.0
TAIL_FIT_END = 15800.0
TAU_SPE = 600.0
SYSTEMATIC_VARIANTS = (
    {"start": 1000.0},
    {"start": 2500.0},
    {"tau_spe": 450.0},
    {"tau_spe": 750.0},
    {"start": 4000.0, "tau_spe": None},
)

# Electron-like band: peak of the fprompt (fast/total) distribution of the
# candidates with charge > charge_min, searched in ER_BAND_BINS = (bins, min,
# max). Nominal window [peak - 2 sigma, peak + 1 sigma], asymmetric because
# alpha-like and long-tail populations lie at higher fprompt; the variants give
# the selection systematic.
ER_BAND_BINS = (90, 0.05, 0.50)
ER_BAND_MIN_SIGMA = 0.005
# The electron-like peak is searched below this value (observed peaks: 0.29-0.38
# depending on the setup; alpha-like events are at ~0.5).
ER_PEAK_SEARCH_MAX = 0.42
SELECTION_VARIANTS = {
    "nominal": {"fprompt_sigmas": (2.0, 1.0)},
    "fp_tight": {"fprompt_sigmas": (1.0, 0.5)},
    "fp_loose": {"fprompt_sigmas": (3.0, 2.0)},
    "high_charge": {"fprompt_sigmas": (2.0, 1.0), "charge_factor": 3.0},
}

# Diagnostic 2D histogram of the event selection.
LOG_CHARGE_BINS = (80, 4.0, 8.0)
FPROMPT_BINS = (50, 0.0, 1.0)

# Summary of all runs versus acquisition date: one panel per data-taking
# period (a new period starts after a gap longer than PERIOD_GAP_DAYS).
MUON_RUNS = ("run11", "run30", "run111")
PERIOD_GAP_DAYS = 3.0


def load_config(config_file):
    config_path = Path(config_file).expanduser().resolve()
    spec = importlib.util.spec_from_file_location("purity_config", config_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load config file: {config_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "PURITY_CONFIG"):
        raise AttributeError(f"{config_path} does not define PURITY_CONFIG")
    return module.PURITY_CONFIG


# ---------------------------------------------------------------- raw data --

def open_channel(path):
    header = np.fromfile(path, dtype=np.uint32, count=HEADER_WORDS)
    if header.size < HEADER_WORDS:
        return None, 0, 0
    event_bytes = int(header[0])
    n_samples = (event_bytes - 4 * HEADER_WORDS) // 2
    n_events = path.stat().st_size // event_bytes
    if n_events == 0:
        return None, 0, 0
    data = np.memmap(path, dtype=np.uint8, mode="r",
                     shape=(n_events, event_bytes))
    return data, n_events, n_samples


def find_runs(dataset_folder, exclude=()):
    runs = []
    for run_dir in Path(dataset_folder).iterdir():
        if not run_dir.is_dir() or not run_dir.name.startswith("run"):
            continue
        if run_dir.name in exclude:
            continue
        if all((run_dir / name).is_file() for name in CHANNEL_FILES):
            runs.append(run_dir)
    return sorted(runs, key=lambda d: (len(d.name), d.name))


def window_sum(cumulative, start, length):
    stop = np.minimum(start + length, cumulative.shape[1] - 1)
    rows = np.arange(cumulative.shape[0])
    return cumulative[rows, stop] - cumulative[rows, start]


def find_er_band(fprompt, minimum_events=50):
    """Electron-like peak in the fprompt distribution of the candidates.

    The peak position is the mode; its width is taken from the half width at
    half maximum on the low side, which is free from other populations.
    """
    if fprompt.size < minimum_events:
        return None
    counts, edges = np.histogram(fprompt, bins=ER_BAND_BINS[0], range=ER_BAND_BINS[1:])
    smooth = np.convolve(counts, np.ones(5) / 5.0, mode="same")
    searchable = 0.5 * (edges[1:] + edges[:-1]) < ER_PEAK_SEARCH_MAX
    peak = int(np.where(searchable, smooth, -1.0).argmax())
    half = 0.5 * smooth[peak]
    low = peak
    while low > 0 and smooth[low] > half:
        low -= 1
    width = ER_BAND_BINS[2] - ER_BAND_BINS[1]
    sigma = max((peak - low) * width / ER_BAND_BINS[0] / 1.1774, ER_BAND_MIN_SIGMA)
    centre = 0.5 * (edges[peak] + edges[peak + 1])
    return centre, sigma


def new_accumulator(n_bins):
    return {"n": 0, "rsum": np.zeros((3, n_bins)), "rsum2": np.zeros((3, n_bins)),
            "rcount": np.zeros((3, n_bins))}


def accumulate(acc, rebinned, rvalid):
    acc["n"] += rebinned.shape[0]
    acc["rsum"] += np.where(rvalid, rebinned, 0.0).sum(axis=0)
    acc["rsum2"] += np.where(rvalid, rebinned ** 2, 0.0).sum(axis=0)
    acc["rcount"] += rvalid.sum(axis=0)


def process_run(run_dir, selection):
    """Average the aligned, charge-normalized waveforms of electron-like events.

    Returns the accumulators of the nominal selection and of its variants
    (SELECTION_VARIANTS), the fprompt window used and the run statistics.
    """
    channels = [open_channel(run_dir / name) for name in CHANNEL_FILES]
    if any(data is None for data, _, _ in channels):
        return None
    n_samples = channels[0][2]
    if any(ns != n_samples for _, _, ns in channels) or n_samples < 4000:
        warnings.warn(f"{run_dir}: unsupported record length; skipping")
        return None
    n_events = min(n for _, n, _ in channels)
    n_bins = n_samples // REBIN
    columns = np.arange(n_samples)

    hist2d = np.zeros((LOG_CHARGE_BINS[0], FPROMPT_BINS[0]))
    stats = {"n_events": n_events, "saturated": 0, "bad_baseline": 0,
             "pileup": 0, "candidates": 0}
    # Candidates: aligned waveforms rebinned to 40 ns, kept in memory until
    # the electron-like band is known.
    cand_wvf, cand_valid, cand_fprompt, cand_charge = [], [], [], []

    for start in range(0, n_events, CHUNK_EVENTS):
        stop = min(start + CHUNK_EVENTS, n_events)
        n = stop - start
        saturated = np.zeros(n, dtype=bool)
        sloped = np.zeros(n, dtype=bool)
        signals = []
        baseline_rms = []
        half = (BASELINE_RANGE[0] + BASELINE_RANGE[1]) // 2
        for data, _, _ in channels:
            raw = np.asarray(data[start:stop, 4 * HEADER_WORDS:]).view(np.uint16)
            wvf = raw.reshape(n, n_samples).astype(np.float32)
            pre = wvf[:, BASELINE_RANGE[0]:BASELINE_RANGE[1]]
            signal = wvf - pre.mean(axis=1, keepdims=True)
            saturated |= (wvf.max(axis=1) >= ADC_MAX) | (
                signal.max(axis=1) >= SATURATION_ABOVE_BASELINE)
            slope = (wvf[:, BASELINE_RANGE[0]:half].mean(axis=1)
                     - wvf[:, half:BASELINE_RANGE[1]].mean(axis=1))
            sloped |= np.abs(slope) > BASELINE_SLOPE_MAX_ADC
            signals.append(signal)
            baseline_rms.append(pre.std(axis=1))
        signals.append(signals[0] + signals[1])
        total = signals[2]
        baseline_rms = np.maximum(*baseline_rms)
        good_baseline = ~sloped & (
            baseline_rms < selection["baseline_rms_factor"] * np.median(baseline_rms))

        window = total[:, T50_SEARCH[0]:T50_SEARCH[1]]
        amplitude = window.max(axis=1)
        t50 = (window >= 0.5 * amplitude[:, None]).argmax(axis=1) + T50_SEARCH[0]
        first = np.maximum(t50 - PRE_T50, 0)

        cumulative = np.zeros((n, n_samples + 1))
        np.cumsum(total, axis=1, out=cumulative[:, 1:])
        # fprompt = fast / total over the whole record: with a finite slow
        # window (Read::slow) the cut would select events by the number of
        # photons in part of the tail and bias the average shape.
        total_sum = window_sum(cumulative, first, n_samples)
        charge = total_sum * DT_NS
        fast = window_sum(cumulative, first, FAST_SAMPLES)
        with np.errstate(divide="ignore", invalid="ignore"):
            fprompt = np.where(total_sum > 0, fast / total_sum, -1.0)

        rise = total[:, RISE_SAMPLES:] - total[:, :-RISE_SAMPLES]
        main_rise = rise[:, T50_SEARCH[0]:T50_SEARCH[1]].max(axis=1)
        tail_first = np.minimum(t50 + PILEUP_START, rise.shape[1] - 1)
        tail_mask = np.arange(rise.shape[1])[None, :] >= tail_first[:, None]
        tail_rise = np.where(tail_mask, rise, -np.inf).max(axis=1)
        pileup = tail_rise > np.maximum(PILEUP_RISE_FRACTION * main_rise,
                                        PILEUP_RISE_MIN_ADC)

        quality = ~saturated & good_baseline & ~pileup & (amplitude > 0)
        stats["saturated"] += int(saturated.sum())
        stats["bad_baseline"] += int((~good_baseline & ~saturated).sum())
        stats["pileup"] += int((pileup & ~saturated & good_baseline).sum())
        with np.errstate(divide="ignore", invalid="ignore"):
            log_charge = np.log10(np.where(charge > 0, charge, np.nan))
        counts, _, _ = np.histogram2d(
            log_charge[quality], fprompt[quality],
            bins=(LOG_CHARGE_BINS[0], FPROMPT_BINS[0]),
            range=(LOG_CHARGE_BINS[1:], FPROMPT_BINS[1:]))
        hist2d += counts

        candidate = (quality & (fprompt > ER_BAND_BINS[1])
                     & (fprompt < ER_BAND_BINS[2])
                     & (charge > selection["charge_min"]))
        index = np.nonzero(candidate)[0]
        if index.size == 0:
            continue
        source = columns[None, :] + (t50[index] - ALIGN_SAMPLE)[:, None]
        valid = (source >= 0) & (source < n_samples)
        source = np.clip(source, 0, n_samples - 1)
        shape = (-1, n_bins, REBIN)
        rvalid = valid[:, :n_bins * REBIN].reshape(shape).all(axis=2)
        rebinned_channels = []
        for signal in signals:
            rows = signal[index]
            cumulative_rows = np.zeros((rows.shape[0], n_samples + 1))
            np.cumsum(rows, axis=1, out=cumulative_rows[:, 1:])
            norm = window_sum(cumulative_rows, first[index], n_samples) * DT_NS
            norm = np.where(norm > 0, norm, np.nan)
            aligned = np.take_along_axis(rows, source, axis=1) / norm[:, None]
            aligned = np.nan_to_num(np.where(valid, aligned, 0.0))
            rebinned_channels.append(
                aligned[:, :n_bins * REBIN].reshape(shape).mean(axis=2))
        cand_wvf.append(np.stack(rebinned_channels, axis=1).astype(np.float32))
        cand_valid.append(rvalid)
        cand_fprompt.append(fprompt[index])
        cand_charge.append(charge[index])

    acquired = datetime.datetime.fromtimestamp(
        (run_dir / CHANNEL_FILES[0]).stat().st_mtime)
    result = {"run": run_dir.name, "acquired": acquired, "stats": stats,
              "hist2d": hist2d, "band": None,
              "acc": {name: new_accumulator(n_bins) for name in SELECTION_VARIANTS}}
    if not cand_wvf:
        return result
    wvf = np.concatenate(cand_wvf)
    rvalid = np.concatenate(cand_valid)[:, None, :]
    fprompt = np.concatenate(cand_fprompt)
    charge = np.concatenate(cand_charge)
    stats["candidates"] = fprompt.size
    band = find_er_band(fprompt)
    if band is None:
        return result
    centre, sigma = band
    result["band"] = band
    for name, variant in SELECTION_VARIANTS.items():
        low_k, high_k = variant["fprompt_sigmas"]
        mask = ((fprompt >= centre - low_k * sigma) & (fprompt < centre + high_k * sigma)
                & (charge > selection["charge_min"] * variant.get("charge_factor", 1.0)))
        index = np.nonzero(mask)[0]
        for first_event in range(0, index.size, CHUNK_EVENTS):
            chunk = index[first_event:first_event + CHUNK_EVENTS]
            accumulate(result["acc"][name], wvf[chunk].astype(np.float64),
                       np.broadcast_to(rvalid[chunk], wvf[chunk].shape))
    return result


def merge(results, name="nominal"):
    merged = None
    for result in results:
        acc = result["acc"][name]
        if acc["n"] == 0:
            continue
        if merged is None:
            merged = {key: (value.copy() if isinstance(value, np.ndarray) else value)
                      for key, value in acc.items()}
            continue
        merged["n"] += acc["n"]
        for key in ("rsum", "rsum2", "rcount"):
            merged[key] += acc[key]
    return merged


def average(acc):
    """Rebinned mean waveform and its error (standard error of the mean)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = acc["rsum"] / acc["rcount"]
        variance = acc["rsum2"] / acc["rcount"] - mean ** 2
        error = np.sqrt(np.maximum(variance, 0.0) / acc["rcount"])
    return mean, error


def bin_times(n_bins):
    """Centre of the 40 ns bins, in ns from the t50."""
    return (np.arange(n_bins) + 0.5) * REBIN * DT_NS - ALIGN_SAMPLE * DT_NS


# --------------------------------------------------------------------- fit --

def fit_tail(mean, error, start=TAIL_FIT_START, end=TAIL_FIT_END, tau_spe=TAU_SPE):
    """Fit A_T exp(-t/tau_T) [+ A_S exp(-t/tau_spe)] + c to the average tail."""
    times = bin_times(mean.size)
    mask = (times >= start) & (times <= end) & (error > 0) & np.isfinite(mean)
    if mask.sum() < 20:
        return None
    x, y, e = times[mask], mean[mask], error[mask]

    def model(p):
        value = p[0] * np.exp(-x / p[1]) + p[2]
        if tau_spe:
            value = value + p[3] * np.exp(-x / tau_spe)
        return value

    p0 = [max(y[0], 1e-9) * math.exp(start / 1250.0), 1250.0, 0.0]
    lower, upper = [0.0, 200.0, -np.inf], [np.inf, 5000.0, np.inf]
    if tau_spe:
        p0.append(0.0)
        lower.append(-np.inf)
        upper.append(np.inf)
    result = least_squares(lambda p: (model(p) - y) / e, p0,
                           bounds=(lower, upper), x_scale="jac")
    chi2 = float((result.fun ** 2).sum())
    ndf = int(mask.sum() - len(p0))
    try:
        covariance = np.linalg.inv(result.jac.T @ result.jac)
    except np.linalg.LinAlgError:
        return None
    # Scale the errors so that the reduced chi2 is ~1 (never scale down).
    scale = math.sqrt(chi2 / ndf) if ndf > 0 and chi2 > ndf else 1.0
    errors = np.sqrt(np.diag(covariance)) * scale
    return {"params": result.x, "errors": errors, "chi2": chi2, "ndf": ndf,
            "tau": result.x[1], "tau_error": errors[1], "start": start,
            "end": end, "tau_spe": tau_spe, "error_scale": scale}


def fit_with_systematics(accumulators, channel_index, min_events=50):
    """Nominal fit plus systematics from the fit variants and the selection
    variants (added in quadrature)."""
    mean, error = average(accumulators["nominal"])
    nominal = fit_tail(mean[channel_index], error[channel_index])
    if nominal is None:
        return None
    fit_deviations = []
    for variant in SYSTEMATIC_VARIANTS:
        options = {"start": TAIL_FIT_START, "tau_spe": TAU_SPE, **variant}
        alternative = fit_tail(mean[channel_index], error[channel_index],
                               start=options["start"], tau_spe=options["tau_spe"])
        if alternative is not None:
            fit_deviations.append(abs(alternative["tau"] - nominal["tau"]))
    selection_deviations = []
    for name, acc in accumulators.items():
        if name == "nominal" or acc is None or acc["n"] < min_events:
            continue
        variant_mean, variant_error = average(acc)
        alternative = fit_tail(variant_mean[channel_index], variant_error[channel_index])
        if alternative is not None:
            selection_deviations.append(abs(alternative["tau"] - nominal["tau"]))
    nominal["tau_syst_fit"] = max(fit_deviations, default=0.0)
    nominal["tau_syst_sel"] = max(selection_deviations, default=0.0)
    nominal["tau_syst"] = math.hypot(nominal["tau_syst_fit"], nominal["tau_syst_sel"])
    return nominal


def contamination(tau, tau_error, tau0, k_q):
    """Equivalent impurity concentration (ppm) from the triplet lifetime."""
    value = (1.0 / tau - 1.0 / tau0) / k_q
    return value, tau_error / tau ** 2 / k_q


# ----------------------------------------------------------------- outputs --

def setup_root():
    import ROOT
    import tdrstyle

    tdrstyle.setTDRStyle()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    return ROOT


def to_histogram(ROOT, name, mean, error):
    n_bins = mean.size
    t0 = -ALIGN_SAMPLE * DT_NS
    histogram = ROOT.TH1D(name, ";Time from t_{50} [ns];Normalized amplitude [1/ns]",
                          n_bins, t0, t0 + n_bins * REBIN * DT_NS)
    for index in range(n_bins):
        if np.isfinite(mean[index]):
            histogram.SetBinContent(index + 1, float(mean[index]))
            histogram.SetBinError(index + 1, float(error[index]))
    histogram.SetDirectory(0)
    return histogram


def draw_waveform_fit(ROOT, label, mean, error, fits, path, n_events):
    canvas = ROOT.TCanvas(f"c_wvf_{label}", "", 800, 600)
    canvas.SetLogy()
    canvas.SetRightMargin(0.05)
    colors = {"ch0": ROOT.kRed, "ch1": ROOT.kBlue, "sum": ROOT.kBlack}
    histograms = {}
    positive = mean[2][np.isfinite(mean[2]) & (mean[2] > 0)]
    frame = canvas.DrawFrame(-1000.0, max(positive.min(), 1e-8) * 0.5,
                             16000.0, positive.max() * 5.0)
    frame.SetTitle(";Time from t_{50} [ns];Normalized amplitude [1/ns]")
    functions = []
    for k, channel in enumerate(CHANNELS):
        histogram = to_histogram(ROOT, f"h_draw_{label}_{channel}",
                                 mean[k], error[k])
        histogram.SetLineColor(colors[channel])
        histogram.SetMarkerColor(colors[channel])
        histogram.SetMarkerStyle(20)
        histogram.SetMarkerSize(0.4)
        histogram.Draw("P SAME")
        histograms[channel] = histogram
        fit = fits.get(channel)
        if fit is None:
            continue
        p = fit["params"]
        formula = "[0]*exp(-x/[1])+[2]"
        if fit["tau_spe"]:
            formula += f"+[3]*exp(-x/{fit['tau_spe']})"
        function = ROOT.TF1(f"f_{label}_{channel}", formula,
                            fit["start"], fit["end"])
        for index, value in enumerate(p):
            function.SetParameter(index, float(value))
        function.SetLineColor(colors[channel] if channel != "sum" else ROOT.kGreen + 2)
        function.SetLineWidth(2)
        function.Draw("SAME")
        functions.append(function)

    legend = ROOT.TLegend(0.36, 0.70, 0.93, 0.92)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.030)
    names = {"ch0": "top", "ch1": "bottom", "sum": "sum"}
    for channel in CHANNELS:
        fit = fits.get(channel)
        text = names[channel]
        if fit is not None:
            text += (f": #tau_{{T}} = {fit['tau']:.0f} #pm {fit['tau_error']:.0f} (stat)"
                     f" #pm {fit['tau_syst']:.0f} (syst) ns")
        legend.AddEntry(histograms[channel], text, "p")
    legend.Draw()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.035)
    latex.SetTextSize(0.028)
    latex.DrawLatex(0.17, 0.20, label)
    latex.DrawLatex(0.17, 0.16, f"N_{{events}} = {n_events}, fit from "
                                f"{TAIL_FIT_START / 1000:.1f} #mus")
    canvas.SaveAs(str(path))


def draw_selection(ROOT, label, hist2d, selection, band, path):
    histogram = ROOT.TH2D(f"h_sel_{label}",
                          ";log_{10}(charge [ADC#upointns]);fprompt",
                          LOG_CHARGE_BINS[0], *LOG_CHARGE_BINS[1:],
                          FPROMPT_BINS[0], *FPROMPT_BINS[1:])
    for i in range(LOG_CHARGE_BINS[0]):
        for j in range(FPROMPT_BINS[0]):
            histogram.SetBinContent(i + 1, j + 1, hist2d[i, j])
    canvas = ROOT.TCanvas(f"c_sel_{label}", "", 800, 600)
    canvas.SetLogz()
    canvas.SetRightMargin(0.15)
    histogram.Draw("COLZ")
    low_k, high_k = SELECTION_VARIANTS["nominal"]["fprompt_sigmas"]
    box = ROOT.TBox(math.log10(selection["charge_min"]), band[0] - low_k * band[1],
                    LOG_CHARGE_BINS[2], band[0] + high_k * band[1])
    box.SetFillStyle(0)
    box.SetLineColor(ROOT.kRed)
    box.SetLineWidth(2)
    box.Draw()
    canvas.SaveAs(str(path))


def draw_tau_vs_run(ROOT, label, rows, path):
    canvas = ROOT.TCanvas(f"c_runs_{label}", "", 800, 600)
    graphs = {}
    colors = {"ch0": ROOT.kRed, "ch1": ROOT.kBlue, "sum": ROOT.kBlack}
    values = []
    for channel in CHANNELS:
        graph = ROOT.TGraphErrors()
        for index, row in enumerate(rows):
            fit = row["fits"].get(channel)
            if fit is None:
                continue
            offset = {"ch0": -0.15, "ch1": 0.15, "sum": 0.0}[channel]
            point = graph.GetN()
            graph.SetPoint(point, index + 1 + offset, fit["tau"])
            graph.SetPointError(point, 0.0, math.hypot(fit["tau_error"], fit["tau_syst"]))
            values.append(fit["tau"])
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(colors[channel])
        graph.SetLineColor(colors[channel])
        graphs[channel] = graph
    if not values:
        return
    frame = canvas.DrawFrame(0.0, min(values) - 100.0, len(rows) + 1.0,
                             max(values) + 150.0)
    frame.SetTitle(";Run;#tau_{T} [ns]")
    axis = frame.GetXaxis()
    axis.SetNdivisions(len(rows) + 1, False)
    ROOT.gPad.SetGridy()
    for graph in graphs.values():
        if graph.GetN():
            graph.Draw("P SAME")
    labels = ROOT.TLatex()
    labels.SetTextFont(42)
    labels.SetTextSize(0.025)
    labels.SetTextAngle(90)
    for index, row in enumerate(rows):
        labels.DrawLatex(index + 1.1, min(values) - 90.0, row["run"])
    legend = ROOT.TLegend(0.70, 0.75, 0.92, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    for channel, text in zip(CHANNELS, ("top", "bottom", "sum")):
        legend.AddEntry(graphs[channel], text, "p")
    legend.Draw()
    canvas.SaveAs(str(path))


def draw_overview(ROOT, summary, config, path):
    """tau_T of every dataset (all runs merged) versus acquisition date."""
    canvas = ROOT.TCanvas("c_overview", "", 1000, 600)
    canvas.SetRightMargin(0.12)
    graph = ROOT.TGraphErrors()
    values = []
    for row in summary:
        fit = row["fits"].get("sum")
        if fit is None:
            continue
        timestamp = row["acquired"].timestamp()
        point = graph.GetN()
        graph.SetPoint(point, timestamp, fit["tau"])
        graph.SetPointError(point, 0.0, math.hypot(fit["tau_error"], fit["tau_syst"]))
        values.append((timestamp, fit["tau"]))
    if not values:
        return
    times = [t for t, _ in values]
    taus = [v for _, v in values]
    margin = 5 * 86400.0
    frame = canvas.DrawFrame(min(times) - margin, min(taus) - 100.0,
                             max(times) + margin, max(taus) + 100.0)
    frame.SetTitle(";Acquisition date;#tau_{T} (sum, all runs) [ns]")
    frame.GetXaxis().SetTimeDisplay(1)
    frame.GetXaxis().SetTimeFormat("%d/%m/%y%F1970-01-01 00:00:00")
    frame.GetXaxis().SetNdivisions(508)
    ROOT.gPad.SetGridy()
    graph.SetMarkerStyle(20)
    graph.Draw("P SAME")

    # Right axis: equivalent N2 contamination.
    tau0 = config["tau0_ns"]
    k_n2 = config["k_n2_per_ns_ppm"]
    low, high = frame.GetMinimum(), frame.GetMaximum()
    ppm_low = (1.0 / high - 1.0 / tau0) / k_n2
    ppm_high = (1.0 / low - 1.0 / tau0) / k_n2
    ROOT.gPad.Update()
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),
                       ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                       ppm_low, ppm_high, 510, "-L")
    axis.SetTitle("N_{2}-equivalent [ppm]")
    axis.SetTitleFont(42)
    axis.SetLabelFont(42)
    axis.SetTitleSize(0.045)
    axis.SetLabelSize(0.04)
    axis.SetTitleOffset(1.1)
    axis.Draw()
    note = ROOT.TLatex()
    note.SetNDC()
    note.SetTextFont(42)
    note.SetTextSize(0.03)
    note.DrawLatex(0.15, 0.92, f"right axis: #tau_{{T}}^{{0}} = {tau0:.0f} ns, "
                               f"k_{{N_{{2}}}} = {k_n2 * 1e3:.2f} #mus^{{-1}}ppm^{{-1}} "
                               "(non-linear scale)")
    canvas.SaveAs(str(path))


def setup_name(label):
    """Detector setup from the dataset label (no date, no bias voltage)."""
    name = re.sub(r"^\d{8}_", "", label)
    name = re.sub(r"_\d+_\d+V$", "", name)
    return name.replace("PDU_", "").replace("_", " ")


def utc_seconds(moment):
    """Naive local time as seconds for a ROOT time axis displayed in GMT."""
    return float(calendar.timegm(moment.timetuple()))


def split_periods(rows, gap_days=PERIOD_GAP_DAYS):
    periods = []
    for row in sorted(rows, key=lambda r: r["acquired"]):
        if periods and (row["acquired"] - periods[-1][-1]["acquired"]).total_seconds() \
                <= gap_days * 86400.0:
            periods[-1].append(row)
        else:
            periods.append([row])
    return periods


def draw_runs_vs_date(ROOT, run_rows, path):
    """tau_T (sum) of every single run versus acquisition time, one panel per
    data-taking period, colour by setup; filled markers for muon runs."""
    run_rows = [row for row in run_rows if "sum" in row["fits"]]
    if not run_rows:
        return
    periods = split_periods(run_rows)
    palette = (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kAzure + 7, ROOT.kRed + 1,
               ROOT.kMagenta + 1, ROOT.kOrange + 7, ROOT.kGray + 2, ROOT.kCyan + 2)
    setups = list(dict.fromkeys(setup_name(row["label"]) for row in
                                sorted(run_rows, key=lambda r: r["acquired"])))
    colors = {setup: palette[i % len(palette)] for i, setup in enumerate(setups)}
    taus = [row["fits"]["sum"]["tau"] for row in run_rows]
    y_low = 50.0 * math.floor((min(taus) - 80.0) / 50.0)
    y_high = 50.0 * math.ceil((max(taus) + 80.0) / 50.0)

    canvas = ROOT.TCanvas("c_runs_vs_date", "", 1600, 650)
    left, right, top = 0.07, 0.02, 0.80
    width = (1.0 - left - right) / len(periods)
    keep = []
    markers = {}
    for index, period in enumerate(periods):
        canvas.cd()
        x_low = 0.0 if index == 0 else left + index * width
        x_high = 1.0 if index == len(periods) - 1 else left + (index + 1) * width
        pad = ROOT.TPad(f"p_period_{index}", "", x_low, 0.0, x_high, top)
        pad_width = x_high - x_low
        pad.SetLeftMargin((left if index == 0 else 0.025) / pad_width)
        pad.SetRightMargin((right if index == len(periods) - 1 else 0.025) / pad_width)
        pad.SetBottomMargin(0.16)
        pad.SetTopMargin(0.07)
        pad.Draw()
        pad.cd()
        keep.append(pad)

        # Whole days (ticks at midnight, label dd/mm) or, for periods within
        # a single day, whole hours (label hh:mm).
        times = [utc_seconds(row["acquired"]) for row in period]
        step = 86400.0 if max(times) - min(times) > 6 * 3600.0 else 3600.0
        t_low = step * math.floor(min(times) / step)
        t_high = step * math.ceil(max(times) / step)
        if step == 86400.0 and t_high - t_low < 2 * step:
            t_high += step
        n_ticks = int(round((t_high - t_low) / step))
        frame = pad.DrawFrame(t_low, y_low, t_high, y_high)
        frame.SetTitle(";;#tau_{T} (sum) [ns]" if index == 0 else ";;")
        for axis in (frame.GetXaxis(), frame.GetYaxis()):
            axis.SetLabelFont(43)
            axis.SetLabelSize(17)
            axis.SetTitleFont(43)
            axis.SetTitleSize(22)
        frame.GetYaxis().SetTitleOffset(1.6)
        frame.GetXaxis().SetLabelOffset(0.015)
        if index > 0:
            frame.GetYaxis().SetLabelSize(0)
        frame.GetXaxis().SetTimeDisplay(1)
        frame.GetXaxis().SetTimeFormat(
            ("%d/%m" if step == 86400.0 else "%H:%M") + "%F1970-01-01 00:00:00s0 GMT")
        frame.GetXaxis().SetNdivisions(n_ticks + (400 if step == 86400.0 else 200), False)
        pad.SetGridy()
        pad.SetGridx()
        keep.append(frame)
        first, last = period[0]["acquired"], period[-1]["acquired"]
        header = (first.strftime("%d/%m/%Y") if first.date() == last.date() else
                  f"{first:%d/%m} - {last:%d/%m/%Y}")
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(43)
        text.SetTextSize(18)
        text.SetTextAlign(21)
        text.DrawLatex(pad.GetLeftMargin() + 0.5 * (1.0 - pad.GetLeftMargin()
                                                     - pad.GetRightMargin()), 0.945, header)
        keep.append(text)

        for row in period:
            fit = row["fits"]["sum"]
            setup = setup_name(row["label"])
            muon = row["run"] in MUON_RUNS
            graph = ROOT.TGraphErrors(1)
            graph.SetPoint(0, utc_seconds(row["acquired"]), fit["tau"])
            graph.SetPointError(0, 0.0, math.hypot(fit["tau_error"], fit["tau_syst"]))
            graph.SetMarkerStyle(20 if muon else 24)
            graph.SetMarkerSize(1.2 if muon else 0.9)
            graph.SetMarkerColor(colors[setup])
            graph.SetLineColor(colors[setup])
            graph.Draw("P SAME")
            keep.append(graph)
            markers.setdefault(setup, graph)
        pad.RedrawAxis()

    canvas.cd()
    legend = ROOT.TLegend(left, top + 0.01, 1.0 - right, 0.99)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(43)
    legend.SetTextSize(18)
    legend.SetNColumns(4)
    for setup in setups:
        marker = ROOT.TMarker(0.0, 0.0, 20)
        marker.SetMarkerColor(colors[setup])
        keep.append(marker)
        legend.AddEntry(marker, setup, "p")
    filled = ROOT.TMarker(0.0, 0.0, 20)
    hollow = ROOT.TMarker(0.0, 0.0, 24)
    keep.extend((filled, hollow))
    legend.AddEntry(filled, "muon run", "p")
    legend.AddEntry(hollow, "alpha run (ER events)", "p")
    legend.Draw()
    canvas.SaveAs(str(path))


# -------------------------------------------------------------------- main --

def analyze_dataset(dataset, config, plot_folder, workers):
    ROOT = setup_root()
    label = dataset["label"]
    selection = {**config["selection"], **dataset.get("selection", {})}
    runs = find_runs(dataset["folder"], dataset.get("exclude_runs", ()))
    if not runs:
        warnings.warn(f"[{label}] no runs with raw data")
        return None
    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(process_run, run, selection): run for run in runs}
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)
                st = result["stats"]
                band = result["band"]
                band_text = (f"fprompt band {band[0]:.3f} (sigma {band[1]:.3f})"
                             if band else "no band")
                print(f"[{label} {result['run']}] events {st['n_events']}, "
                      f"saturated {st['saturated']}, bad baseline "
                      f"{st['bad_baseline']}, pile-up {st['pileup']}, "
                      f"candidates {st['candidates']}, {band_text}, selected "
                      f"{result['acc']['nominal']['n']}", flush=True)
    order = {run.name: index for index, run in enumerate(runs)}
    results.sort(key=lambda r: order[r["run"]])
    used = [r for r in results if r["acc"]["nominal"]["n"] > 0]
    if not used:
        warnings.warn(f"[{label}] no electron-like events selected")
        return None

    out_dir = plot_folder / label
    out_dir.mkdir(parents=True, exist_ok=True)
    rootfile = ROOT.TFile.Open(str(out_dir / f"{label}_lar_purity.root"), "RECREATE")
    min_events = config["min_events_per_run"]
    rows = []
    for result in used:
        accumulators = result["acc"]
        mean, error = average(accumulators["nominal"])
        for k, channel in enumerate(CHANNELS):
            to_histogram(ROOT, f"h_{result['run']}_{channel}", mean[k], error[k]).Write()
        fits = {}
        if accumulators["nominal"]["n"] >= min_events:
            for k, channel in enumerate(CHANNELS):
                fit = fit_with_systematics(accumulators, k)
                if fit is not None:
                    fits[channel] = fit
        rows.append({"run": result["run"], "acquired": result["acquired"],
                     "n": accumulators["nominal"]["n"], "band": result["band"],
                     "fits": fits})

    merged = {name: merge(used, name) for name in SELECTION_VARIANTS}
    mean, error = average(merged["nominal"])
    merged_fits = {}
    for k, channel in enumerate(CHANNELS):
        to_histogram(ROOT, f"h_all_{channel}", mean[k], error[k]).Write()
        fit = fit_with_systematics(merged, k)
        if fit is not None:
            merged_fits[channel] = fit
    rootfile.Close()

    draw_waveform_fit(ROOT, label, mean, error, merged_fits,
                      out_dir / f"{label}_avg_waveform_fit.png", merged["nominal"]["n"])
    bands = np.array([r["band"] for r in used if r["band"] is not None])
    band = (float(np.median(bands[:, 0])), float(np.median(bands[:, 1])))
    draw_selection(ROOT, label, sum(r["hist2d"] for r in results), selection, band,
                   out_dir / f"{label}_selection.png")
    fitted_rows = [row for row in rows if row["fits"]]
    if fitted_rows:
        draw_tau_vs_run(ROOT, label, fitted_rows, out_dir / f"{label}_tau_vs_run.png")
    acquired = max(row["acquired"] for row in rows)
    summary_row = {"label": label, "run": "all", "acquired": acquired,
                   "n": merged["nominal"]["n"], "band": band, "fits": merged_fits,
                   "n_runs": len(used)}
    write_table(out_dir / f"{label}_tau_vs_run.csv", label, rows + [summary_row], config)
    for channel in CHANNELS:
        fit = merged_fits.get(channel)
        if fit:
            print(f"[{label}] {channel}: tau_T = {fit['tau']:.0f} +- "
                  f"{fit['tau_error']:.0f} (stat) +- {fit['tau_syst_fit']:.0f} (fit) +- "
                  f"{fit['tau_syst_sel']:.0f} (sel) ns, "
                  f"chi2/ndf = {fit['chi2']:.0f}/{fit['ndf']}")
    return summary_row, rows


TABLE_HEADER = [
    "label", "run", "acquired", "n_events", "fprompt_band", "fprompt_sigma",
    "channel", "tau_T_ns", "tau_T_stat_ns", "tau_T_syst_fit_ns", "tau_T_syst_sel_ns",
    "chi2", "ndf", "N2_eq_ppm", "N2_eq_err_ppm", "O2_eq_ppm", "O2_eq_err_ppm",
]


def table_rows(label, row, config):
    for channel in CHANNELS:
        fit = row["fits"].get(channel)
        if fit is None:
            continue
        total = math.hypot(fit["tau_error"], fit["tau_syst"])
        n2, n2_error = contamination(fit["tau"], total, config["tau0_ns"],
                                     config["k_n2_per_ns_ppm"])
        o2, o2_error = contamination(fit["tau"], total, config["tau0_ns"],
                                     config["k_o2_per_ns_ppm"])
        yield [
            label, row["run"], row["acquired"].strftime("%Y-%m-%d %H:%M"),
            row["n"], f"{row['band'][0]:.3f}", f"{row['band'][1]:.3f}",
            channel, f"{fit['tau']:.1f}", f"{fit['tau_error']:.1f}",
            f"{fit['tau_syst_fit']:.1f}", f"{fit['tau_syst_sel']:.1f}",
            f"{fit['chi2']:.1f}", fit["ndf"],
            f"{n2:.3f}", f"{n2_error:.3f}", f"{o2:.3f}", f"{o2_error:.3f}",
        ]


def read_table_rows(plot_folder, datasets, merged=True):
    """Rebuild the results from the per-dataset CSV tables: the "all" rows
    (merged=True) or the single runs (merged=False)."""
    rows = []
    for dataset in datasets:
        label = dataset["label"]
        path = plot_folder / label / f"{label}_tau_vs_run.csv"
        if not path.is_file():
            continue
        by_run = {}
        with open(path, newline="") as table:
            for entry in csv.DictReader(table):
                if (entry["run"] == "all") != merged:
                    continue
                row = by_run.get(entry["run"])
                if row is None:
                    row = by_run[entry["run"]] = {
                        "label": label, "run": entry["run"],
                        "acquired": datetime.datetime.strptime(
                            entry["acquired"], "%Y-%m-%d %H:%M"),
                        "n": int(entry["n_events"]),
                        "band": (float(entry["fprompt_band"]),
                                 float(entry["fprompt_sigma"])),
                        "fits": {}}
                fit_syst = float(entry["tau_T_syst_fit_ns"])
                sel_syst = float(entry["tau_T_syst_sel_ns"])
                row["fits"][entry["channel"]] = {
                    "tau": float(entry["tau_T_ns"]),
                    "tau_error": float(entry["tau_T_stat_ns"]),
                    "tau_syst_fit": fit_syst, "tau_syst_sel": sel_syst,
                    "tau_syst": math.hypot(fit_syst, sel_syst),
                    "chi2": float(entry["chi2"]), "ndf": int(entry["ndf"]),
                }
        rows.extend(by_run.values())
    return rows


def write_table(path, label, rows, config):
    with open(path, "w", newline="") as table:
        writer = csv.writer(table)
        writer.writerow(TABLE_HEADER)
        for row in rows:
            writer.writerows(table_rows(label, row, config))
    print(f"Written {path}")


def main():
    parser = argparse.ArgumentParser(description="LAr purity from the triplet lifetime")
    parser.add_argument("config", nargs="?", default=DEFAULT_CONFIG)
    parser.add_argument("--only", nargs="*", default=None,
                        help="labels of the datasets to process")
    parser.add_argument("--summary", action="store_true",
                        help="only rebuild the summary from the per-dataset tables")
    args = parser.parse_args()

    config = load_config(args.config)
    plot_folder = Path(config["plot_folder"]).expanduser()
    plot_folder.mkdir(parents=True, exist_ok=True)
    workers = config.get("workers", 4)

    if not args.summary:
        for dataset in config["datasets"]:
            if args.only and dataset["label"] not in args.only:
                continue
            analyze_dataset(dataset, config, plot_folder, workers)
        if args.only:
            return
    summary = read_table_rows(plot_folder, config["datasets"])
    ROOT = setup_root()
    summary.sort(key=lambda row: row["acquired"])
    with open(plot_folder / "lar_purity_summary.csv", "w", newline="") as table:
        writer = csv.writer(table)
        writer.writerow(TABLE_HEADER)
        for row in summary:
            writer.writerows(table_rows(row["label"], row, config))
    print(f"Written {plot_folder / 'lar_purity_summary.csv'}")
    draw_overview(ROOT, summary, config, plot_folder / "lar_purity_vs_date.png")
    draw_runs_vs_date(ROOT, read_table_rows(plot_folder, config["datasets"], merged=False),
                      plot_folder / "lar_purity_runs_vs_date.png")


if __name__ == "__main__":
    main()
