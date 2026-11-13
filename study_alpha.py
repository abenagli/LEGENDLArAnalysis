#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np
import ROOT
import tdrstyle


SCRIPT_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = SCRIPT_DIR / "config_analysis_202604_PDU_ModuLYNO_ptp_blu.py"

# Default geometry: run name -> source position in cm, measured from the
# bottom end. A config can override it with "run_positions" (globally or per
# dataset).
RUN_POSITIONS = {
    "run1": 4.0,
    "run2": 8.0,
    "run3": 12.0,
    "run4": 20.0,
    "run5": 24.0,
    "run6": 28.0,
    "run7": 36.0,
    "run8": 40.0,
    "run9": 44.0,
}

# Fixed peak-search and fit parameters.
MAX_PEAKS = 2
PEAK_SEARCH_SIGMA = 10
PEAK_SEARCH_THRESHOLD = 0.01
GAUSSIAN_SEED_HALF_WIDTH = 10000.0
GAUSSIAN_FIT_HALF_WIDTH = 100000.0
CRYSTALBALL_MIN_FACTOR = 0.75
CRYSTALBALL_MAX_FACTOR = 1.25
CRYSTALBALL_ALPHA = 1.0
CRYSTALBALL_N = 0.5
ALPHA_INPUT_FILENAMES = ("analyzed_6000ns.root", "analyzed.root")
DEFAULT_GAIN_GRAPHS = {"ch0": "gain_vs_Vov_Ch0", "ch1": "gain_vs_Vov_Ch1"}

# ln(S_bottom/S_top) fit: the first and last positions are excluded because
# they are affected by saturation from (quasi-)direct light.
RATIO_FIT_EXCLUDED_POINTS = 1


tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def load_config(config_file):
    config_path = Path(config_file).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    spec = importlib.util.spec_from_file_location("alpha_config", config_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load config file: {config_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "ALPHA_ANALYSIS_CONFIG"):
        raise AttributeError(f"{config_path} does not define ALPHA_ANALYSIS_CONFIG")
    common_names = ("V_BREAKDOWN", "SERIES_VOLTAGE_DROP", "PLOT_FOLDER")
    missing = [name for name in common_names if not hasattr(module, name)]
    if missing:
        raise AttributeError(
            f"{config_path} does not define the common parameter(s): "
            + ", ".join(missing)
        )
    common_config = {
        "v_breakdown": module.V_BREAKDOWN,
        "series_voltage_drop": module.SERIES_VOLTAGE_DROP,
        "plot_folder": module.PLOT_FOLDER,
    }
    return module.ALPHA_ANALYSIS_CONFIG, common_config


def prepare_plot_folder(folder):
    plot_folder = Path(folder).expanduser()
    if plot_folder.exists():
        if not plot_folder.is_dir():
            raise NotADirectoryError(f"Plot path is not a directory: {plot_folder}")
        warnings.warn(f"Plot folder already exists: {plot_folder}", stacklevel=2)
    else:
        plot_folder.mkdir(parents=True)
        print(f"Created plot folder: {plot_folder}")
    return plot_folder


def find_alpha_input(run_folder, filenames):
    for filename in filenames:
        input_path = run_folder / filename
        if input_path.is_file():
            return input_path
    return None


def load_channel_gains(dataset, config, v_ov):
    default_gain = config.get("charge_1pe_default", 1.0)
    gains = {"ch0": default_gain, "ch1": default_gain}
    gain_file_name = dataset.get("gain_file", config.get("gain_file"))
    if gain_file_name is None:
        print(
            f"[{dataset['label']}] No gain file provided. "
            "Keeping alpha results in ADC."
        )
        return gains, False

    gain_path = Path(gain_file_name).expanduser()
    if not gain_path.is_file():
        warnings.warn(
            f"[{dataset['label']}] Gain file not found: {gain_path}; "
            "keeping results in ADC",
            stacklevel=2,
        )
        return gains, False

    gain_file = ROOT.TFile.Open(str(gain_path), "READ")
    if not gain_file or gain_file.IsZombie():
        warnings.warn(
            f"[{dataset['label']}] Cannot open gain file: {gain_path}; "
            "keeping results in ADC",
            stacklevel=2,
        )
        return gains, False

    # Old gain files store a single graph (e.g. gain_vs_Vov_20260211) shared
    # by both channels; "gain_graphs" maps each channel to its graph name.
    gain_graphs = dataset.get(
        "gain_graphs", config.get("gain_graphs", DEFAULT_GAIN_GRAPHS)
    )
    values = {}
    try:
        for channel_key in ("ch0", "ch1"):
            graph_name = gain_graphs[channel_key]
            graph = gain_file.Get(graph_name)
            if not graph:
                warnings.warn(f"{graph_name} not found in {gain_path}", stacklevel=2)
                continue

            if not config.get("allow_gain_extrapolation", False):
                x_values = [graph.GetPointX(index) for index in range(graph.GetN())]
                if not x_values or v_ov < min(x_values) or v_ov > max(x_values):
                    warnings.warn(
                        f"[{dataset['label']}] V_OV = {v_ov:.2f} V is outside "
                        f"the range of {graph_name}; extrapolation is disabled",
                        stacklevel=2,
                    )
                    continue

            value = graph.Eval(v_ov)
            if not math.isfinite(value) or value <= 0:
                warnings.warn(
                    f"[{dataset['label']}] Invalid gain {value:.3f} from "
                    f"{graph_name}",
                    stacklevel=2,
                )
                continue
            values[channel_key] = value
    finally:
        gain_file.Close()

    if len(values) != 2:
        warnings.warn(
            f"[{dataset['label']}] Both channel gains are required; "
            "keeping all results in ADC",
            stacklevel=2,
        )
        return gains, False

    gains.update(values)
    print(
        f"[{dataset['label']}] Gain at V_OV = {v_ov:.2f} V: "
        f"Ch0 = {gains['ch0']:.3f} ADC, Ch1 = {gains['ch1']:.3f} ADC"
    )
    return gains, True


def fit_alpha_peak(histogram, label, run, channel):
    spectrum = ROOT.TSpectrum(MAX_PEAKS)
    n_peaks = spectrum.Search(
        histogram, PEAK_SEARCH_SIGMA, "", PEAK_SEARCH_THRESHOLD
    )
    positions = spectrum.GetPositionX()
    positions = np.array([positions[index] for index in range(n_peaks)])
    positions.sort()
    print(f"[{label} {run} {channel}] peaks: {positions}")
    if n_peaks == 0:
        warnings.warn(
            f"[{label} {run} {channel}] No peak found; skipping fit",
            stacklevel=2,
        )
        return None

    peak_position = positions[-1]
    gaussian = ROOT.TF1(
        f"fit_gaus_{label}_{channel}_{run}",
        "gaus",
        peak_position - GAUSSIAN_SEED_HALF_WIDTH,
        peak_position + GAUSSIAN_SEED_HALF_WIDTH,
    )
    histogram.Fit(
        gaussian,
        "NSRQ",
        "",
        peak_position - GAUSSIAN_FIT_HALF_WIDTH,
        peak_position + GAUSSIAN_FIT_HALF_WIDTH,
    )

    fit_min = CRYSTALBALL_MIN_FACTOR * peak_position
    fit_max = CRYSTALBALL_MAX_FACTOR * peak_position
    crystal_ball = ROOT.TF1(
        f"fit_{label}_{channel}_{run}", "crystalball", fit_min, fit_max
    )
    crystal_ball.SetParameters(
        gaussian.GetParameter(0),
        gaussian.GetParameter(1),
        gaussian.GetParameter(2),
        CRYSTALBALL_ALPHA,
        CRYSTALBALL_N,
    )
    histogram.Fit(crystal_ball, "SRQ", "", fit_min, fit_max)
    return crystal_ball


def analyze_dataset(dataset, config, common_config, plot_folder):
    label = dataset["label"]
    base_folder = Path(dataset["base_folder"]).expanduser()
    v_ov = (
        dataset["v_bias"]
        - common_config["v_breakdown"]
        - common_config["series_voltage_drop"]
    )
    gains, use_gain = load_channel_gains(dataset, config, v_ov)
    y_axis_label = "#alpha peak [N_{ph.e.}]" if use_gain else "#alpha peak [ADC]"

    histogram_config = config["histogram"]
    plot_ranges = config["plot_ranges"]
    channel_labels = config["channel_labels"]
    fprompt_cut = config["fprompt_cut"]
    run_positions = dataset.get(
        "run_positions", config.get("run_positions", RUN_POSITIONS)
    )
    input_filenames = tuple(config.get("input_filenames", ALPHA_INPUT_FILENAMES))

    histograms = {}
    fits = {}
    peaks = {}
    peak_errors = {}
    input_files = []
    existing_runs = set()

    output_path = plot_folder / f"{label}_alpha.root"
    outfile = ROOT.TFile.Open(str(output_path), "RECREATE")
    if not outfile or outfile.IsZombie():
        raise OSError(f"Cannot create ROOT file: {output_path}")

    try:
        for run in run_positions:
            run_folder = base_folder / run
            if not run_folder.is_dir():
                continue
            existing_runs.add(run)

            input_path = find_alpha_input(run_folder, input_filenames)
            if input_path is None:
                expected = ", ".join(input_filenames)
                warnings.warn(
                    f"Input file not found in {run_folder} "
                    f"(expected one of: {expected}); skipping",
                    stacklevel=2,
                )
                continue

            input_file = ROOT.TFile.Open(str(input_path), "READ")
            if not input_file or input_file.IsZombie():
                warnings.warn(f"Cannot open ROOT file; skipping: {input_path}", stacklevel=2)
                continue
            input_files.append(input_file)
            tree = input_file.Get("t1")
            if not tree:
                warnings.warn(f"Tree t1 not found; skipping: {input_path}", stacklevel=2)
                continue

            definitions = {
                "ch0": ("Ch0.charge", f"Ch0.fprompt>{fprompt_cut:.6g}"),
                "ch1": ("Ch1.charge", f"Ch1.fprompt>{fprompt_cut:.6g}"),
                "sum": (
                    "(Ch0.charge + Ch1.charge)",
                    f"Ch0.fprompt>{fprompt_cut:.6g} && "
                    f"Ch1.fprompt>{fprompt_cut:.6g}",
                ),
            }

            for channel, (expression, cut) in definitions.items():
                histogram_name = f"h_{label}_{channel}_{run}"
                histogram = ROOT.TH1F(
                    histogram_name,
                    "",
                    histogram_config["n_bins"],
                    histogram_config["x_min"],
                    histogram_config["x_max"],
                )
                tree.Draw(f"{expression} >> {histogram_name}", cut, "goff")
                histogram.SetDirectory(0)
                histograms[(run, channel)] = histogram
                fit = fit_alpha_peak(histogram, label, run, channel)
                if fit is None:
                    continue
                fits[(run, channel)] = fit
                peaks[(run, channel)] = fit.GetParameter(1)
                peak_errors[(run, channel)] = fit.GetParError(1)

            outfile.cd()
            for channel in definitions:
                histogram = histograms[(run, channel)]
                histogram.Write()
                fit = fits.get((run, channel))
                if fit:
                    fit.Write()

            if any((run, channel) not in fits for channel in definitions):
                warnings.warn(
                    f"[{label} {run}] Incomplete fits; skipping charge plot",
                    stacklevel=2,
                )
                continue

            canvas = ROOT.TCanvas(f"c_charge_{label}_{run}", "", 800, 600)
            maximum = 1.2 * max(
                histograms[(run, channel)].GetMaximum() for channel in definitions
            )
            colors = {"ch0": ROOT.kRed, "ch1": ROOT.kBlue, "sum": ROOT.kBlack}
            for channel in definitions:
                histograms[(run, channel)].SetLineColor(colors[channel])
                histograms[(run, channel)].SetLineWidth(2)
                fits[(run, channel)].SetLineColor(colors[channel])
                fits[(run, channel)].SetLineWidth(2)

            histograms[(run, "sum")].SetMaximum(maximum)
            histograms[(run, "sum")].SetTitle(";Charge [ADC];Counts")
            histograms[(run, "sum")].Draw()
            histograms[(run, "ch0")].Draw("SAME")
            histograms[(run, "ch1")].Draw("SAME")
            for channel in ("ch0", "ch1", "sum"):
                fits[(run, channel)].Draw("SAME")

            legend = ROOT.TLegend(0.60, 0.62, 0.88, 0.80)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(42)
            legend.SetTextSize(0.04)
            for channel in ("ch0", "ch1", "sum"):
                legend.AddEntry(
                    histograms[(run, channel)], channel_labels[channel], "l"
                )
            legend.Draw()

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.05)
            latex.SetTextFont(42)
            latex.DrawLatex(0.65, 0.84, f"V_{{OV}} = {v_ov:.2f} V")
            canvas.SaveAs(str(plot_folder / f"{label}_alpha_charge_{run}.png"))

        canvas = ROOT.TCanvas(f"c_position_{label}", "", 800, 600)
        y_axis_max = (
            plot_ranges["alpha_peak_max_pe"]
            if use_gain
            else plot_ranges["alpha_peak_max_adc"]
        )
        frame = ROOT.gPad.DrawFrame(
            plot_ranges["position_min"],
            0.0,
            plot_ranges["position_max"],
            y_axis_max,
        )
        ROOT.gPad.SetGridy()
        frame.SetTitle(f";Position [cm];{y_axis_label}")

        graphs = {
            "ch0": ROOT.TGraphErrors(),
            "ch1": ROOT.TGraphErrors(),
            "sum": ROOT.TGraphErrors(),
        }
        rows = []
        for run, position in sorted(run_positions.items(), key=lambda item: item[1]):
            if run not in existing_runs:
                continue
            if any((run, channel) not in peaks for channel in graphs):
                warnings.warn(
                    f"[{label} {run}] Incomplete peak results; skipping point",
                    stacklevel=2,
                )
                continue

            ch0_peak = peaks[(run, "ch0")] / gains["ch0"]
            ch0_error = peak_errors[(run, "ch0")] / gains["ch0"]
            ch1_peak = peaks[(run, "ch1")] / gains["ch1"]
            ch1_error = peak_errors[(run, "ch1")] / gains["ch1"]
            if use_gain:
                sum_peak = ch0_peak + ch1_peak
                sum_error = math.hypot(ch0_error, ch1_error)
            else:
                sum_peak = peaks[(run, "sum")]
                sum_error = peak_errors[(run, "sum")]

            values = {
                "ch0": (ch0_peak, ch0_error),
                "ch1": (ch1_peak, ch1_error),
                "sum": (sum_peak, sum_error),
            }
            for channel, (value, error) in values.items():
                graph = graphs[channel]
                index = graph.GetN()
                graph.SetPoint(index, position, value)
                graph.SetPointError(index, 0.0, error)
            rows.append({"run": run, "position": position, **values})

        if not any(graph.GetN() for graph in graphs.values()):
            warnings.warn(
                f"[{label}] No complete alpha-peak results; "
                "position plot will not be written",
                stacklevel=2,
            )
            return None

        styles = {
            "ch0": (ROOT.kRed, 20),
            "ch1": (ROOT.kBlue, 20),
            "sum": (ROOT.kBlack, 20),
        }
        for channel, graph in graphs.items():
            color, marker = styles[channel]
            graph.SetMarkerStyle(marker)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.Draw("PL SAME")

        legend = ROOT.TLegend(0.60, 0.62, 0.88, 0.80)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        for channel in ("ch0", "ch1", "sum"):
            legend.AddEntry(graphs[channel], channel_labels[channel], "pl")
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.05)
        latex.SetTextFont(42)
        latex.DrawLatex(0.65, 0.84, f"V_{{OV}} = {v_ov:.2f} V")
        canvas.SaveAs(str(plot_folder / f"{label}_alpha_peak_vs_position.png"))

        outfile.cd()
        for channel, graph in graphs.items():
            graph.Write(f"g_{channel}")

        unit = "pe" if use_gain else "ADC"
        write_position_table(
            plot_folder / f"{label}_alpha_peak_vs_position.csv", rows, unit
        )
        ratio_fit = draw_ratio_vs_position(
            label, rows, v_ov, plot_ranges, plot_folder, outfile
        )
        if ratio_fit is None:
            return None
        return {
            "label": label,
            "v_bias": dataset["v_bias"],
            "v_ov": v_ov,
            "unit": unit,
            **ratio_fit,
        }
    finally:
        outfile.Close()
        for input_file in input_files:
            input_file.Close()


def write_position_table(path, rows, unit):
    """Write S_top, S_bottom, S_sum and ln(S_bottom/S_top) versus position."""
    with open(path, "w", newline="") as table:
        writer = csv.writer(table)
        writer.writerow([
            "run", "position_cm",
            f"top_{unit}", f"top_err_{unit}",
            f"bottom_{unit}", f"bottom_err_{unit}",
            f"sum_{unit}", f"sum_err_{unit}",
            "ln_bottom_over_top", "ln_bottom_over_top_err",
        ])
        for row in rows:
            ln_ratio, ln_ratio_error = ln_bottom_over_top(row)
            writer.writerow([
                row["run"], f"{row['position']:g}",
                f"{row['ch0'][0]:.6g}", f"{row['ch0'][1]:.6g}",
                f"{row['ch1'][0]:.6g}", f"{row['ch1'][1]:.6g}",
                f"{row['sum'][0]:.6g}", f"{row['sum'][1]:.6g}",
                f"{ln_ratio:.6g}", f"{ln_ratio_error:.6g}",
            ])
    print(f"Written {path}")


def ln_bottom_over_top(row):
    top, top_error = row["ch0"]
    bottom, bottom_error = row["ch1"]
    if top <= 0 or bottom <= 0:
        return math.nan, math.nan
    return (
        math.log(bottom / top),
        math.hypot(bottom_error / bottom, top_error / top),
    )


def draw_ratio_vs_position(label, rows, v_ov, plot_ranges, plot_folder, outfile):
    """Plot ln(S_bottom/S_top) versus position and fit it with a straight line.

    With S_B ~ c_B exp(-x/lambda_eff) and S_T ~ c_T exp(-(L-x)/lambda_eff),
    ln(S_B/S_T) = ln(c_B/c_T) + L/lambda_eff - (2/lambda_eff) x, hence
    lambda_eff = -2/slope. The first and last points are excluded from the fit.
    """
    points = []
    for row in rows:
        ln_ratio, ln_ratio_error = ln_bottom_over_top(row)
        if math.isfinite(ln_ratio):
            points.append((row["position"], ln_ratio, ln_ratio_error))

    n_excluded = RATIO_FIT_EXCLUDED_POINTS
    fit_points = points[n_excluded:len(points) - n_excluded]
    excluded_points = points[:n_excluded] + points[len(points) - n_excluded:]
    if len(fit_points) < 2:
        warnings.warn(
            f"[{label}] Only {len(fit_points)} point(s) left after excluding "
            "the first and last positions; ratio plot will not be written",
            stacklevel=2,
        )
        return None
    if len(fit_points) == 2:
        warnings.warn(
            f"[{label}] Ratio fit with only 2 points (ndf = 0)", stacklevel=2
        )

    graph_fit = ROOT.TGraphErrors()
    for position, value, error in fit_points:
        index = graph_fit.GetN()
        graph_fit.SetPoint(index, position, value)
        graph_fit.SetPointError(index, 0.0, error)
    graph_excluded = ROOT.TGraphErrors()
    for position, value, error in excluded_points:
        index = graph_excluded.GetN()
        graph_excluded.SetPoint(index, position, value)
        graph_excluded.SetPointError(index, 0.0, error)

    fit_min = fit_points[0][0]
    fit_max = fit_points[-1][0]
    line = ROOT.TF1(f"fit_ratio_{label}", "pol1", fit_min, fit_max)
    graph_fit.Fit(line, "SRQ", "", fit_min, fit_max)
    # The peak-position errors are purely statistical and chi2/ndf >> 1: scale
    # the parameter errors so that the reduced chi2 is ~1 (never scale down).
    chi2 = line.GetChisquare()
    ndf = line.GetNDF()
    error_scale = math.sqrt(chi2 / ndf) if ndf > 0 and chi2 > ndf else 1.0
    intercept = line.GetParameter(0)
    intercept_error = line.GetParError(0) * error_scale
    slope = line.GetParameter(1)
    slope_error = line.GetParError(1) * error_scale
    lambda_eff = -2.0 / slope if slope != 0 else math.inf
    lambda_eff_error = 2.0 * slope_error / slope**2 if slope != 0 else math.inf
    print(
        f"[{label}] ln(S_B/S_T) fit: intercept = {intercept:.4f} +- "
        f"{intercept_error:.4f}, slope = {slope:.5f} +- {slope_error:.5f} cm^-1, "
        f"lambda_eff = {lambda_eff:.1f} +- {lambda_eff_error:.1f} cm "
        f"(errors scaled by {error_scale:.2f})"
    )

    y_values = [value for _, value, _ in points]
    y_errors = [error for _, _, error in points]
    y_low = min(value - error for value, error in zip(y_values, y_errors))
    y_high = max(value + error for value, error in zip(y_values, y_errors))
    y_span = max(y_high - y_low, 0.1)

    canvas = ROOT.TCanvas(f"c_ratio_{label}", "", 800, 600)
    frame = ROOT.gPad.DrawFrame(
        plot_ranges["position_min"],
        y_low - 0.15 * y_span,
        plot_ranges["position_max"],
        y_high + 0.9 * y_span,
    )
    ROOT.gPad.SetGridy()
    frame.SetTitle(";Position [cm];ln(S_{bottom}/S_{top})")

    graph_fit.SetMarkerStyle(20)
    graph_fit.SetMarkerColor(ROOT.kBlack)
    graph_fit.SetLineColor(ROOT.kBlack)
    graph_excluded.SetMarkerStyle(24)
    graph_excluded.SetMarkerColor(ROOT.kGray + 2)
    graph_excluded.SetLineColor(ROOT.kGray + 2)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    graph_fit.Draw("P SAME")
    if graph_excluded.GetN():
        graph_excluded.Draw("P SAME")
    line.Draw("SAME")

    legend = ROOT.TLegend(0.18, 0.18, 0.48, 0.32)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.AddEntry(graph_fit, "fitted points", "pe")
    if graph_excluded.GetN():
        legend.AddEntry(graph_excluded, "excluded (saturation)", "pe")
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.65, 0.84, f"V_{{OV}} = {v_ov:.2f} V")
    latex.SetTextSize(0.035)
    text_lines = [
        "fit: ln(S_{B}/S_{T}) = p_{0} + p_{1} x",
        f"p_{{0}} = {intercept:.3f} #pm {intercept_error:.3f}",
        f"p_{{1}} = ({1e3 * slope:.2f} #pm {1e3 * slope_error:.2f})"
        "#times10^{-3} cm^{-1}",
        f"#lambda_{{eff}} = -2/p_{{1}} = {lambda_eff:.1f} #pm "
        f"{lambda_eff_error:.1f} cm",
        f"#chi^{{2}}/ndf = {chi2:.1f}/{ndf}",
    ]
    if error_scale > 1.0:
        text_lines.append(
            f"errors scaled by #sqrt{{#chi^{{2}}/ndf}} = {error_scale:.2f}"
        )
    for index, text in enumerate(text_lines):
        latex.DrawLatex(0.50, 0.77 - 0.05 * index, text)
    canvas.SaveAs(str(plot_folder / f"{label}_alpha_ratio_vs_position.png"))

    outfile.cd()
    graph_fit.Write("g_ln_ratio_fit")
    graph_excluded.Write("g_ln_ratio_excluded")
    line.Write()

    return {
        "n_fit_points": len(fit_points),
        "intercept": intercept,
        "intercept_error": intercept_error,
        "slope": slope,
        "slope_error": slope_error,
        "lambda_eff": lambda_eff,
        "lambda_eff_error": lambda_eff_error,
        "chi2": chi2,
        "ndf": ndf,
        "error_scale": error_scale,
    }


def write_lambda_summary(path, results):
    with open(path, "w", newline="") as table:
        writer = csv.writer(table)
        writer.writerow([
            "label", "v_bias_V", "v_ov_V", "unit", "n_fit_points",
            "intercept", "intercept_err", "slope_per_cm", "slope_err_per_cm",
            "lambda_eff_cm", "lambda_eff_err_cm", "chi2", "ndf", "error_scale",
        ])
        for result in results:
            writer.writerow([
                result["label"], f"{result['v_bias']:g}",
                f"{result['v_ov']:.3f}", result["unit"],
                result["n_fit_points"],
                f"{result['intercept']:.6g}", f"{result['intercept_error']:.6g}",
                f"{result['slope']:.6g}", f"{result['slope_error']:.6g}",
                f"{result['lambda_eff']:.6g}", f"{result['lambda_eff_error']:.6g}",
                f"{result['chi2']:.6g}", result["ndf"],
                f"{result['error_scale']:.4g}",
            ])
    print(f"Written {path}")


def main():
    parser = argparse.ArgumentParser(description="Analyze alpha data versus position")
    parser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG,
        help=f"Python config file (default: {DEFAULT_CONFIG.name})",
    )
    args = parser.parse_args()

    config, common_config = load_config(args.config)
    plot_folder = prepare_plot_folder(common_config["plot_folder"])
    ratio_results = []
    for dataset in config["datasets"]:
        result = analyze_dataset(dataset, config, common_config, plot_folder)
        if result is not None:
            ratio_results.append(result)
    if ratio_results:
        write_lambda_summary(plot_folder / "lambda_eff_summary.csv", ratio_results)


if __name__ == "__main__":
    main()
