#!/usr/bin/env python3

import argparse
import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np
import ROOT
import tdrstyle


SCRIPT_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = (
    SCRIPT_DIR / "config_analysis_202604_PDU_ModuLYNO_ptp_blu.py"
)
INPUT_FILE_PATTERN = "analyzed_{window}ns.root"
TREE_NAME = "t1"


tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def load_config(config_file):
    """Load GAIN_VS_WINDOW_CONFIG from a Python configuration file."""
    config_path = Path(config_file).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    spec = importlib.util.spec_from_file_location("gain_vs_window_config", config_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load config file: {config_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "GAIN_VS_WINDOW_CONFIG"):
        raise AttributeError(
            f"{config_path} does not define GAIN_VS_WINDOW_CONFIG"
        )
    return module.GAIN_VS_WINDOW_CONFIG


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


def open_input_files(dataset):
    base_folder = Path(dataset["base_folder"]).expanduser()
    inputs = {}

    for window in dataset["windows"]:
        input_path = base_folder / INPUT_FILE_PATTERN.format(window=window)
        if not input_path.is_file():
            raise FileNotFoundError(f"Input file not found: {input_path}")

        root_file = ROOT.TFile.Open(str(input_path), "READ")
        if not root_file or root_file.IsZombie():
            raise OSError(f"Cannot open ROOT file: {input_path}")

        tree = root_file.Get(TREE_NAME)
        if not tree:
            root_file.Close()
            raise KeyError(f"Tree '{TREE_NAME}' not found in {input_path}")
        inputs[window] = (root_file, tree)

    return inputs


def analyze_channel(dataset, channel, channel_config, inputs, outfile, plot_folder, config):
    label = dataset["label"]
    windows = dataset["windows"]
    n_bins = channel_config["n_bins"]
    x_min = channel_config["x_min"]
    x_max = channel_config["x_max"]

    peak_search_sigma = dataset.get(
        "peak_search_sigma", config.get("peak_search_sigma", 10)
    )
    peak_search_threshold = dataset.get(
        "peak_search_threshold", config.get("peak_search_threshold", 0.01)
    )
    fit_half_width = dataset.get(
        "fit_half_width", config.get("fit_half_width", 0.25)
    )
    gain_peak_indices = dataset.get(
        "gain_peak_indices", config.get("gain_peak_indices", [1, 2])
    )
    snr_peak_index = dataset.get(
        "snr_peak_index", config.get("snr_peak_index", 1)
    )
    required_peaks = max(*gain_peak_indices, snr_peak_index) + 1

    histograms = {}
    fits = {}
    gains = {}
    gain_errors = {}
    snrs = {}

    v_ov = (
        dataset["v_bias"]
        - config["v_breakdown"]
        - config.get("series_voltage_drop", 0.0)
    )

    for window in windows:
        tree = inputs[window][1]
        histogram_name = f"h_{label}_{channel}_{window}ns"
        histogram = ROOT.TH1F(histogram_name, "", n_bins, x_min, x_max)
        histograms[window] = histogram
        tree.Draw(f"{channel}.charge >> {histogram_name}", "", "goff")

        spectrum = ROOT.TSpectrum(dataset.get("max_peaks", config.get("max_peaks", 5)))
        n_peaks = spectrum.Search(
            histogram, peak_search_sigma, "", peak_search_threshold
        )
        positions = spectrum.GetPositionX()
        positions = np.array([positions[i] for i in range(n_peaks)])
        positions.sort()
        print(f"[{label} {channel}, {window} ns] peaks: {positions}")

        if n_peaks < required_peaks:
            warnings.warn(
                f"{label} {channel}, {window} ns: found {n_peaks} peaks, "
                f"but at least {required_peaks} are required; skipping this window",
                stacklevel=2,
            )
            continue

        peak_spacing = positions[1] - positions[0]
        for peak_index, peak_position in enumerate(positions):
            fit_name = f"fit_{label}_{channel}_{window}ns_{peak_index}"
            fit_min = peak_position - fit_half_width * peak_spacing
            fit_max = peak_position + fit_half_width * peak_spacing
            fit = ROOT.TF1(fit_name, "gaus", fit_min, fit_max)
            fits[(window, peak_index)] = fit
            histogram.Fit(fit, "SRQ+", "", fit_min, fit_max)

        low_peak, high_peak = gain_peak_indices
        low_fit = fits[(window, low_peak)]
        high_fit = fits[(window, high_peak)]
        snr_fit = fits[(window, snr_peak_index)]
        gains[window] = high_fit.GetParameter(1) - low_fit.GetParameter(1)
        gain_errors[window] = math.hypot(
            high_fit.GetParError(1), low_fit.GetParError(1)
        )
        snrs[window] = snr_fit.GetParameter(1) / snr_fit.GetParameter(2)
        print(
            f"[{label} {channel}] Window: {window} ns, "
            f"Gain: {gains[window]:.2f} +/- {gain_errors[window]:.2f}, "
            f"SNR: {snrs[window]:.2f}"
        )

        canvas = ROOT.TCanvas(f"c_charge_{label}_{channel}_{window}ns", "", 800, 600)
        histogram.SetTitle(";Charge [ADC];Counts")
        histogram.Draw()
        outfile.cd()
        histogram.Write()
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.05)
        latex.SetTextFont(42)
        latex.DrawLatex(0.60, 0.84, f"{channel}, V_{{OV}} = {v_ov:.2f} V")
        canvas.SaveAs(str(plot_folder / f"{label}_{channel}_charge_{window:04d}ns.png"))

    valid_windows = [window for window in windows if window in gains]
    if not valid_windows:
        warnings.warn(
            f"{label} {channel}: no window has enough peaks; skipping this channel",
            stacklevel=2,
        )
        return

    canvas = ROOT.TCanvas(f"c_gain_snr_{label}_{channel}", "", 800, 600)
    plot_ranges = config["plot_ranges"] | dataset.get("plot_ranges", {})
    frame = ROOT.gPad.DrawFrame(
        plot_ranges.get("window_min", 0.0),
        plot_ranges.get("gain_min", 0.0),
        plot_ranges["window_max"],
        plot_ranges["gain_max"],
    )
    frame.SetTitle(";Integration window [ns];Gain [a.u.]")
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(0.90)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleColor(ROOT.kRed)
    frame.GetYaxis().SetLabelColor(ROOT.kRed)
    ROOT.gPad.SetGridy()
    ROOT.gPad.Update()

    gain_graph = ROOT.TGraphErrors(len(valid_windows))
    snr_graph = ROOT.TGraph(len(valid_windows))
    for index, window in enumerate(valid_windows):
        gain_graph.SetPoint(index, window, gains[window])
        gain_graph.SetPointError(index, 0, gain_errors[window])
        snr_graph.SetPoint(index, window, snrs[window])

    gain_graph.SetMarkerStyle(20)
    gain_graph.SetMarkerColor(ROOT.kRed)
    gain_graph.SetLineColor(ROOT.kRed)
    gain_graph.Draw("PL SAME")

    gain_min = frame.GetYaxis().GetXmin()
    gain_max = frame.GetYaxis().GetXmax()
    snr_min = plot_ranges.get("snr_min", 0.0)
    snr_max = plot_ranges["snr_max"]
    scale = (gain_max - gain_min) / (snr_max - snr_min)

    scaled_snr_graph = ROOT.TGraph(snr_graph.GetN())
    for index in range(snr_graph.GetN()):
        x_value = snr_graph.GetPointX(index)
        y_value = snr_graph.GetPointY(index)
        scaled_snr_graph.SetPoint(
            index, x_value, gain_min + (y_value - snr_min) * scale
        )

    scaled_snr_graph.SetMarkerStyle(20)
    scaled_snr_graph.SetMarkerColor(ROOT.kBlue)
    scaled_snr_graph.SetLineColor(ROOT.kBlue)
    scaled_snr_graph.Draw("PL SAME")
    ROOT.gPad.Update()

    axis = ROOT.TGaxis(
        ROOT.gPad.GetUxmax(), gain_min,
        ROOT.gPad.GetUxmax(), gain_max,
        snr_min, snr_max,
        510, "+L",
    )
    axis.SetLineColor(ROOT.kBlue)
    axis.SetLabelColor(ROOT.kBlue)
    axis.SetTitle("S / N")
    axis.SetTitleColor(ROOT.kBlue)
    axis.SetTitleSize(0.05)
    axis.SetTitleFont(42)
    axis.SetTitleOffset(0.80)
    axis.SetLabelSize(0.04)
    axis.SetLabelFont(42)
    axis.Draw()

    vov_label = ROOT.TLatex()
    vov_label.SetNDC()
    vov_label.SetTextSize(0.05)
    vov_label.SetTextFont(42)
    vov_label.DrawLatex(0.60, 0.84, f"{channel}, V_{{OV}} = {v_ov:.2f} V")

    canvas.SaveAs(str(plot_folder / f"{label}_{channel}_gain_and_SNR_vs_window.png"))

    outfile.cd()
    gain_graph.Write(f"gain_vs_window_{channel}")
    snr_graph.Write(f"SNR_vs_window_{channel}")


def analyze_dataset(dataset, config, plot_folder):
    label = dataset["label"]
    inputs = open_input_files(dataset)
    output_path = plot_folder / f"{label}_gain_vs_window.root"
    outfile = ROOT.TFile.Open(str(output_path), "RECREATE")
    if not outfile or outfile.IsZombie():
        for root_file, _ in inputs.values():
            root_file.Close()
        raise OSError(f"Cannot create ROOT file: {output_path}")

    try:
        for channel, channel_config in dataset["channels"].items():
            analyze_channel(
                dataset,
                channel,
                channel_config,
                inputs,
                outfile,
                plot_folder,
                config,
            )
    finally:
        outfile.Close()
        for root_file, _ in inputs.values():
            root_file.Close()


def main():
    parser = argparse.ArgumentParser(
        description="Study gain and SNR versus integration window"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG,
        help=f"Python config file (default: {DEFAULT_CONFIG.name})",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    plot_folder = prepare_plot_folder(config["plot_folder"])
    for dataset in config["datasets"]:
        analyze_dataset(dataset, config, plot_folder)


if __name__ == "__main__":
    main()
