#!/usr/bin/env python3

import argparse
import importlib.util
import math
import warnings
from pathlib import Path

import ROOT
import tdrstyle


SCRIPT_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = (
    SCRIPT_DIR / "config_analysis_202604_PDU_ModuLYNO_ptp_blu.py"
)
CHANNEL_STYLES = {
    "Ch0": {"marker": 20, "line_style": 1},
    "Ch1": {"marker": 22, "line_style": 2},
}


tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def load_config(config_file):
    config_path = Path(config_file).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    spec = importlib.util.spec_from_file_location("gain_config", config_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load config file: {config_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "GAIN_VS_WINDOW_CONFIG"):
        raise AttributeError(
            f"{config_path} does not define GAIN_VS_WINDOW_CONFIG"
        )
    if not hasattr(module, "GAIN_ANALYSIS_CONFIG"):
        raise AttributeError(f"{config_path} does not define GAIN_ANALYSIS_CONFIG")

    return module.GAIN_VS_WINDOW_CONFIG, module.GAIN_ANALYSIS_CONFIG


def find_point(graph, x_target):
    """Return the point index and y value at exactly x_target."""
    for index in range(graph.GetN()):
        x_value = graph.GetPointX(index)
        if math.isclose(x_value, x_target, rel_tol=0.0, abs_tol=1e-9):
            return index, graph.GetPointY(index)
    return None


def configured_channels(datasets):
    channels = []
    for dataset in datasets:
        for channel in dataset["channels"]:
            if channel not in channels:
                channels.append(channel)
    return channels


def select_datasets(window_config, gain_config):
    datasets_by_label = {
        dataset["label"]: dataset for dataset in window_config["datasets"]
    }
    selected = []
    for label in gain_config["datasets"]:
        dataset = datasets_by_label.get(label)
        if dataset is None:
            warnings.warn(
                f"Dataset '{label}' requested by GAIN_ANALYSIS_CONFIG is not "
                "defined in GAIN_VS_WINDOW_CONFIG; skipping it",
                stacklevel=2,
            )
            continue
        selected.append(dataset)
    return selected


def build_graphs(window_config, gain_config):
    plot_folder = Path(window_config["plot_folder"]).expanduser()
    datasets = select_datasets(window_config, gain_config)
    if not datasets:
        raise RuntimeError("No valid datasets selected for gain analysis")
    channels = configured_channels(datasets)
    gain_graphs = {channel: ROOT.TGraphErrors() for channel in channels}
    snr_graphs = {channel: ROOT.TGraphErrors() for channel in channels}
    input_files = []

    gain_window = gain_config["gain_window"]
    correction_window = gain_config["correction_window"]
    v_breakdown = window_config["v_breakdown"]
    voltage_drop = window_config.get("series_voltage_drop", 0.0)
    records = []
    correction_candidates = {channel: [] for channel in channels}

    for dataset in datasets:
        input_path = plot_folder / f'{dataset["label"]}_gain_vs_window.root'
        if not input_path.is_file():
            warnings.warn(f"Input file not found; skipping: {input_path}", stacklevel=2)
            continue

        input_file = ROOT.TFile.Open(str(input_path), "READ")
        if not input_file or input_file.IsZombie():
            warnings.warn(f"Cannot open ROOT file; skipping: {input_path}", stacklevel=2)
            continue
        input_files.append(input_file)

        v_ov = dataset["v_bias"] - v_breakdown - voltage_drop
        for channel in dataset["channels"]:
            gain_name = f"gain_vs_window_{channel}"
            snr_name = f"SNR_vs_window_{channel}"
            gain_vs_window = input_file.Get(gain_name)
            snr_vs_window = input_file.Get(snr_name)
            if not gain_vs_window or not snr_vs_window:
                warnings.warn(
                    f"{input_path}: missing {gain_name} or {snr_name}; "
                    f"skipping {channel}",
                    stacklevel=2,
                )
                continue

            gain_point = find_point(gain_vs_window, gain_window)
            correction_point = find_point(gain_vs_window, correction_window)
            snr_point = find_point(snr_vs_window, gain_window)
            if gain_point is None or snr_point is None:
                warnings.warn(
                    f"{dataset['label']} {channel}: required windows "
                    f"{gain_window} ns is not available for gain or SNR; "
                    "skipping this point",
                    stacklevel=2,
                )
                continue

            gain_index, gain = gain_point
            _, snr = snr_point
            if gain == 0 or not math.isfinite(gain):
                warnings.warn(
                    f"{dataset['label']} {channel}: gain is zero or invalid at "
                    f"{gain_window} ns; skipping this point",
                    stacklevel=2,
                )
                continue

            record = {
                "label": dataset["label"],
                "channel": channel,
                "v_ov": v_ov,
                "gain": gain,
                "gain_error": gain_vs_window.GetErrorY(gain_index),
                "snr": snr,
                "correction": None,
            }
            if correction_point is not None:
                _, correction_gain = correction_point
                correction = correction_gain / gain
                if math.isfinite(correction) and correction > 0:
                    record["correction"] = correction
                    correction_candidates[channel].append(
                        {
                            "label": dataset["label"],
                            "v_ov": v_ov,
                            "correction": correction,
                        }
                    )
                else:
                    warnings.warn(
                        f"{dataset['label']} {channel}: invalid correction "
                        f"{correction:.4g}; looking for a fallback",
                        stacklevel=2,
                    )
            records.append(record)

    for record in records:
        channel = record["channel"]
        correction = record["correction"]
        if correction is None:
            candidates = correction_candidates[channel]
            if not candidates:
                warnings.warn(
                    f"{record['label']} {channel}: {correction_window} ns is "
                    "not available and no fallback correction exists for this "
                    "channel; skipping this point",
                    stacklevel=2,
                )
                continue
            fallback = min(
                candidates,
                key=lambda candidate: abs(candidate["v_ov"] - record["v_ov"]),
            )
            correction = fallback["correction"]
            warnings.warn(
                f"{record['label']} {channel}: {correction_window} ns is not "
                f"available; using correction {correction:.4f} from "
                f"{fallback['label']}",
                stacklevel=2,
            )

        gain = record["gain"]
        gain_error = record["gain_error"]
        snr = record["snr"]
        v_ov = record["v_ov"]
        gain_graph = gain_graphs[channel]
        point_index = gain_graph.GetN()
        gain_graph.SetPoint(point_index, v_ov, gain * correction)
        gain_graph.SetPointError(point_index, 0.0, gain_error * correction)

        snr_graph = snr_graphs[channel]
        point_index = snr_graph.GetN()
        snr_graph.SetPoint(point_index, v_ov, snr)
        snr_graph.SetPointError(point_index, 0.0, 0.0)

        print(
            f"[{record['label']} {channel}] V_OV: {v_ov:.2f} V, "
            f"Gain: {gain * correction:.2f} +/- "
            f"{gain_error * correction:.2f}, SNR: {snr:.2f}"
        )

    return gain_graphs, snr_graphs, input_files


def style_graph(graph, channel, color):
    style = CHANNEL_STYLES.get(channel, {"marker": 20, "line_style": 1})
    graph.SetMarkerStyle(style["marker"])
    graph.SetMarkerSize(1.3)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineStyle(style["line_style"])
    graph.SetLineWidth(2)


def draw_results(gain_graphs, snr_graphs, gain_config, plot_folder):
    ranges = gain_config["plot_ranges"]
    canvas = ROOT.TCanvas("c_gain_SNR_vs_Vov", "c_gain_SNR_vs_Vov", 800, 700)
    frame = ROOT.gPad.DrawFrame(
        ranges["vov_min"],
        ranges["gain_min"],
        ranges["vov_max"],
        ranges["gain_max"],
    )
    frame.SetTitle(";V_{OV} [V];Gain [a.u.]")
    frame.GetYaxis().SetTitleColor(ROOT.kRed)
    frame.GetYaxis().SetLabelColor(ROOT.kRed)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicky(0)
    ROOT.gPad.Update()

    fits = {}
    equations = []
    drawable_channels = []
    for channel, graph in gain_graphs.items():
        if graph.GetN() == 0:
            warnings.warn(f"{channel}: no valid gain points to draw", stacklevel=2)
            continue

        style_graph(graph, channel, ROOT.kRed)
        graph.Draw("PL SAME")
        drawable_channels.append(channel)

        if graph.GetN() >= 2:
            style = CHANNEL_STYLES.get(channel, {"line_style": 1})
            fit = ROOT.TF1(
                f"fit_gain_{channel}",
                "[0]*(x+[1])",
                gain_config["fit_min"],
                gain_config["fit_max"],
            )
            fit.SetLineColor(ROOT.kRed)
            fit.SetLineWidth(2)
            fit.SetLineStyle(style["line_style"])
            graph.Fit(fit, "QNRS+")
            fit.Draw("SAME")
            fits[channel] = fit

            equation = ROOT.TLatex()
            equation.SetNDC()
            equation.SetTextSize(0.035)
            equation.SetTextColor(ROOT.kRed)
            equation.DrawLatex(
                0.15,
                0.84 - 0.05 * len(equations),
                f"{channel}: Gain = {fit.GetParameter(0):.1f} "
                f"#times (V_{{OV}} + {fit.GetParameter(1):.2f} V)",
            )
            equations.append(equation)

    gain_min = frame.GetYaxis().GetXmin()
    gain_max = frame.GetYaxis().GetXmax()
    snr_min = ranges["snr_min"]
    snr_max = ranges["snr_max"]
    scale = (gain_max - gain_min) / (snr_max - snr_min)

    scaled_snr_graphs = {}
    for channel, graph in snr_graphs.items():
        if graph.GetN() == 0:
            warnings.warn(f"{channel}: no valid SNR points to draw", stacklevel=2)
            continue

        scaled_graph = ROOT.TGraph(graph.GetN())
        for index in range(graph.GetN()):
            x_value = graph.GetPointX(index)
            y_value = graph.GetPointY(index)
            scaled_graph.SetPoint(
                index, x_value, gain_min + (y_value - snr_min) * scale
            )
        style_graph(scaled_graph, channel, ROOT.kBlue)
        scaled_graph.Draw("PL SAME")
        scaled_snr_graphs[channel] = scaled_graph

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
    axis.SetLabelSize(0.05)
    axis.SetLabelFont(42)
    axis.Draw()

    legend = ROOT.TLegend(0.62, 0.18, 0.86, 0.30)
    legend.SetBorderSize(0)
    legend.SetFillStyle(1001)
    legend.SetTextSize(0.045)
    legend.SetTextFont(42)
    for channel in drawable_channels:
        legend.AddEntry(gain_graphs[channel], channel, "PL")
    legend.Draw()

    canvas.SaveAs(str(plot_folder / "gain_and_SNR_vs_Vov.png"))
    return canvas, fits, scaled_snr_graphs, axis, legend, equations


def main():
    parser = argparse.ArgumentParser(
        description="Plot gain and SNR versus overvoltage for Ch0 and Ch1"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG,
        help=f"Python config file (default: {DEFAULT_CONFIG.name})",
    )
    args = parser.parse_args()

    window_config, gain_config = load_config(args.config)
    plot_folder = Path(window_config["plot_folder"]).expanduser()
    if not plot_folder.is_dir():
        raise NotADirectoryError(f"Plot folder does not exist: {plot_folder}")

    gain_graphs, snr_graphs, input_files = build_graphs(window_config, gain_config)
    output_path = plot_folder / "gain_vs_Vov.root"
    outfile = ROOT.TFile.Open(str(output_path), "RECREATE")
    if not outfile or outfile.IsZombie():
        for input_file in input_files:
            input_file.Close()
        raise OSError(f"Cannot create ROOT file: {output_path}")

    try:
        draw_results(gain_graphs, snr_graphs, gain_config, plot_folder)
        outfile.cd()
        for channel, graph in gain_graphs.items():
            graph.Write(f"gain_vs_Vov_{channel}")
        for channel, graph in snr_graphs.items():
            graph.Write(f"SNR_vs_Vov_{channel}")
    finally:
        outfile.Close()
        for input_file in input_files:
            input_file.Close()


if __name__ == "__main__":
    main()
