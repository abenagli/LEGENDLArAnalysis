"""Configuration for the alpha analysis of the February 2026 PDU cold test.

Only ALPHA_ANALYSIS_CONFIG is provided: the gain files of this campaign
(<date>_gain.root) were produced by an older version of study_gain.py and
store one gain graph per day, shared by both channels.
"""

V_BREAKDOWN = 42.09  # V, extrapolated from LN2 to LAr with 44 mV/K
SERIES_VOLTAGE_DROP = 1.0
PLOT_FOLDER = "/var/www/html/LEGEND/202602_PDU_cold_test"

_DATA = "/home/carla/Documents/ADC_data/Legend"


# Configuration for study_alpha.py
ALPHA_ANALYSIS_CONFIG = {
    "allow_gain_extrapolation": True,
    "charge_1pe_default": 1.0,
    "fprompt_cut": 0.48,
    # Alpha runs of this campaign were reconstructed into analyzed.root.
    "input_filenames": ["analyzed.root"],

    "datasets": [
        {
            "label": "20260211_44_7V",
            "v_bias": 44.7,
            "base_folder": f"{_DATA}/20260211_PDU_cold_test_44_7V",
            "gain_file": f"{PLOT_FOLDER}/20260211_gain.root",
            "gain_graphs": {"ch0": "gain_vs_Vov_20260211", "ch1": "gain_vs_Vov_20260211"},
        },
        {
            "label": "20260211_45_3V",
            "v_bias": 45.3,
            "base_folder": f"{_DATA}/20260211_PDU_cold_test_45_3V",
            "gain_file": f"{PLOT_FOLDER}/20260211_gain.root",
            "gain_graphs": {"ch0": "gain_vs_Vov_20260211", "ch1": "gain_vs_Vov_20260211"},
        },
        {
            "label": "20260211_45_7V",
            "v_bias": 45.7,
            "base_folder": f"{_DATA}/20260211_PDU_cold_test_45_7V",
            "gain_file": f"{PLOT_FOLDER}/20260211_gain.root",
            "gain_graphs": {"ch0": "gain_vs_Vov_20260211", "ch1": "gain_vs_Vov_20260211"},
        },
        {
            "label": "20260213_45_3V",
            "v_bias": 45.3,
            "base_folder": f"{_DATA}/20260213_PDU_cold_test_45_3V",
            "gain_file": f"{PLOT_FOLDER}/20260213_gain.root",
            "gain_graphs": {"ch0": "gain_vs_Vov_20260213", "ch1": "gain_vs_Vov_20260213"},
        },
        {
            "label": "20260213_45_7V",
            "v_bias": 45.7,
            "base_folder": f"{_DATA}/20260213_PDU_cold_test_45_7V",
            "gain_file": f"{PLOT_FOLDER}/20260213_gain.root",
            "gain_graphs": {"ch0": "gain_vs_Vov_20260213", "ch1": "gain_vs_Vov_20260213"},
        },
    ],

    "histogram": {
        "n_bins": 500,
        "x_min": 0.0,
        "x_max": 3000000.0,
    },
    "plot_ranges": {
        "position_min": 0.0,
        "position_max": 50.0,
        "alpha_peak_max_adc": 1500000.0,
        "alpha_peak_max_pe": 1000.0,
    },
    "channel_labels": {
        "ch0": "top",
        "ch1": "bottom",
        "sum": "sum",
    },
}
