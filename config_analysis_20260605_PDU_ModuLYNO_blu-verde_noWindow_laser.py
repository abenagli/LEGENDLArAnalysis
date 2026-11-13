"""Configuration for the alpha analysis of the 2026-06-05 blu-verde module
(no window, laser). No gain measurement: results stay in ADC."""

V_BREAKDOWN = 42.09
SERIES_VOLTAGE_DROP = 1.0
PLOT_FOLDER = "/var/www/html/LEGEND/20260605_PDU_ModulYNO_blu-verde_noWindow_laser"


# Configuration for study_alpha.py
ALPHA_ANALYSIS_CONFIG = {
    "gain_file": None,
    "charge_1pe_default": 1.0,
    "fprompt_cut": 0.48,
    "input_filenames": ["analyzed.root"],
    # The 20, 24 and 28 cm positions were retaken as run41, run51, run61.
    "run_positions": {
        "run1": 4.0,
        "run2": 8.0,
        "run3": 12.0,
        "run41": 20.0,
        "run51": 24.0,
        "run61": 28.0,
        "run7": 36.0,
        "run8": 40.0,
        "run9": 44.0,
    },

    "datasets": [
        {
            "label": "20260605_55_0V",
            "v_bias": 55.0,
            "base_folder": "/home/carla/Documents/ADC_data/Legend/20260605_PDU_ModuLYNO_blu-verde_noWindow_laser_55_0V",
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
