"""Configuration for the gain and alpha analyses."""

V_BREAKDOWN = 41.83
SERIES_VOLTAGE_DROP = 1.0
PLOT_FOLDER = "/var/www/html/LEGEND/202606_PDU_MLYNO_verdenew_PEN2"

GAIN_VS_WINDOW_CONFIG = {
    "plot_folder": PLOT_FOLDER,
    "v_breakdown": V_BREAKDOWN,
    "series_voltage_drop": SERIES_VOLTAGE_DROP,

    # Peak search and fit defaults. A dataset can override any of these values.
    "max_peaks": 5,
    "peak_search_sigma": 10,
    "peak_search_threshold": 0.01,
    "fit_half_width": 0.25,
    "gain_peak_indices": [1, 2],
    "snr_peak_index": 1,

    "plot_ranges": {
        "window_min": 0.0,
        "window_max": 5000.0,
        "gain_min": 0.0,
        "gain_max": 6000.0,
        "snr_min": 0.0,
        "snr_max": 15.0,
    },

    "datasets": [

        {
            "label": "20260611_45_7",
            "v_bias": 45.7,
            "base_folder": "/home/carla/Documents/ADC_data/Legend/20260611_PDU_MLYNO_verdenew_PEN2/run12",
            "windows": [400, 500, 750, 1000, 1500, 2000, 3000],
            "channels": {
                "Ch0": {
                    "n_bins": 1000,
                    "x_min": -5000.0,
                    "x_max": 20000.0,
                },
                "Ch1": {
                    "n_bins": 1000,
                    "x_min": -5000.0,
                    "x_max": 20000.0,
                },
            },
        },

    ],
}


# Configuration for study_gain.py
GAIN_ANALYSIS_CONFIG = {
    # Labels of the datasets to include in the V_OV plot.
    "datasets": [
        "20260611_45_7",
    ],
    # Gain and SNR are evaluated at this integration window.
    "gain_window": 1000,
    # Preserve the original correction to the gain measured at this window.
    "correction_window": 3000,
    "fit_min": 0.0,
    "fit_max": 5.0,
    "plot_ranges": {
        "vov_min": 0.0,
        "vov_max": 5.0,
        "gain_min": 0.0,
        "gain_max": 6000.0,
        "snr_min": 0.0,
        "snr_max": 20.0,
    },
}


# Configuration for study_alpha.py
ALPHA_ANALYSIS_CONFIG = {
    # Set to the output of study_gain.py to express results in photoelectrons.
    # Both gain_vs_Vov_Ch0 and gain_vs_Vov_Ch1 must be present.
    "gain_file": "/var/www/html/LEGEND/202606_PDU_MLYNO_verdenew_PEN2/gain_vs_Vov.root",
    "allow_gain_extrapolation": True,
    "charge_1pe_default": 1.0,
    "fprompt_cut": 0.48,

    "datasets": [
        {
            "label": "20260611_45_7",
            "v_bias": 45.7,
            "base_folder": "/home/carla/Documents/ADC_data/Legend/20260611_PDU_MLYNO_verdenew_PEN2/",
        },

        {
            "label": "20260612_45_7",
            "v_bias": 45.7,
            "base_folder": "/home/carla/Documents/ADC_data/Legend/20260612_PDU_MLYNO_verdenew_PEN2/",
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
