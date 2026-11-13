"""Configuration for study_lar_purity.py (LAr purity from the triplet lifetime).

All runs with raw data (0_wave0.dat, 0_wave1.dat) of every dataset are
processed; runs without electron-like events (SPE, laser) are skipped.
"""

_DATA = "/home/carla/Documents/ADC_data/Legend"

_DATASETS = [
    "20260205_test",
    "20260210_PDU_cold_test_45_3V",
    "20260211_PDU_cold_test_42_1V",
    "20260211_PDU_cold_test_43_7V",
    "20260211_PDU_cold_test_44_7V",
    "20260211_PDU_cold_test_45_3V",
    "20260211_PDU_cold_test_45_7V",
    "20260211_PDU_cold_test_47_2V",
    "20260212_PDU_cold_test_45_3V",
    "20260212_PDU_cold_test_45_7V",
    "20260213_PDU_cold_test_42_2V",
    "20260213_PDU_cold_test_43_7V",
    "20260213_PDU_cold_test_44_7V",
    "20260213_PDU_cold_test_45_3V",
    "20260213_PDU_cold_test_45_7V",
    "20260213_PDU_cold_test_47_2V",
    "20260409_PDU_LNGS_PEN_cold_test_44_7V",
    "20260409_PDU_LNGS_PEN_cold_test_45_3V",
    "20260409_PDU_LNGS_PEN_cold_test_45_7V",
    "20260409_PDU_LNGS_PEN_cold_test_47_2V",
    "20260415_PDU_ModuLYNO_ptp_blu_44_7V",
    "20260416_PDU_ModuLYNO_ptp_blu_45_7V",
    "20260416_PDU_ModuLYNO_ptp_blu_45_83V",
    "20260416_PDU_ModuLYNO_ptp_blu_47_2V",
    "20260417_PDU_ModuLYNO_ptp_blu_44_7V",
    "20260417_PDU_ModuLYNO_ptp_blu_45_83V",
    "20260605_PDU_ModuLYNO_blu-verde_noWindow_laser_55_0V",
    "20260609_PDU_MLYNO_verdenew_PEN_TPB",
    "20260610_PDU_MLYNO_verdenew_PEN_TPB",
    "20260611_PDU_MLYNO_verdenew_PEN2",
    "20260612_PDU_MLYNO_verdenew_PEN2",
    "20260616_PDU_MLYNO_verdenew_TPB_ThinPen",
    "20260617_PDU_MLYNO_verdenew_TPB_ThinPen",
]

# 20260617: in the alpha runs (acquired after run11) the top channel (Ch0) has a
# distorted pulse shape (tau_T 1.6-2.8 us on Ch0 vs ~1.3 us on Ch1, alpha fprompt
# 0.46 instead of 0.57): only the muon run11 is used.
_EXCLUDE_RUNS = {
    "20260617_PDU_MLYNO_verdenew_TPB_ThinPen": (
        "run1", "run2", "run3", "run4", "run5", "run6", "run21", "run34", "run67"),
}

PURITY_CONFIG = {
    "plot_folder": "/var/www/html/LEGEND/LAr_purity",
    "workers": 4,
    # Runs with fewer selected events are merged into the dataset average but
    # not fitted on their own.
    "min_events_per_run": 200,

    # Electron-like selection, on the sum of the two channels. The charge
    # (ADC*ns, ~3400 ADC*ns per photoelectron) is integrated from 40 ns before
    # the t50 to the end of the record; fprompt = first 600 ns / total. The
    # fprompt window is found run by run around the electron-like peak.
    "selection": {
        "charge_min": 5.0e5,
        "baseline_rms_factor": 3.0,
    },

    # 1/tau_T = 1/tau_T0 + k_Q [X]. tau_T0 = 1.6 us for pure LAr; quenching
    # constants from WArP: N2 0.11 /us/ppm (JINST 5 P06003),
    # O2 0.54 /us/ppm (JINST 5 P05003). Units: ns and ppm.
    "tau0_ns": 1600.0,
    "k_n2_per_ns_ppm": 0.11e-3,
    "k_o2_per_ns_ppm": 0.54e-3,

    "datasets": [
        {"label": name, "folder": f"{_DATA}/{name}",
         "exclude_runs": _EXCLUDE_RUNS.get(name, ())}
        for name in _DATASETS
    ],
}
