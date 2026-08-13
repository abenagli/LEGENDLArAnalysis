"""
ADC Data Processing Configuration

This file defines which datasets and runs to process, along with their parameters.
Modify the CONFIG dictionary below to match your data and processing needs.
"""

# Base path where all dataset folders are located
BASE_PATH = "/home/carla/Documents/ADC_data/Legend"

# Number of parallel workers for processing
N_WORKERS = 8

# Path to the ROOT macro (relative to script directory)
MACRO_FILE = "adc_read_all_data.C"

# Dataset configurations
CONFIG = {
    "datasets": [
        {
            "name": "20260616_PDU_MLYNO_verdenew_TPB_ThinPen",
            "runs": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            "charge_times": [6000],  # Integration times (ns)
            "startCharge": 3800,     # Start of charge window (ns)
        },
        {
            "name": "20260616_PDU_MLYNO_verdenew_TPB_ThinPen",
            "runs": [12, 13, 14],
            "charge_times": [400, 500, 750, 1000, 1500, 2000, 3000],
            "startCharge": 4000,
        },

        {
            "name": "20260617_PDU_MLYNO_verdenew_TPB_ThinPen",
            "runs": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            "charge_times": [6000],  # Integration times (ns)
            "startCharge": 3800,     # Start of charge window (ns)
        },
        {
            "name": "20260617_PDU_MLYNO_verdenew_TPB_ThinPen",
            "runs": [12, 13, 14],
            "charge_times": [400, 500, 750, 1000, 1500, 2000, 3000],
            "startCharge": 4000,
        },
    ]
}
