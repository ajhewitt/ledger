#!/bin/bash

DATA_DIR="./data/raw/planck"
mkdir -p $DATA_DIR

# The NERSC NPIPE Portal path for Single-frequency
BASE_URL="https://portal.nersc.gov/project/cmb/planck2020/frequency_maps/Single-frequency"

echo "--- Fetching Verified Planck NPIPE (R4.00) 143 GHz ---"

# 1. Full Map (Signal + Noise + Hits)
curl -L -o $DATA_DIR/npipe_143_full.fits "$BASE_URL/HFI_SkyMap_143_2048_R4.00_full.fits"

# 2. Half-Ring Splits (Differential Audit)
curl -L -o $DATA_DIR/npipe_143_hr1.fits "$BASE_URL/HFI_SkyMap_143_2048_R4.00_full-ringhalf-1.fits"
curl -L -o $DATA_DIR/npipe_143_hr2.fits "$BASE_URL/HFI_SkyMap_143_2048_R4.00_full-ringhalf-2.fits"

echo "Download complete. Verify integrity: head -c 80 $DATA_DIR/npipe_143_full.fits"
