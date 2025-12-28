#!/bin/bash

DATA_DIR="./data/raw/planck"
mkdir -p $DATA_DIR

# Target the IRSA NASA mirror for NPIPE (PR4)
# Frequency 143 GHz is the gold standard for CMB polarization
BASE_URL="https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps"

echo "--- Fetching Planck NPIPE from NASA IRSA ---"

# Download the 143GHz Map (Includes I, Q, U)
curl -L -o $DATA_DIR/npipe_143.fits "$BASE_URL/HFI_SkyMap_143_2048_R3.01_full.fits"

# Download the Hits Map (The 'Context' of the scan)
curl -L -o $DATA_DIR/npipe_hits.fits "$BASE_URL/HFI_SkyMap_143_2048_R3.01_hits.fits"

echo "IRSA Planck download complete."
