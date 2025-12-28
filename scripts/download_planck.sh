#!/bin/bash
# scripts/download_planck.sh

# Ensure directory exists
mkdir -p data/raw

# Base URL for Planck Legacy Archive (ESA)
BASE_URL="http://pla.esac.esa.int/pla/aio/product-action?PRODUCT_ID"

echo "--- Starting Ledger Data Ingestion ---"

# 1. The Record: SMICA CMB Map (PR3 / 2018)
# This contains the I/Q/U signal we will audit.
echo "[1/3] Downloading Planck SMICA Map..."
wget -c -O data/raw/COM_CMB_IQU-smica_2048_R3.00_full.fits \
  "$BASE_URL=COM_CMB_IQU-smica_2048_R3.00_full.fits"

# 2. The Context Source A: Scan Strategy
# We download the 143 GHz frequency map. Field 3 (N_obs) contains 
# the hit-count map, which encodes the "cost" of the scan strategy.
echo "[2/3] Downloading HFI 143 GHz Map (for Scan History)..."
wget -c -O data/raw/HFI_SkyMap_143_2048_R3.01_full.fits \
  "$BASE_URL=HFI_SkyMap_143_2048_R3.01_full.fits"

# 3. The Context Source B: Zodiacal Light
# This model defines the local solar system foreground context.
echo "[3/3] Downloading Zodiacal Model..."
wget -c -O data/raw/COM_CompMap_Zodi-Model_2048_R2.00.fits \
  "$BASE_URL=COM_CompMap_Zodi-Model_2048_R2.00.fits"

echo "--- Ingestion Complete. Data stored in data/raw/ ---"
