#!/bin/bash
# scripts/download_p6_jwst.sh

DATA_DIR="./data/raw/jwst"
mkdir -p $DATA_DIR

# The exact IRSA NPIPE (PR4) path we found
BASE_URL="https://irsa.ipac.caltech.edu/data/COSMOS/tables/cosmos2020"

echo "--- Fetching COSMOS2020 Classic Ledger (3.9G) ---"

# This is the 2.2 release you saw in the directory listing
curl -L -o $DATA_DIR/cosmos2020_classic.fits.gz \
"$BASE_URL/COSMOS2020_CLASSIC_R1_v2.2_p3.fits.gz"

echo "JWST/COSMOS Ledger Entry acquired."
