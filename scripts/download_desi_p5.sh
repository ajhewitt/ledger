#!/bin/bash

# Exact URL from verified index (12/28/2025)
BASE_URL="https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5"
DATA_DIR="./data/raw/desi"
mkdir -p $DATA_DIR

echo "--- Fetching DESI DR1 v1.5 LRG Catalogs (The Record) ---"

# 1. Download LRG Clustering Catalogs (NGC and SGC)
curl -L -o $DATA_DIR/LRG_NGC_clustering.dat.fits "$BASE_URL/LRG_NGC_clustering.dat.fits"
curl -L -o $DATA_DIR/LRG_SGC_clustering.dat.fits "$BASE_URL/LRG_SGC_clustering.dat.fits"

# 2. Download LRG Full Catalog (Contains context/target bitmasks)
curl -L -o $DATA_DIR/LRG_full.dat.fits "$BASE_URL/LRG_full.dat.fits"

# 3. Download a single Random catalog for NGC to enable basic clustering proxies
curl -L -o $DATA_DIR/LRG_NGC_0_clustering.ran.fits "$BASE_URL/LRG_NGC_0_clustering.ran.fits"

echo "Download complete. Files are in $DATA_DIR"
