#!/usr/bin/env python3
import numpy as np
import pandas as pd
import json
import os
import sys

def rotate_coords(ra, dec, theta_deg):
    theta_rad = np.radians(theta_deg)
    ra_c, dec_c = np.mean(ra), np.mean(dec)
    x, y = ra - ra_c, dec - dec_c
    x_rot = x * np.cos(theta_rad) - y * np.sin(theta_rad)
    y_rot = x * np.sin(theta_rad) + y * np.cos(theta_rad)
    return x_rot + ra_c, y_rot + dec_c

def get_variance_ratio(ra, dec, bins=10):
    ra_bins = np.linspace(np.min(ra), np.max(ra), bins + 1)
    dec_bins = np.linspace(np.min(dec), np.max(dec), bins + 1)
    hist, _, _ = np.histogram2d(ra, dec, bins=[ra_bins, dec_bins])
    inner_counts = hist[1:-1, 1:-1].flatten()
    mu = np.mean(inner_counts)
    return np.var(inner_counts, ddof=1) / mu if mu > 0 else 0

if __name__ == "__main__":
    path = "data/raw/COSMOS2020_subset.csv"
    
    if not os.path.exists(path):
        print(f"Error: {path} not found.")
        sys.exit(1)

    print(f"[*] REAL DATA AUDIT STARTING: Loading {path}...")
    df = pd.read_csv(path)
    
    # Standardize column names if they aren't already
    if 'ra' not in df.columns and 'ALPHA_J2000' in df.columns:
        df.rename(columns={'ALPHA_J2000': 'ra', 'DELTA_J2000': 'dec'}, inplace=True)
    
    # Filter for the Structuring Phase (3.0 < z < 6.0)
    subset = df[(df['lp_zPDF'] > 3.0) & (df['lp_zPDF'] < 6.0)]
    print(f"[*] Auditing {len(subset)} high-redshift sources...")
    
    results = []
    for theta in [0, 10, 20, 45, 90]:
        print(f"[*] Analyzing theta = {theta} degrees...")
        # Use the updated names here
        r_ra, r_dec = rotate_coords(subset['ra'], subset['dec'], theta)
        v = get_variance_ratio(r_ra, r_dec)
        results.append({"angle": theta, "v_ratio": v})
        print(f"    Measured Variance Ratio (V): {v:.4f}")

    v0, v20 = results[0]['v_ratio'], results[2]['v_ratio']
    a_score = (v0 - v20) / v0 if v0 > 0 else 0
    verdict = "PAROCHIAL" if a_score > 0.3 else "TERRITORY"
    
    final_report = {
        "anisotropy_score": a_score,
        "verdict": verdict,
        "data_points": results
    }

    os.makedirs("data/processed", exist_ok=True)
    with open("data/processed/p6_real_data_report.json", "w") as f:
        json.dump(final_report, f, indent=4)
        
    print(f"\n[FINAL VERDICT] Anisotropy Score: {a_score:.4f}")
    print(f"[*] Result: {verdict}")
