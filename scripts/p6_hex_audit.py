#!/usr/bin/env python3
import numpy as np
import pandas as pd
import json
import os

def rotate_coords(ra, dec, theta_deg):
    theta_rad = np.radians(theta_deg)
    ra_c, dec_c = np.mean(ra), np.mean(dec)
    x, y = ra - ra_c, dec - dec_c
    x_rot = x * np.cos(theta_rad) - y * np.sin(theta_rad)
    y_rot = x * np.sin(theta_rad) + y * np.cos(theta_rad)
    return x_rot + ra_c, y_rot + dec_c

def get_hex_variance(ra, dec, gridsize=10):
    """
    Bins data into hexagons and returns the Variance Ratio (V).
    """
    # Use numpy's fast histogram logic but adapt for a hexagonal 'shorthand' 
    # Or use a simplified approach: Hexagons can be thought of as offset rectangles.
    # For a high-rigor check, we'll use the counts from a hexbin collection.
    import matplotlib.pyplot as plt
    
    # We create a dummy figure to use the hexbin calculation engine
    fig, ax = plt.subplots()
    hb = ax.hexbin(ra, dec, gridsize=gridsize, mincnt=0)
    counts = hb.get_array()
    plt.close(fig) # Don't actually show the plot
    
    # Exclude empty space outside the survey area
    valid_counts = counts[counts > np.percentile(counts, 10)]
    
    mu = np.mean(valid_counts)
    return np.var(valid_counts, ddof=1) / mu if mu > 0 else 0

if __name__ == "__main__":
    path = "data/raw/COSMOS2020_subset.csv"
    if not os.path.exists(path):
        print("Error: Run extraction first.")
        exit(1)

    df = pd.read_csv(path)
    print(f"[*] HEXAGONAL AUDIT: Processing {len(df)} sources...")

    results = []
    # Test a wider range of angles to see if the "Square Spike" vanishes
    for theta in [0, 15, 30, 45, 60, 75, 90]:
        r_ra, r_dec = rotate_coords(df['ra'], df['dec'], theta)
        v = get_hex_variance(r_ra, r_dec)
        results.append({"angle": theta, "v_ratio": v})
        print(f"Angle {theta:02d}° | Hex-V: {v:.4f}")

    # Interpretation
    v_vals = [r['v_ratio'] for r in results]
    fluctuation = (np.max(v_vals) - np.min(v_vals)) / np.mean(v_vals)
    
    print(f"\n[*] Rotational Fluctuation: {fluctuation:.4f}")
    if fluctuation < 0.2:
        print("[VERDICT] TERRITORY: The structure is isotropic and real.")
    else:
        print("[VERDICT] CONTEXT: The structure is geometrically coupled.")
