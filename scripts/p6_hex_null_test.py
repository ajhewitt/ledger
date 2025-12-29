#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def rotate_coords(ra, dec, theta_deg):
    theta_rad = np.radians(theta_deg)
    ra_c, dec_c = np.mean(ra), np.mean(dec)
    x, y = ra - ra_c, dec - dec_c
    x_rot = x * np.cos(theta_rad) - y * np.sin(theta_rad)
    y_rot = x * np.sin(theta_rad) + y * np.cos(theta_rad)
    return x_rot + ra_c, y_rot + dec_c

def get_hex_variance(ra, dec, gridsize=10):
    fig, ax = plt.subplots()
    hb = ax.hexbin(ra, dec, gridsize=gridsize, mincnt=0)
    counts = hb.get_array()
    plt.close(fig)
    # Use only non-zero bins to represent the survey area
    valid_counts = counts[counts > 0]
    mu = np.mean(valid_counts)
    return np.var(valid_counts, ddof=1) / mu if mu > 0 else 0

if __name__ == "__main__":
    # Parameters from your real data run
    N_SOURCES = 205141
    RA_RANGE = [149.0, 151.0] # Approx COSMOS bounds
    DEC_RANGE = [1.5, 3.5]
    
    print(f"[*] RANDOM NULL TEST: Generating {N_SOURCES} Poisson sources...")
    np.random.seed(1337)
    ra_null = np.random.uniform(RA_RANGE[0], RA_RANGE[1], N_SOURCES)
    dec_null = np.random.uniform(DEC_RANGE[0], DEC_RANGE[1], N_SOURCES)

    angles = [0, 15, 30, 45, 60, 75, 90]
    print(f"{'Angle':<10} | {'Null Hex-V':<15}")
    print("-" * 30)

    null_results = []
    for theta in angles:
        r_ra, r_dec = rotate_coords(ra_null, dec_null, theta)
        v = get_hex_variance(r_ra, r_dec)
        null_results.append(v)
        print(f"{theta:02d}°        | {v:.4f}")

    avg_v = np.mean(null_results)
    fluctuation = (np.max(null_results) - np.min(null_results)) / avg_v
    
    print("-" * 30)
    print(f"[*] Null Rotational Fluctuation: {fluctuation:.4f}")
    
    if fluctuation < 0.1:
        print("[RESULT] The Hex-Audit is robust. Your real-data spike is likely ANOMALOUS.")
    else:
        print("[RESULT] Warning: Hex-binning shows geometric bias even in random data.")
