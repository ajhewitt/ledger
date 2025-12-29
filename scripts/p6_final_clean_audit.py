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

def get_clean_v(ra, dec, bins=100):
    # We use a high-resolution square grid for the final "clean" baseline
    # to avoid the complex Moire patterns of low-res hexbins.
    hist, _, _ = np.histogram2d(ra, dec, bins=bins)
    counts = hist.flatten()
    # Filter out bins that fall outside the survey footprint
    counts = counts[counts > np.percentile(counts, 5)] 
    mu = np.mean(counts)
    return np.var(counts, ddof=1) / mu if mu > 0 else 0

if __name__ == "__main__":
    path = "data/raw/COSMOS2020_subset.csv"
    df = pd.read_csv(path)
    
    # Generate identical Random Null
    np.random.seed(42)
    ra_null = np.random.uniform(df['ra'].min(), df['ra'].max(), len(df))
    dec_null = np.random.uniform(df['dec'].min(), df['dec'].max(), len(df))

    angles = np.linspace(0, 90, 10)
    real_v, null_v = [], []

    print(f"[*] RUNNING FINAL CLEAN AUDIT (Resolution: 100x100)")
    for theta in angles:
        # Process Real
        r_ra, r_dec = rotate_coords(df['ra'], df['dec'], theta)
        real_v.append(get_clean_v(r_ra, r_dec))
        
        # Process Null
        n_ra, n_dec = rotate_coords(ra_null, dec_null, theta)
        null_v.append(get_clean_v(n_ra, n_dec))
        print(f"Angle {theta:2.0f}° | Real V: {real_v[-1]:.4f} | Null V: {null_v[-1]:.4f}")

    # Plotting the "Decoherence Verdict"
    plt.figure(figsize=(10, 6))
    plt.plot(angles, real_v, 'bo-', label='COSMOS2020 (Structuring Phase)')
    plt.plot(angles, null_v, 'r--', label='Poisson Null (Random)')
    plt.axhline(y=np.mean(null_v), color='gray', linestyle=':', label='Isotropic Floor')
    
    plt.title('P6: Rotational Invariance Audit (Cleaned)', fontsize=14)
    plt.xlabel('Rotation Angle (degrees)')
    plt.ylabel('Variance Ratio (V)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('data/processed/p6_final_clean_curve.png')
    print("[*] Final plot saved to data/processed/p6_final_clean_curve.png")
