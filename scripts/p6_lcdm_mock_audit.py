#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter

def rotate_coords(ra, dec, theta_deg):
    theta_rad = np.radians(theta_deg)
    ra_c, dec_c = np.mean(ra), np.mean(dec)
    x, y = ra - ra_c, dec - dec_c
    x_rot = x * np.cos(theta_rad) - y * np.sin(theta_rad)
    y_rot = x * np.sin(theta_rad) + y * np.cos(theta_rad)
    return x_rot + ra_c, y_rot + dec_c

def get_clean_v(ra, dec, bins=100):
    hist, _, _ = np.histogram2d(ra, dec, bins=bins)
    counts = hist.flatten()
    # Filter for survey footprint (interior 90% of bins)
    counts = counts[counts > np.percentile(counts, 5)] 
    mu = np.mean(counts)
    return np.var(counts, ddof=1) / mu if mu > 0 else 0

def generate_lcdm_mock(n_sources, ra_range, dec_range):
    """Generates a Gaussian Random Field (GRF) to simulate LCDM clustering."""
    grid_res = 256
    # Create random field in Fourier space
    shape = (grid_res, grid_res)
    k = np.sqrt(np.sum(np.mgrid[:grid_res, :grid_res]**2, axis=0))
    # Power spectrum P(k) ~ k^-1.5 (approximate for galaxy clustering)
    ps = np.zeros_like(k)
    ps[k > 0] = k[k > 0]**-1.5
    
    noise = np.random.normal(size=shape) + 1j * np.random.normal(size=shape)
    field = np.fft.ifft2(np.fft.fft2(noise) * np.sqrt(ps)).real
    field = (field - field.min()) / (field.max() - field.min())
    
    # Sample sources based on the density field
    ra_samples = []
    dec_samples = []
    while len(ra_samples) < n_sources:
        ra_test = np.random.uniform(ra_range[0], ra_range[1], n_sources)
        dec_test = np.random.uniform(dec_range[0], dec_range[1], n_sources)
        
        # Map RA/Dec to grid indices
        ix = ((ra_test - ra_range[0]) / (ra_range[1] - ra_range[0]) * (grid_res-1)).astype(int)
        iy = ((dec_test - dec_range[0]) / (dec_range[1] - dec_range[0]) * (grid_res-1)).astype(int)
        
        probs = field[ix, iy]
        mask = np.random.random(n_sources) < probs
        ra_samples.extend(ra_test[mask])
        dec_samples.extend(dec_test[mask])
        
    return np.array(ra_samples[:n_sources]), np.array(dec_samples[:n_sources])

if __name__ == "__main__":
    # 1. Load Real Data Params
    path = "data/raw/COSMOS2020_subset.csv"
    real_df = pd.read_csv(path)
    N = len(real_df)
    RA_R = [real_df['ra'].min(), real_df['ra'].max()]
    DEC_R = [real_df['dec'].min(), real_df['dec'].max()]

    print(f"[*] Generating LCDM Physical Null ({N} sources)...")
    ra_lcdm, dec_lcdm = generate_lcdm_mock(N, RA_R, DEC_R)

    angles = np.linspace(0, 90, 10)
    lcdm_v = []

    print(f"[*] Auditing Rotation...")
    for theta in angles:
        r_ra, r_dec = rotate_coords(ra_lcdm, dec_lcdm, theta)
        v = get_clean_v(r_ra, r_dec)
        lcdm_v.append(v)
        print(f"Angle {theta:2.0f}° | LCDM V: {v:.4f}")

    # Results Comparison (based on your previous Real V data)
    # Real V data was: [2.23, 3.48, 4.03, 4.45, 4.62, 4.65, 4.38, 4.05, 3.50, 2.23]
    real_v = [2.2336, 3.4895, 4.0305, 4.4523, 4.6212, 4.6559, 4.3894, 4.0532, 3.5057, 2.2336]
    
    plt.figure(figsize=(10, 6))
    plt.plot(angles, real_v, 'bo-', label='COSMOS2020 (Real)')
    plt.plot(angles, lcdm_v, 'g--', label='LCDM Physical Null')
    plt.title('P6 Forensic: Real Data vs Physical Null (LCDM)')
    plt.xlabel('Rotation Angle')
    plt.ylabel('Variance Ratio (V)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('data/processed/p6_lcdm_comparison.png')
    print("[*] Comparison plot saved.")
