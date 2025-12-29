#!/usr/bin/env python3
"""
P6-B Audit: The Rotational Null Test
------------------------------------
Author: The Ledger Project
Context: Distinguishing Gravity (Isotropic) from PbC (Anisotropic)

Methodology:
    1. Load/Generate Galaxy Catalog (RA, DEC).
    2. Define the fixed Survey Mask (10x10 Grid).
    3. Rotate the Galaxy coordinates by angle theta relative to the center.
    4. Re-bin the rotated galaxies into the fixed mask.
    5. Calculate Variance Ratio V(theta).
    6. If V(theta) collapses as theta increases, the structure is anisotropic (PbC).
"""

import numpy as np
import pandas as pd
import argparse
import json
import os
import sys

# --- Constants ---
GRID_SIZE = 10
OCCUPANCY_THRESHOLD = 0.90

def generate_grid_locked_data():
    """
    Generates synthetic data where clusters are ALIGNED with the grid.
    This simulates a 'PbC' signal (structure locked to the window).
    """
    print("[*] Generating Grid-Locked Synthetic Data (PbC Simulation)...")
    np.random.seed(42)
    n_patches = 100
    mu = 861
    
    # High variance distribution (Negative Binomial)
    # Variance ~ 15.4 * Mean
    r_param = mu / 14.39
    p_param = r_param / (mu + r_param)
    
    patch_counts = np.random.negative_binomial(r_param, p_param, n_patches)
    
    ra_list = []
    dec_list = []
    
    # Define Field Bounds (centered at 150.0, 2.5)
    ra_edges = np.linspace(149.0, 151.0, 11)
    dec_edges = np.linspace(1.5, 3.5, 11)
    
    for i in range(10):
        for j in range(10):
            idx = i * 10 + j
            count = patch_counts[idx]
            # Place galaxies INSIDE the grid patch
            r = np.random.uniform(ra_edges[i], ra_edges[i+1], count)
            d = np.random.uniform(dec_edges[j], dec_edges[j+1], count)
            ra_list.append(r)
            dec_list.append(d)
            
    return pd.DataFrame({
        'ra': np.concatenate(ra_list),
        'dec': np.concatenate(dec_list)
    })

def rotate_coordinates(ra, dec, theta_deg):
    """
    Rotates RA/Dec coordinates around the field center.
    """
    theta_rad = np.radians(theta_deg)
    
    # 1. Shift to origin
    ra_center = np.mean(ra)
    dec_center = np.mean(dec)
    
    x = ra - ra_center
    y = dec - dec_center
    
    # 2. Rotate
    # x' = x cos - y sin
    # y' = x sin + y cos
    x_rot = x * np.cos(theta_rad) - y * np.sin(theta_rad)
    y_rot = x * np.sin(theta_rad) + y * np.cos(theta_rad)
    
    # 3. Shift back
    return x_rot + ra_center, y_rot + dec_center

def calculate_variance_ratio(ra, dec):
    """
    Calculates the Mask-Aware Variance Ratio for a given set of coordinates.
    """
    # Define bins (Fixed Window)
    # Bounds must match the generation bounds to represent the "Camera"
    ra_bins = np.linspace(149.0, 151.0, GRID_SIZE + 1)
    dec_bins = np.linspace(1.5, 3.5, GRID_SIZE + 1)
    
    # Bin data
    hist, _, _ = np.histogram2d(ra, dec, bins=[ra_bins, dec_bins])
    
    # Mask Logic (Simulated Edge Mask)
    valid_counts = []
    
    # Define a simple mask: Inner 8x8 is valid (Occupancy=1.0), Edge is invalid (<0.9)
    # This simplifies the previous probability mask to a binary selection for clarity
    for i in range(GRID_SIZE):
        for j in range(GRID_SIZE):
            # Check if inner patch
            if 0 < i < 9 and 0 < j < 9:
                valid_counts.append(hist[i, j])
                
    valid_counts = np.array(valid_counts)
    
    # Statistics
    if len(valid_counts) == 0: return 0.0
    
    mean = np.mean(valid_counts)
    var = np.var(valid_counts, ddof=1)
    
    if mean == 0: return 0.0
    return var / mean

def run_rotation_audit(df):
    print("--- P6-B: Rotational Null Test ---")
    results = []
    
    # Sweep from 0 to 90 degrees
    angles = np.arange(0, 95, 5)
    
    v_0 = None
    
    for theta in angles:
        # 1. Rotate
        ra_rot, dec_rot = rotate_coordinates(df['ra'].values, df['dec'].values, theta)
        
        # 2. Measure
        v_ratio = calculate_variance_ratio(ra_rot, dec_rot)
        
        if theta == 0:
            v_0 = v_ratio
            
        # 3. Store
        decoherence = 0.0
        if v_0 > 0:
            decoherence = 1.0 - (v_ratio / v_0)
            
        print(f"Angle: {theta:02d} deg | V_eff: {v_ratio:.4f} | Decoherence: {decoherence*100:.1f}%")
        
        results.append({
            "angle": int(theta),
            "variance_ratio": float(v_ratio),
            "decoherence_fraction": float(decoherence)
        })
        
    return results

if __name__ == "__main__":
    # Load Data (Defaults to Synthetic Grid-Locked)
    data = generate_grid_locked_data()
    
    # Run Audit
    audit_data = run_rotation_audit(data)
    
    # Save
    output_dir = "data/processed"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "p6_rotation_results.json")
    
    with open(output_path, 'w') as f:
        json.dump(audit_data, f, indent=4)
        
    print(f"\n[*] Rotation Audit Complete. Saved to {output_path}")
