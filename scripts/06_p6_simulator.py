#!/usr/bin/env python3
import numpy as np

def run_standard_model_p6():
    print("--- P6 Simulation: Calibrating Standard Model Baseline ---")
    n_sources = 86108
    n_side = 10 
    
    # WEAKENING THE CLUSTERING:
    # Instead of 50 dense clusters, we use 2,000 'seeds' to mimic the Cosmic Web
    n_seeds = 2000 
    seeds = np.random.random((n_seeds, 2))
    
    # We assign galaxies to seeds with a wider spread (0.1 instead of 0.05)
    # This creates a 'softer' filamentary structure
    coords = []
    for _ in range(n_sources):
        center = seeds[np.random.randint(0, n_seeds)]
        coords.append(center + np.random.normal(0, 0.1, 2))
    
    coords = np.clip(np.array(coords), 0, 1)

    # Calculate Variance Ratio (V)
    ra_bins = np.linspace(0, 1, n_side + 1)
    dec_bins = np.linspace(0, 1, n_side + 1)
    counts = []
    for i in range(n_side):
        for j in range(n_side):
            m = (coords[:,0] >= ra_bins[i]) & (coords[:,0] < ra_bins[i+1]) & \
                (coords[:,1] >= dec_bins[j]) & (coords[:,1] < dec_bins[j+1])
            counts.append(np.sum(m))
    
    v_ratio = np.var(counts) / np.mean(counts)
    sigma_v = np.sqrt(2 / (n_side**2))
    significance = (v_ratio - 1) / sigma_v

    print(f"Standard Model Variance Ratio (V): {v_ratio:.4f}")
    print(f"Standard Model Significance:      {significance:.2f} sigma")

if __name__ == "__main__":
    run_standard_model_p6()
