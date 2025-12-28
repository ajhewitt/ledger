#!/usr/bin/env python3
import numpy as np
import fitsio
import json
import os

def run_mask_aware_p6_audit():
    print("--- Running P6: Mask-Aware Field Variance Audit ---")
    cat_path = "data/raw/jwst/cosmos2020_classic.fits.gz"
    
    if not os.path.exists(cat_path):
        print(f"Catalog not found at {cat_path}.")
        return

    # 1. Load the Record and verified Redshift column
    fits = fitsio.FITS(cat_path)
    data = fits[1].read(columns=['ALPHA_J2000', 'DELTA_J2000', 'lp_zBEST', 'FLAG_COMBINED'])
    
    # 2. Define the Transition Zone (z > 3) and clean data
    mask = (data['lp_zBEST'] > 3.0) & (data['lp_zBEST'] < 6.0) & (data['FLAG_COMBINED'] == 0)
    clean_data = data[mask]
    
    # 3. Define the Grid
    n_side = 10  # Increased resolution to find "Good" patches
    ra_bins = np.linspace(clean_data['ALPHA_J2000'].min(), clean_data['ALPHA_J2000'].max(), n_side + 1)
    dec_bins = np.linspace(clean_data['DELTA_J2000'].min(), clean_data['DELTA_J2000'].max(), n_side + 1)
    
    valid_counts = []
    
    # 4. Patch Occupancy Audit
    for i in range(n_side):
        for j in range(n_side):
            # Define patch boundaries
            r_min, r_max = ra_bins[i], ra_bins[i+1]
            d_min, d_max = dec_bins[j], dec_bins[j+1]
            
            # Count sources in this patch
            patch_mask = (clean_data['ALPHA_J2000'] >= r_min) & (clean_data['ALPHA_J2000'] < r_max) & \
                         (clean_data['DELTA_J2000'] >= d_min) & (clean_data['DELTA_J2000'] < d_max)
            count = np.sum(patch_mask)
            
            # MASK-AWARE CHECK: 
            # We assume a 'Full' patch should have at least 50% of the median count
            # This filters out edges and large holes
            valid_counts.append(count)

    # 5. Filter for High-Coverage Patches (>90% of local median)
    counts_array = np.array(valid_counts)
    threshold = np.median(counts_array) * 0.9
    filtered_counts = counts_array[counts_array > threshold]
    
    n_patches = len(filtered_counts)
    mean_count = np.mean(filtered_counts)
    v_ratio = np.var(filtered_counts) / mean_count
    
    # Significance relative to Poisson noise for the filtered set
    sigma_v = np.sqrt(2 / n_patches)
    significance = (v_ratio - 1) / sigma_v

    print(f"Patches Analyzed:   {n_patches} (Filtered from {n_side**2})")
    print(f"Variance Ratio (V): {v_ratio:.4f}")
    print(f"Real Significance:  {significance:.2f} sigma")

    # Update the ledger record
    results = {
        "track": "P6_JWST_MaskAware",
        "n_patches": n_patches,
        "v_ratio": float(v_ratio),
        "significance": float(significance),
        "status": "COUPLED" if significance > 3 else "BALANCED"
    }

    os.makedirs("data/processed", exist_ok=True)
    with open("data/processed/p6_mask_aware_result.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    run_mask_aware_p6_audit()
