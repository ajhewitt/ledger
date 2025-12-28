#!/usr/bin/env python3
import numpy as np
import fitsio
import json
import os

def run_robust_p6_audit():
    print("--- Running P6: COSMOS Field Variance Audit (Final Pillar) ---")
    cat_path = "data/raw/jwst/cosmos2020_classic.fits.gz"
    
    if not os.path.exists(cat_path):
        print(f"Catalog not found at {cat_path}.")
        return

    # 1. Open the file to inspect all columns
    fits = fitsio.FITS(cat_path)
    all_cols = fits[1].get_colnames()
    
    # Identify the best redshift column (LePhare is preferred for COSMOS2020)
    z_col = None
    for candidate in ['lp_zBEST', 'LP_ZBEST', 'z_phot', 'ZPDF']:
        if candidate in all_cols:
            z_col = candidate
            break
    
    if not z_col:
        # Fallback: search for anything with 'z' and 'best' or 'phot'
        z_col = [c for c in all_cols if 'z' in c.lower() and ('best' in c.lower() or 'phot' in c.lower())][0]

    print(f"Using Columns: ALPHA_J2000, DELTA_J2000, {z_col}")

    # 2. Load the Record
    data = fits[1].read(columns=['ALPHA_J2000', 'DELTA_J2000', z_col])
    
    # 3. Filter for Transition Zone (z > 3)
    # This is the "Goldilocks" zone between the coupled CMB and balanced LSS
    mask = (data[z_col] > 3.0) & (data[z_col] < 6.0)
    deep_data = data[mask]
    print(f"Auditing {len(deep_data)} high-redshift sources.")

    # 4. Spatial Patching (Field Variance)
    n_side = 4
    ra_bins = np.linspace(deep_data['ALPHA_J2000'].min(), deep_data['ALPHA_J2000'].max(), n_side + 1)
    dec_bins = np.linspace(deep_data['DELTA_J2000'].min(), deep_data['DELTA_J2000'].max(), n_side + 1)
    
    counts = []
    for i in range(n_side):
        for j in range(n_side):
            patch_mask = (deep_data['ALPHA_J2000'] >= ra_bins[i]) & (deep_data['ALPHA_J2000'] < ra_bins[i+1]) & \
                         (deep_data['DELTA_J2000'] >= dec_bins[j]) & (deep_data['DELTA_J2000'] < dec_bins[j+1])
            counts.append(np.sum(patch_mask))
    
    counts = np.array(counts)
    
    # 5. Calculate Variance Ratio (V)
    mean_count = np.mean(counts)
    v_ratio = np.var(counts) / mean_count
    sigma_v = np.sqrt(2 / (n_side**2))
    significance = (v_ratio - 1) / sigma_v

    print(f"Variance Ratio (V):   {v_ratio:.4f}")
    print(f"P6 Significance:      {significance:.2f} sigma")

    results = {
        "track": "P6_JWST",
        "redshift_col": z_col,
        "v_ratio": float(v_ratio),
        "significance": float(significance),
        "status": "COUPLED" if significance > 3 else "BALANCED"
    }

    os.makedirs("data/processed", exist_ok=True)
    with open("data/processed/p6_audit_result.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    run_robust_p6_audit()
