#!/usr/bin/env python3
import numpy as np
import fitsio
import os

def verify_p6_geometry():
    print("--- P6 Verification: Grid & Mask Audit ---")
    cat_path = "data/raw/jwst/cosmos2020_classic.fits.gz"
    
    # We load RA/Dec and a Flag column to see what's 'Good' data
    # FLAG_COMBINED = 0 usually means clean data area
    data = fitsio.read(cat_path, columns=['ALPHA_J2000', 'DELTA_J2000', 'lp_zBEST', 'FLAG_COMBINED'])
    
    mask = (data['lp_zBEST'] > 3.0) & (data['lp_zBEST'] < 6.0)
    clean_data = data[mask & (data['FLAG_COMBINED'] == 0)]
    
    print(f"Clean sources: {len(clean_data)} (out of {np.sum(mask)})")

    # Audit the Grid
    n_side = 4
    ra_bins = np.linspace(clean_data['ALPHA_J2000'].min(), clean_data['ALPHA_J2000'].max(), n_side + 1)
    dec_bins = np.linspace(clean_data['DELTA_J2000'].min(), clean_data['DELTA_J2000'].max(), n_side + 1)
    
    counts = []
    for i in range(n_side):
        for j in range(n_side):
            p_mask = (clean_data['ALPHA_J2000'] >= ra_bins[i]) & (clean_data['ALPHA_J2000'] < ra_bins[i+1]) & \
                     (clean_data['DELTA_J2000'] >= dec_bins[j]) & (clean_data['DELTA_J2000'] < dec_bins[j+1])
            c = np.sum(p_mask)
            print(f"Patch ({i},{j}) Count: {c}")
            counts.append(c)
    
    counts = np.array(counts)
    v_ratio = np.var(counts) / np.mean(counts)
    print(f"\nCorrected Variance Ratio: {v_ratio:.4f}")

if __name__ == "__main__":
    verify_p6_geometry()
