import numpy as np
import fitsio

def check_nz_consistency():
    print("--- P5 Robustness: N(z) Consistency Check ---")
    data = fitsio.read('data/raw/desi/LRG_NGC_clustering.dat.fits', columns=['Z'])
    rand = fitsio.read('data/raw/desi/LRG_NGC_0_clustering.ran.fits', columns=['Z'])
    
    # Check if mean redshifts are identical before weighting
    print(f"Data Mean Z (Raw): {np.mean(data['Z']):.6f}")
    print(f"Rand Mean Z (Raw): {np.mean(rand['Z']):.6f}")
    
    # If these differ by > 1e-4, the randoms aren't a true null for redshift.
    diff = abs(np.mean(data['Z']) - np.mean(rand['Z']))
    print(f"Baseline Mismatch: {diff:.6f} ({'FAIL' if diff > 1e-4 else 'PASS'})")

if __name__ == "__main__":
    check_nz_consistency()
