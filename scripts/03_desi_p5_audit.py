import numpy as np
import fitsio
import json
import os
from tqdm import tqdm

def run_robust_p5_audit(n_jackknife=20):
    print("--- Running Robust P5: Self-Consistent Commutator Audit ---")
    cat_path = "data/raw/desi/LRG_NGC_clustering.dat.fits"
    
    if not os.path.exists(cat_path):
        print("Error: LRG_NGC_clustering.dat.fits not found.")
        return

    # 1. Load Data
    # 'RA' and 'DEC' are needed for the Jackknife spatial splits
    cols = ['Z', 'WEIGHT_SYS', 'WEIGHT_COMP', 'RA', 'DEC']
    data = fitsio.read(cat_path, columns=cols)
    print(f"Loaded {len(data)} objects.")

    # 2. Global Commutator Calculation
    # We compare Systematic weights against Selection (Completeness) weights
    def get_dz(d):
        z_sys = np.average(d['Z'], weights=d['WEIGHT_SYS'])
        z_comp = np.average(d['Z'], weights=d['WEIGHT_COMP'])
        return z_sys - z_comp

    global_dz = get_dz(data)

    # 3. Jackknife Resampling for Significance
    # We split the sky into 'n_jackknife' RA-bins to estimate the variance
    print(f"Performing {n_jackknife} Jackknife spatial splits...")
    ra_sorted_idx = np.argsort(data['RA'])
    splits = np.array_split(ra_sorted_idx, n_jackknife)
    
    jk_shifts = []
    for i in tqdm(range(n_jackknife)):
        # Create a "leave-one-out" mask
        mask = np.ones(len(data), dtype=bool)
        mask[splits[i]] = False
        
        # Calculate dz on the remaining (n-1) samples
        jk_dz = get_dz(data[mask])
        jk_shifts.append(jk_dz)

    # 4. Statistical Analysis
    jk_shifts = np.array(jk_shifts)
    mean_jk = np.mean(jk_shifts)
    
    # Jackknife error formula: sqrt((n-1)/n * sum((jk_i - mean_jk)^2))
    variance = ((n_jackknife - 1) / n_jackknife) * np.sum((jk_shifts - mean_jk)**2)
    error_bar = np.sqrt(variance)
    
    # Significance (Z-score)
    significance = abs(global_dz) / error_bar

    # 5. Output and Logging
    results = {
        "tracer": "LRG_NGC",
        "global_dz": float(global_dz),
        "error_bar": float(error_bar),
        "significance_sigma": float(significance),
        "n_samples": len(data),
        "status": "DETECTION" if significance > 3 else "NULL"
    }

    os.makedirs("data/processed", exist_ok=True)
    with open("data/processed/p5_robust_audit.json", "w") as f:
        json.dump(results, f, indent=4)

    print("\n" + "="*40)
    print(f"P5 Commutator Shift: {global_dz:.8f}")
    print(f"Jackknife Error:     {error_bar:.8f}")
    print(f"Significance:        {significance:.2f} sigma")
    print("="*40)

if __name__ == "__main__":
    run_robust_p5_audit()
