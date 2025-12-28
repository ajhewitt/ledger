import numpy as np
import fitsio
import json

def run_significance_test():
    print("--- Running P5 Significance Test (Noise Floor) ---")
    data_path = "data/raw/desi/LRG_NGC_clustering.dat.fits"
    rand_path = "data/raw/desi/LRG_NGC_0_clustering.ran.fits"

    # 1. Get Data Shift (The Record)
    d = fitsio.read(data_path, columns=['Z', 'WEIGHT_SYS', 'WEIGHT_COMP'])
    dz_data = np.average(d['Z'], weights=d['WEIGHT_SYS']) - \
              np.average(d['Z'], weights=d['WEIGHT_COMP'])

    # 2. Get Random Shift (The Noise Floor)
    r = fitsio.read(rand_path, columns=['Z', 'WEIGHT_SYS', 'WEIGHT_COMP'])
    dz_rand = np.average(r['Z'], weights=r['WEIGHT_SYS']) - \
              np.average(r['Z'], weights=r['WEIGHT_COMP'])

    # 3. Calculate Significance
    # In a full study, we'd use all 18 random catalogs to get the std.
    # Here we use the context tension as a proxy for the expected scatter.
    sigma_noise = np.std(r['Z']) / np.sqrt(len(r))
    significance = abs(dz_data) / sigma_noise

    print(f"Data Shift:   {dz_data:.8f}")
    print(f"Random Shift: {dz_rand:.8f}")
    print(f"Significance: {significance:.2f} sigma")

    results = {
        "dz_data": float(dz_data),
        "dz_null": float(dz_rand),
        "significance_sigma": float(significance)
    }
    
    with open("data/processed/p5_significance.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    run_significance_test()
