import numpy as np
import fitsio
import json
import os

def run_p5_audit():
    print("--- Running P5: Target Selection Commutator Audit ---")
    cat_path = "data/raw/desi/LRG_NGC_clustering.dat.fits"
    
    if not os.path.exists(cat_path):
        print("Data not found. Run scripts/download_desi_p5.sh first.")
        return

    # 1. Load the specific columns from the iron v1.5 ledger
    data = fitsio.read(cat_path, columns=['Z', 'WEIGHT', 'WEIGHT_SYS', 'WEIGHT_COMP'])
    
    # 2. Compute the Commutator shift: 
    # Difference between 'Selection-only' and 'Systematic-only' reconstructions
    # This is a proxy for the P5 commutator [Selection, Inference]
    z_sys = np.average(data['Z'], weights=data['WEIGHT_SYS'])
    z_comp = np.average(data['Z'], weights=data['WEIGHT_COMP'])
    
    # The P5 shift: Does the bookkeeping of fiber assignment shift the physics?
    shift = z_sys - z_comp
    
    # 3. Calculate "Context Tension"
    # How much extra variance is the fiber assignment adding?
    tension = np.std(data['WEIGHT_COMP'] / data['WEIGHT_SYS'])
    
    results = {
        "tracer": "LRG_NGC",
        "commutator_z_shift": float(shift),
        "context_tension": float(tension),
        "status": "PASS" if abs(shift) < 1e-4 else "INVESTIGATE"
    }
    
    os.makedirs("data/processed", exist_ok=True)
    with open("data/processed/p5_audit_result.json", "w") as f:
        json.dump(results, f, indent=4)
    
    print(f"P5 Commutator Shift (z): {shift:.8f}")
    print(f"Context Tension (Selection/Sys): {tension:.4f}")

if __name__ == "__main__":
    run_p5_audit()
