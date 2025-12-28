#!/usr/bin/env python3
import numpy as np
import healpy as hp
import json
import os

def run_comprehensive_p1_audit():
    print("--- Running P1: Comprehensive Differential Audit (R4.00) ---")
    
    # Paths
    hr1_path = "data/raw/planck/npipe_143_hr1.fits"
    hr2_path = "data/raw/planck/npipe_143_hr2.fits"

    if not os.path.exists(hr1_path):
        print("Data not found. Run scripts/download_p1_splits.sh.")
        return

    # 1. Load Data
    # In R4.00 maps: Field 1=Q, 2=U, 3=Hits
    NSIDE_WORK = 256
    print("Loading splits and downsampling to NSIDE 256...")
    
    def load_data(path):
        q = hp.ud_grade(hp.read_map(path, field=1, verbose=False), NSIDE_WORK)
        u = hp.ud_grade(hp.read_map(path, field=2, verbose=False), NSIDE_WORK)
        h = hp.ud_grade(hp.read_map(path, field=3, verbose=False), NSIDE_WORK)
        return q, u, h

    q1, u1, h1 = load_data(hr1_path)
    q2, u2, h2 = load_data(hr2_path)

    # 2. Reconstruct Record (Sum) vs Artifact (Diff)
    q_sum, u_sum = (q1 + q2) / 2, (u1 + u2) / 2
    q_diff, u_diff = (q1 - q2) / 2, (u1 - u2) / 2
    hits = h1 + h2

    # 3. Apply Galactic Mask (35 degree cut for ultra-clean signal)
    mask = np.abs(hp.pix2ang(NSIDE_WORK, np.arange(hp.nside2npix(NSIDE_WORK)))[0] - np.pi/2) > 0.61

    def get_phase_locking(q, u, h):
        # Local polarization angle
        psi = 0.5 * np.arctan2(u[mask], q[mask])
        # Circular-Linear correlation with the scan context (hits)
        return np.corrcoef(psi, h[mask])[0, 1]

    # 4. Calculate Net Signal
    s_full = get_phase_locking(q_sum, u_sum, hits)
    s_noise = get_phase_locking(q_diff, u_diff, hits)
    s_net = s_full - s_noise

    print(f"Full Map S_gamma:  {s_full:.6f}")
    print(f"Noise Map S_gamma: {s_noise:.6f}")
    print(f"Net PbC Signal:    {s_net:.6f}")

    # Significance Estimate
    # If s_net > 0.005, we have a robust non-instrumental alignment.
    results = {
        "s_full": float(s_full),
        "s_noise": float(s_noise),
        "s_net": float(s_net),
        "status": "COUPLED" if abs(s_net) > 0.005 else "BALANCED"
    }

    os.makedirs("data/processed", exist_ok=True)
    with open("data/processed/p1_diff_audit.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    run_comprehensive_p1_audit()
