#!/usr/bin/env python
import argparse
import numpy as np
import healpy as hp
from pathlib import Path
import json

# Import from ledger library
from pbc.stats import calc_phase_alignment
from pbc.utils.io import read_map, save_json

def main():
    parser = argparse.ArgumentParser(description="P1 Audit: Phase-Locking Diagnostic")
    parser.add_argument("--nside", type=int, default=512, help="Analysis resolution")
    parser.add_argument("--lmin", type=int, default=2)
    parser.add_argument("--lmax", type=int, default=20, help="Max multipole for phase test")
    args = parser.parse_args()

    # Define Paths
    data_dir = Path("data")
    context_path = data_dir / "processed/context_vector_c.npy"
    cmb_path = data_dir / "raw/COM_CMB_IQU-smica_2048_R3.00_full.fits"
    
    print(f"--- Running P1 Audit (Lmin={args.lmin}, Lmax={args.lmax}) ---")

    # 1. Load the Context Vector (The 'Cost')
    if not context_path.exists():
        raise FileNotFoundError(f"Context vector not found at {context_path}. Run 01_context_builder.py first.")
    
    c_alm = np.load(context_path)
    print(f"Loaded Context Vector: {c_alm.shape} modes")

    # 2. Load the Record (The 'History')
    # We read the I_STOKES (Intensity) column from our mock/real map
    print(f"Loading CMB Record: {cmb_path}")
    map_cmb = read_map(cmb_path, quick_nside=args.nside)
    
    # Convert to harmonic space (alm)
    print("Transforming Record to Harmonic Space...")
    alm_cmb = hp.map2alm(map_cmb, lmax=args.lmax)

    # 3. Compute the Diagnostic
    # This checks: "Did the cost of observing pin down the history?"
    s_gamma = calc_phase_alignment(alm_cmb, c_alm, lmin=args.lmin, lmax=args.lmax)
    
    print(f"\n>>> RESULT: P1 Phase Alignment S_gamma = {s_gamma:.5f} <<<")

    # 4. Save the Ledger Entry
    result = {
        "diagnostic": "P1_Phase_Locking",
        "s_gamma": s_gamma,
        "lmin": args.lmin,
        "lmax": args.lmax,
        "map": str(cmb_path),
        "context": str(context_path)
    }
    
    out_path = data_dir / "processed/p1_audit_result.json"
    save_json(result, out_path)
    print(f"Audit recorded in {out_path}")

if __name__ == "__main__":
    main()
