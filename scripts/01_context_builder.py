#!/usr/bin/env python
import argparse
import numpy as np
import healpy as hp
from pathlib import Path
import matplotlib.pyplot as plt

# Import from your new library
from pbc.context import build_context_template
from pbc.utils.io import save_npy

def plot_context_vector(c_alm, nside, title, outfile):
    """Optional: Visualizes the context template in map space."""
    try:
        # Inverse harmonic transform to see what the 'Cost' looks like
        m = hp.alm2map(c_alm, nside)
        plt.figure(figsize=(10, 5))
        hp.mollview(m, title=title, hold=True)
        plt.savefig(outfile)
        plt.close()
        print(f"Saved visualization to {outfile}")
    except Exception as e:
        print(f"Skipping visualization (matplotlib error): {e}")

def main():
    parser = argparse.ArgumentParser(description="Construct the Context Template (T_ctx)")
    parser.add_argument("--nside", type=int, default=512, help="Analysis resolution (use 512 for mock data)")
    parser.add_argument("--lmax", type=int, default=20, help="Maximum multipole for Context Definition")
    args = parser.parse_args()

    # Paths to the 'Record' and 'Context' files
    # (These match the filenames from download_planck.sh / generate_mock_data.py)
    raw_dir = Path("data/raw")
    exposure_path = raw_dir / "HFI_SkyMap_143_2048_R3.01_full.fits"
    zodi_path = raw_dir / "COM_CompMap_Zodi-Model_2048_R2.00.fits"
    
    output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"--- Building Context Template (NSIDE={args.nside}, LMAX={args.lmax}) ---")
    
    # 1. Execute the Protocol (Section 5)
    # This rotates the scan strategy to Ecliptic coords and finds the shared structure with Zodi
    print("Running SVD on Exposure/Zodi alignment...")
    c_alm = build_context_template(
        exposure_path=str(exposure_path),
        zodi_path=str(zodi_path),
        nside=args.nside,
        lmax=args.lmax
    )
    
    # 2. Save the Vector
    # This 'c_alm' is the "fingerprint" of the observer.
    out_path = output_dir / "context_vector_c.npy"
    save_npy(c_alm, out_path)
    print(f"Context Vector 'c' saved to: {out_path}")
    
    # 3. Visualize
    # We want to see what the "Observer's Bias" looks like on the sky
    plot_path = output_dir / "context_map_view.png"
    plot_context_vector(c_alm, args.nside, "The Cost of Observation (T_ctx)", plot_path)

if __name__ == "__main__":
    main()
