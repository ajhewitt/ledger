#!/usr/bin/env python
import argparse
import numpy as np
import healpy as hp
from pathlib import Path
from tqdm import tqdm

from pbc.stats import calc_phase_alignment
from pbc.utils.io import save_json

def generate_dummy_spectrum(lmax):
    """Generates a flat Cl spectrum for testing physics-agnostic statistics."""
    ell = np.arange(lmax + 1)
    cl = np.ones_like(ell, dtype=float)
    cl[0:2] = 0  # No monopole/dipole
    return cl

def main():
    parser = argparse.ArgumentParser(description="P1 Injection Test Suite")
    parser.add_argument("--nsims", type=int, default=100, help="Number of Monte Carlo realizations")
    parser.add_argument("--nside", type=int, default=64, help="Healpix resolution (keep low for speed)")
    parser.add_argument("--lmax", type=int, default=32, help="Maximum multipole for P1 test")
    parser.add_argument("--lambda-inj", type=float, default=0.5, help="Amplitude of injected context signal")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=Path, default="artifacts/injection_results.json")
    args = parser.parse_args()

    np.random.seed(args.seed)
    
    # 1. Setup Context Template (T_ctx)
    # In a real run, load this from 'pbc.context.build_context_template'
    # Here, we generate a random 'context' vector to test the statistic against
    print(f"Generating synthetic context (Lmax={args.lmax})...")
    cl_ctx = generate_dummy_spectrum(args.lmax)
    # T_ctx is usually geometric, so let's make it static (fixed across sims)
    alm_ctx = hp.synalm(cl_ctx, lmax=args.lmax, new=True)

    # 2. Simulation Loop
    results = []
    print(f"Running {args.nsims} injections with lambda={args.lambda_inj}...")
    
    for i in tqdm(range(args.nsims)):
        # A. Generate Null CMB (random phases)
        alm_cmb = hp.synalm(cl_ctx, lmax=args.lmax, new=True)
        
        # B. Inject Signal (The "Cost" of Observation)
        # a_obs = a_cmb + lambda * a_ctx
        # Note: This simulates the shift in mean (mu_CE - mu_TE)
        alm_obs = alm_cmb + (args.lambda_inj * alm_ctx)
        
        # C. Measure P1 Statistic (S_gamma)
        # We check phase alignment between the *Observed* sky and the *Context*
        s_gamma = calc_phase_alignment(alm_obs, alm_ctx, lmin=2, lmax=args.lmax)
        
        results.append(s_gamma)

    # 3. Analysis
    results = np.array(results)
    mean_score = np.mean(results)
    std_score = np.std(results)
    
    print(f"\n--- Results (Lambda_inj = {args.lambda_inj}) ---")
    print(f"Mean S_gamma: {mean_score:.4f} +/- {std_score/np.sqrt(args.nsims):.4f}")
    
    # 4. Save
    output = {
        "nsims": args.nsims,
        "lambda_inj": args.lambda_inj,
        "nside": args.nside,
        "mean_s_gamma": mean_score,
        "std_s_gamma": std_score,
        "raw_scores": results.tolist()
    }
    save_json(output, args.out)
    print(f"Saved artifacts to {args.out}")

if __name__ == "__main__":
    main()
