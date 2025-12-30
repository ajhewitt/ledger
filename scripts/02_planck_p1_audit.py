import numpy as np
import healpy as hp
import os
from scipy.stats import pearsonr

# CONFIG
# Matches your download_planck_p1.sh paths
PLANCK_MAP_PATH = 'data/raw/planck/npipe_143.fits' 
HITS_MAP_PATH = 'data/raw/planck/npipe_hits.fits'

def robust_p1_audit():
    print(f"--- STARTING P1 AUDIT (Dipole-Corrected) ---")
    
    if not os.path.exists(PLANCK_MAP_PATH):
        print(f"ERROR: File not found at {PLANCK_MAP_PATH}")
        print("Run 'bash scripts/download_planck_p1.sh' first.")
        return

    # 1. Load Temperature Map (Field 0)
    print("1. Loading Planck NPIPE Temperature Map...")
    # Read field 0 (I), ignore warnings about other fields
    m = hp.read_map(PLANCK_MAP_PATH, field=0, verbose=False)
    nside = hp.npix2nside(len(m))
    
    # 2. Masking the Galaxy (Standard practice)
    print("2. Masking Galactic Plane (|b| < 20 deg)...")
    th, ph = hp.pix2ang(nside, np.arange(len(m)))
    lat = 90 - np.degrees(th)
    # Create mask: 1 where valid, 0 where masked
    mask_bool = np.abs(lat) > 20
    
    # 3. MONOPOLE & DIPOLE REMOVAL
    # We copy the map and set masked pixels to UNSEEN so remove_dipole ignores them
    m_clean = np.copy(m)
    m_clean[~mask_bool] = hp.UNSEEN
    
    print("3. Subtracting Monopole and Dipole (Kinematic Doppler)...")
    # fitval=True returns (monopole, dipole_vector)
    # This ensures we are testing the STRUCTURE (l >= 2), not the motion.
    m_no_dipole = hp.remove_dipole(m_clean, fitval=False, verbose=True)
    
    # 4. Load Context (Hits Map)
    print("4. Loading Scan Strategy (Hits Map)...")
    hits = hp.read_map(HITS_MAP_PATH, verbose=False)
    
    if hp.npix2nside(len(hits)) != nside:
        print(f"   Resampling Hits map to NSIDE {nside}...")
        hits = hp.ud_grade(hits, nside)
        
    # 5. The Correlation Test
    # PbC Hypothesis: The Record (R) is coupled to the Context (N).
    # Since N ~ 1/Hits, we check correlation between |Map| and 1/Hits.
    
    # Select only valid pixels
    valid_idx = np.where(mask_bool)[0]
    
    # We look at the magnitude of the signal vs the noise level
    signal_mag = np.abs(m_no_dipole[valid_idx])
    noise_inv = 1.0 / (hits[valid_idx] + 1.0) # Proxy for N
    
    corr, p_val = pearsonr(signal_mag, noise_inv)
    
    print("\n--- RESULTS ---")
    print(f"Map-Context Correlation: {corr:.5f}")
    print(f"P-Value:                 {p_val:.5e}")
    
    if p_val < 0.05 and abs(corr) > 0.001:
        print(">> DETECTION: The CMB structure is statistically coupled to the Scan Strategy.")
        print(">> This supports the PbC 'Phase-Lock' hypothesis.")
    else:
        print(">> NULL RESULT: The CMB structure is independent of the Scan Strategy.")

if __name__ == "__main__":
    robust_p1_audit()
