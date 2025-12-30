import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# CONFIG
PLANCK_MAP_PATH = 'data/raw/planck/npipe_143.fits'
HITS_MAP_PATH = 'data/raw/planck/npipe_hits.fits'
N_ROTATIONS = 100  # Number of random rotations (The "Null Universe")

def robust_audit():
    print("--- STARTING MONTE CARLO AUDIT (Spatial Correlation Fix) ---")
    
    # 1. Load Data
    print("Loading Maps...")
    m = hp.read_map(PLANCK_MAP_PATH, field=0, verbose=False) # Temp
    hits = hp.read_map(HITS_MAP_PATH, verbose=False)
    nside = hp.npix2nside(len(m))
    
    # Resample Hits if needed
    if hp.npix2nside(len(hits)) != nside:
        hits = hp.ud_grade(hits, nside)
        
    # 2. Mask & Dipole Removal (Clean the Signal)
    print("Cleaning Signal (Removing Galaxy + Dipole)...")
    th, ph = hp.pix2ang(nside, np.arange(len(m)))
    mask_gal = np.abs(90 - np.degrees(th)) > 20 # Keep |b| > 20
    
    # Create a clean map for analysis
    m_clean = np.copy(m)
    m_clean[~mask_gal] = hp.UNSEEN
    m_no_dipole = hp.remove_dipole(m_clean, fitval=False, verbose=False)
    
    # Prepare Vectors for Correlation
    # We only care about the VALID pixels common to all rotations (simplified)
    # Ideally, we rotate the mask too, but for speed we compare valid-to-valid.
    
    # Define the "Target" (Inverse Hits)
    # We fix the hits map.
    valid_mask = mask_gal # Base mask
    target_noise = 1.0 / (hits + 1.0)
    
    # 3. Measure Real Correlation
    # We only correlate pixels that are valid in the unrotated frame
    r_obs, _ = pearsonr(np.abs(m_no_dipole[valid_mask]), target_noise[valid_mask])
    print(f"\nOBSERVED CORRELATION (r_obs): {r_obs:.5f}")
    
    # 4. Monte Carlo Rotations (The Null Test)
    print(f"Running {N_ROTATIONS} random rotations to build Null Distribution...")
    null_corrs = []
    
    for i in range(N_ROTATIONS):
        if i % 10 == 0: print(f"  Rotation {i}/{N_ROTATIONS}...")
        
        # Random rotation around Z-axis (Longitude shift)
        # This breaks the alignment with the Scan Rings (which are fixed in Ecliptic)
        # while preserving the internal statistics of the CMB.
        rot_ang = np.random.uniform(10, 350) 
        
        # Rotate the MAP (keep mask fixed relative to map, so we rotate the full field)
        # Note: Rotating ALM is better but slower. rotate_alm is precise.
        # Pixel space rotation is faster for this audit.
        r = hp.Rotator(rot=[rot_ang, 0, 0], deg=True)
        m_rot = r.rotate_map_pixel(m_no_dipole)
        
        # Correlate the ROTATED map with the FIXED hits/noise
        # We use the same mask logic
        r_null, _ = pearsonr(np.abs(m_rot[valid_mask]), target_noise[valid_mask])
        null_corrs.append(r_null)
        
    # 5. Analysis
    null_mean = np.mean(null_corrs)
    null_std = np.std(null_corrs)
    z_score = (r_obs - null_mean) / null_std
    
    print("\n--- MONTE CARLO RESULTS ---")
    print(f"Null Mean: {null_mean:.5f} +/- {null_std:.5f}")
    print(f"Z-Score:   {z_score:.2f} sigma")
    
    # Plot
    plt.figure(figsize=(8, 5))
    plt.hist(null_corrs, bins=20, color='gray', alpha=0.7, label='Random Rotations (Null)')
    plt.axvline(r_obs, color='red', linewidth=3, linestyle='--', label=f'Observed ({z_score:.1f}$\sigma$)')
    plt.title("Is the P1 Correlation Real?\n(Comparing Real Data to Random Rotations)")
    plt.xlabel("Correlation Coefficient (r)")
    plt.legend()
    plt.savefig('p1_monte_carlo_result.png')
    print("Plot saved to p1_monte_carlo_result.png")

    if abs(z_score) > 3.0:
        print("VERDICT: REAL. The correlation survives rotation.")
    else:
        print("VERDICT: FAKE. The observed value falls within random chance.")

if __name__ == "__main__":
    robust_audit()
