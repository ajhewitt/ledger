import numpy as np
import healpy as hp
import os
import subprocess
from scipy.stats import pearsonr

# CONFIG
# The correct URL you provided
SMICA_URL = "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica_2048_R3.00_full.fits"
SMICA_PATH = "data/raw/planck/smica_cmb.fits"
HITS_PATH = "data/raw/planck/npipe_hits.fits"

def download_smica_robust():
    # 1. Check if file exists and is valid
    if os.path.exists(SMICA_PATH):
        try:
            with open(SMICA_PATH, 'rb') as f:
                header = f.read(6)
            if header == b'SIMPLE':
                print("Valid SMICA map found locally.")
                return
            else:
                print("Found corrupted file (likely HTML error page). Deleting...")
                os.remove(SMICA_PATH)
        except Exception:
            os.remove(SMICA_PATH)

    # 2. Download if missing
    print(f"Downloading SMICA Map from: {SMICA_URL}")
    print("This is a large file (~600MB). Please wait...")
    try:
        # -L follows redirects, -o output file
        subprocess.run(["curl", "-L", "-o", SMICA_PATH, SMICA_URL], check=True)
    except Exception as e:
        print(f"Download failed: {e}")
        exit()

def run_smica_audit():
    download_smica_robust()
    
    print("Loading Maps...")
    try:
        # Field 0 is Temperature (I)
        m_smica = hp.read_map(SMICA_PATH, field=0, verbose=False)
        hits = hp.read_map(HITS_PATH, verbose=False)
    except Exception as e:
        print(f"Error reading FITS: {e}")
        return

    nside = hp.npix2nside(len(m_smica))
    
    # 1. Masking
    # We apply a strict Galactic cut (|b| > 30) to be absolutely sure
    print("Masking Galaxy (|b| < 30 deg)...")
    th, ph = hp.pix2ang(nside, np.arange(len(m_smica)))
    lat_gal = 90 - np.degrees(th)
    mask = np.abs(lat_gal) > 30 
    
    # 2. Remove Dipole/Monopole
    print("Removing Monopole & Dipole (Kinematics)...")
    m_clean = np.copy(m_smica)
    m_clean[~mask] = hp.UNSEEN
    m_no_dipole = hp.remove_dipole(m_clean, fitval=False, verbose=False)
    
    # 3. Correlation Audit
    valid_idx = np.where(mask)[0]
    target_noise = 1.0 / (hits[valid_idx] + 1.0)
    
    r_obs, _ = pearsonr(np.abs(m_no_dipole[valid_idx]), target_noise)
    print(f"\nSMICA (Clean) Correlation: {r_obs:.5f}")
    
    # 4. Monte Carlo Null Test
    print("Running Null Test (N=20)...")
    nulls = []
    for i in range(20):
        rot_ang = np.random.uniform(10, 350)
        r = hp.Rotator(rot=[rot_ang, 0, 0], deg=True)
        m_rot = r.rotate_map_pixel(m_no_dipole)
        
        r_null, _ = pearsonr(np.abs(m_rot[valid_idx]), target_noise)
        nulls.append(r_null)
        print(f"  Rot {i+1}: r={r_null:.5f}")
        
    z_score = (r_obs - np.mean(nulls)) / np.std(nulls)
    print(f"\nZ-Score on Clean Data: {z_score:.2f} sigma")
    
    if abs(z_score) < 3.0:
        print("VERDICT: ARTIFACT. The signal was Zodiacal Dust.")
    else:
        print("VERDICT: PERSISTENT. The signal is cosmological.")

if __name__ == "__main__":
    run_smica_audit()
