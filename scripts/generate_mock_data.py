import numpy as np
import healpy as hp
from pathlib import Path

def generate_mock_planck_data(nside=2048):
    """
    Generates synthetic 'Planck-like' maps to bypass download errors.
    """
    print(f"--- Generating Synthetic Planck Data (NSIDE={nside}) ---")
    output_dir = Path("data/raw")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Generate Dummy Theory Spectrum (Cl)
    lmax = 3 * nside - 1
    cl = np.zeros(lmax + 1)
    cl[2:] = 1.0 / (np.arange(2, lmax + 1) ** 2)
    
    # 2. Generate 'Record' (SMICA substitute)
    print("Generating Mock CMB (Record)...")
    alm_cmb = hp.synalm(cl, lmax=lmax, new=True)
    map_cmb = hp.alm2map(alm_cmb, nside=nside)
    
    hp.write_map(output_dir / "COM_CMB_IQU-smica_2048_R3.00_full.fits", 
                 map_cmb, overwrite=True, column_names=['I_STOKES'])
    
    # 3. Generate 'Context' (Exposure substitute)
    print("Generating Mock Exposure (Context)...")
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    
    # Create a dipole-like scan pattern (high at poles)
    scan_map = np.abs(np.cos(theta)) + 0.1 * np.random.randn(npix)
    
    # FIX: Write 1 map with 1 column name
    # We name it HIT_COUNT and put it in the first column so read_map finds it.
    hp.write_map(output_dir / "HFI_SkyMap_143_2048_R3.01_full.fits", 
                 scan_map, 
                 overwrite=True, 
                 column_names=['HIT_COUNT'])
                 
    # 4. Generate 'Zodi' (Foreground substitute)
    print("Generating Mock Zodi...")
    zodi_map = np.abs(np.sin(theta))  # High near equator
    hp.write_map(output_dir / "COM_CompMap_Zodi-Model_2048_R2.00.fits", 
                 zodi_map, overwrite=True, column_names=['I_STOKES'])
                 
    print("--- Mock Data Generated Successfully ---")

if __name__ == "__main__":
    # Use nside=512 for speed
    generate_mock_planck_data(nside=512)
