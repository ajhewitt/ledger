import numpy as np
import healpy as hp
from pathlib import Path
from typing import Optional, List, Tuple

from pbc.utils.io import read_map

def build_mask(
    map_in: np.ndarray, 
    threshold_sigma: float = 10.0, 
    apod_arcmin: float = 30.0
) -> np.ndarray:
    """
    Constructs a binary or apodized mask from a data map.
    
    Salvaged from 'make_mask.py'. Used to ensure the context template 
    is defined on the same valid sky region as the CMB analysis.
    
    Args:
        map_in: Input map (e.g., Planck SMICA or Lensing).
        threshold_sigma: Cutoff for outlier removal (point sources).
        apod_arcmin: Apodization scale in arcminutes (C2 window).
        
    Returns:
        np.ndarray: The processed mask (float 0..1).
    """
    # 1. Thresholding (Z-score cut)
    mu = np.median(map_in)
    sig = np.std(map_in)
    binary_mask = np.abs(map_in - mu) < (threshold_sigma * sig)
    
    mask_float = binary_mask.astype(float)
    
    # 2. Apodization
    # Note: We use healpy's smoothing as a lightweight alternative 
    # if NaMaster is not strictly required for this step, 
    # but strictly this should align with the paper's window function.
    if apod_arcmin > 0:
        # Simple Gaussian smoothing approximation for apodization edges
        # In a full pipeline, use nmt.mask_apodization here.
        fwhm_rad = np.radians(apod_arcmin / 60.0)
        mask_apod = hp.smoothing(mask_float, fwhm=fwhm_rad)
        # Re-clip to ensure strictly [0,1]
        mask_apod = np.clip(mask_apod, 0.0, 1.0)
        # Hard zero out the original bad pixels
        mask_apod *= mask_float 
        return mask_apod
    
    return mask_float

def effective_f_sky(mask: np.ndarray) -> float:
    """Returns the effective sky fraction f_sky."""
    return np.mean(mask)

def build_context_template(
    exposure_path: str,
    zodi_path: str,
    mask_path: Optional[str] = None,
    systematics_paths: List[str] = [],
    nside: int = 2048,
    lmax: int = 20
) -> np.ndarray:
    """
    Implements the 'Context Template Construction' algorithm (Section 5).
    
    This builds the basis vectors that define the 'Parochial' prior penalty.
    It projects housekeeping data into Ecliptic coordinates and extracts
    the low-ell structural modes.
    
    Args:
        exposure_path: Path to N_obs or hitcount map.
        zodi_path: Path to Zodiacal emission template.
        mask_path: Path to analysis mask (optional, generates one if None).
        systematics_paths: List of paths to known systematics to project out 
                           (Step 2 of Algorithm [cite: 232]).
        nside: Processing resolution.
        lmax: Maximum multipole for the feature basis (default 20 for low-ell).
        
    Returns:
        np.ndarray: The primary context vector 'c' (alm coefficients), 
                    normalized such that P_ctx = c c.T / |c|^2.
    """
    # --- Step 1: Feature Extraction  ---
    # Load raw housekeeping maps
    exp_map = read_map(exposure_path, quick_nside=nside)
    zodi_map = read_map(zodi_path, quick_nside=nside)
    
    if mask_path:
        mask = read_map(mask_path, quick_nside=nside)
    else:
        # Fallback: build mask from exposure map stability
        mask = build_mask(exp_map)

    # Coordinate Rotation: Galactic -> Ecliptic
    # The "bookkeeping" (scan strategy) is locked to Ecliptic coordinates.
    # We must rotate the maps to capture the geometry correctly.
    r = hp.Rotator(coord=['G', 'E'])
    exp_ecl = r.rotate_map_pixel(exp_map)
    zodi_ecl = r.rotate_map_pixel(zodi_map)
    
    # Extract low-ell harmonics on the masked sky
    # We effectively treat the masked region as zeros for feature definition
    alm_exp = hp.map2alm(exp_ecl * mask, lmax=lmax)
    alm_zodi = hp.map2alm(zodi_ecl * mask, lmax=lmax)
    
    # Construct Design Matrix X (Rows = modes, Cols = features)
    # We stack the alm vectors as columns.
    # Note: hp.map2alm returns complex numbers. We treat Real/Imag 
    # as distinct degrees of freedom for the SVD if we want a real basis,
    # but typically we work in the complex alm space for the projector P_ctx.
    # For this implementation, we simply stack the vectors.
    X = np.vstack([alm_exp, alm_zodi]).T  # Shape: (N_modes, 2)
    
    # --- Step 2: Orthogonalization [cite: 232] ---
    # "Orthogonalize against known systematics (beam asymmetries...)"
    if systematics_paths:
        sys_alms = []
        for sys_path in systematics_paths:
            m_sys = read_map(sys_path, quick_nside=nside)
            # Systematics are usually in Galactic, so rotate them too
            m_sys_ecl = r.rotate_map_pixel(m_sys)
            sys_alms.append(hp.map2alm(m_sys_ecl * mask, lmax=lmax))
        
        S = np.vstack(sys_alms).T
        
        # Projection: X_clean = (I - S(S^T S)^-1 S^T) X
        # Using least squares to subtract the systematics component
        for i in range(X.shape[1]):
            coeffs, _, _, _ = np.linalg.lstsq(S, X[:, i], rcond=None)
            X[:, i] -= S @ coeffs

    # --- Step 3 & 4: Normalization and Basis Selection  ---
    # "Take the leading mode c" via SVD
    U, s, Vh = np.linalg.svd(X, full_matrices=False)
    
    # The first left-singular vector is our primary context direction 'c'
    # This vector encapsulates the dominant correlation between Exposure and Zodi
    c_alm = U[:, 0]
    
    return c_alm
