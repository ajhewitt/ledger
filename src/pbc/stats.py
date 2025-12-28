import numpy as np
import healpy as hp

def calc_phase_alignment(alm_cmb, alm_ctx, lmin=2, lmax=20):
    """
    P1 Diagnostic: Quantify phase alignment between CMB and Context.
    
    The metric S_gamma is defined as the mean cosine of the phase difference 
    weighted by the context power.
    
    References:
        PbC Dual Construction, Section 6, Diagnostic P1.
    
    Args:
        alm_cmb (np.ndarray): Complex alms of the CMB (TE or E-mode).
        alm_ctx (np.ndarray): Complex alms of the Context Template (T_ctx).
        lmin (int): Minimum multipole to consider.
        lmax (int): Maximum multipole to consider.
        
    Returns:
        float: The phase locking statistic S_gamma.
    """
    # Ensure inputs are the same size
    if alm_cmb.size != alm_ctx.size:
        raise ValueError("CMB and Context alms must have the same size.")

    # Extract phases: phi = arctan2(Im, Re)
    # Note: We assume alms are in standard healpy ordering
    phase_cmb = np.angle(alm_cmb)
    phase_ctx = np.angle(alm_ctx)
    
    # Calculate phase difference
    delta_phi = phase_cmb - phase_ctx
    
    # We weight by the magnitude of the Context to prioritize 
    # modes where the "bookkeeping" (scan strategy) is strongest.
    weights = np.abs(alm_ctx)
    
    # Get ell values for the alm array
    lmax_arr = hp.Alm.getlmax(alm_cmb.size)
    ell = hp.Alm.getlm(lmax_arr)[0]
    
    # Filter for target range [lmin, lmax]
    mask = (ell >= lmin) & (ell <= lmax)
    
    # Safety check for empty mask
    if np.sum(mask) == 0:
        return 0.0
    
    # Compute weighted mean cosine similarity
    # If phases are random (Standard Model), sum(w * cos) -> 0.
    # If phases are "locked" (Context Coupling), this sums constructively.
    numerator = np.sum(weights[mask] * np.cos(delta_phi[mask]))
    denominator = np.sum(weights[mask])
    
    if denominator == 0:
        return 0.0
        
    return numerator / denominator

def calc_galaxy_anisotropy(density_map, context_map):
    """
    P5 Diagnostic: Galaxy 2-point Anisotropy (Placeholder).
    
    Intended for use with DESI altmtl data.
    """
    raise NotImplementedError("P5 not yet implemented. Requires DESI altmtl processing.")
