import numpy as np
import healpy as hp
import pytest
from pbc.stats import calc_phase_alignment

def test_p1_null_convergence():
    """
    Scientific Validation: P1 statistic must average to zero on uncorrelated skies.
    Paper Reference: Section 8, 'Null' protocol.
    """
    nsims = 50
    lmax = 16
    nside = 32
    
    # 1. Create a fixed arbitrary Context (e.g., a dipole or random map)
    # The exact shape doesn't matter for the null test, only that it's fixed.
    cl = np.ones(lmax + 1)
    cl[0:2] = 0
    alm_ctx = hp.synalm(cl, lmax=lmax, new=True)
    
    scores = []
    
    # 2. Run Null Suite (Lambda = 0)
    for _ in range(nsims):
        # Generate random CMB (uncorrelated with Context)
        alm_cmb = hp.synalm(cl, lmax=lmax, new=True)
        
        # Measure alignment
        score = calc_phase_alignment(alm_cmb, alm_ctx, lmin=2, lmax=lmax)
        scores.append(score)
        
    scores = np.array(scores)
    
    # 3. Assertions
    mean_score = np.mean(scores)
    std_error = np.std(scores) / np.sqrt(nsims)
    
    # The mean should be within 3 sigma of 0
    # (If this fails, the estimator is biased)
    assert np.abs(mean_score) < (3 * std_error), \
        f"P1 Null Test Failed: Mean {mean_score:.4f} is significantly non-zero (SE={std_error:.4f})"

def test_p1_perfect_alignment():
    """
    Sanity Check: If CMB == Context, score should be exactly 1.0.
    """
    lmax = 10
    cl = np.ones(lmax + 1)
    cl[0:2] = 0
    
    alm_ctx = hp.synalm(cl, lmax=lmax, new=True)
    
    # Test perfect alignment
    score = calc_phase_alignment(alm_ctx, alm_ctx, lmin=2, lmax=lmax)
    
    assert np.isclose(score, 1.0), f"P1 failed identity test. Expected 1.0, got {score}"

def test_p1_anti_alignment():
    """
    Sanity Check: If CMB == -Context, score should be exactly -1.0 (phase diff = pi).
    """
    lmax = 10
    cl = np.ones(lmax + 1)
    cl[0:2] = 0
    
    alm_ctx = hp.synalm(cl, lmax=lmax, new=True)
    
    # Test perfect anti-alignment
    # Note: In harmonic space, -alm means phase shift of pi
    score = calc_phase_alignment(-alm_ctx, alm_ctx, lmin=2, lmax=lmax)
    
    assert np.isclose(score, -1.0), f"P1 failed inversion test. Expected -1.0, got {score}"
