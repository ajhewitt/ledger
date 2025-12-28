import math
import numpy as np
import pymaster as nmt
from typing import Any, Mapping

def create_bins(nside: int, nlb: int = 50, lmin: int | None = None, lmax: int | None = None) -> nmt.NmtBin:
    """
    Construct linear binning scheme.
    """
    if nlb <= 0:
        raise ValueError("nlb must be positive")
        
    kwargs = {"nside": nside, "nlb": nlb}
    if lmin is not None: kwargs["lmin"] = int(lmin)
    if lmax is not None: kwargs["lmax"] = int(lmax)

    try:
        return nmt.NmtBin.from_nside_linear(**kwargs)
    except TypeError:
        # Fallback for older NaMaster versions
        lmin_val = 0 if lmin is None else int(lmin)
        lmax_val = 3 * nside - 1 if lmax is None else int(lmax)
        
        ell_min = np.arange(lmin_val, lmax_val + 1, nlb, dtype=int)
        ell_max = np.minimum(ell_min + nlb - 1, lmax_val)
        return nmt.NmtBin.from_edges(ell_min, ell_max)

def create_field(m: np.ndarray, mask: np.ndarray, lmax: int | None = None) -> nmt.NmtField:
    """Creates a NaMaster scalar field."""
    # Ensure mask is applied
    m_masked = m * mask
    return nmt.NmtField(mask, [m_masked], lmax=lmax)

def compute_bandpowers(f1: nmt.NmtField, f2: nmt.NmtField, b: nmt.NmtBin) -> np.ndarray:
    """Compute workspace-based bandpowers."""
    w = nmt.NmtWorkspace()
    w.compute_coupling_matrix(f1, f2, b)
    cl_coupled = nmt.compute_coupled_cell(f1, f2)
    return w.decouple_cell(cl_coupled)
