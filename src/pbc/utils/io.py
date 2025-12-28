import json
import numpy as np
import healpy as hp
from pathlib import Path
from dataclasses import dataclass

@dataclass
class MapBundle:
    name: str
    map: np.ndarray
    mask: np.ndarray
    nside: int

def read_map(path: str | Path, quick_nside: int | None = None) -> np.ndarray:
    """Reads a FITS map and optionally downgrades it."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Map not found: {p}")
        
    m = hp.read_map(p.as_posix(), verbose=False)
    
    # Handle multi-dimensional maps (e.g. TQU), take T if it's the first dim
    if m.ndim != 1:
        arr = np.asarray(m)
        if arr.ndim == 2 and arr.shape[0] >= 1:
            m = arr[0]
        else:
            raise ValueError(f"Unexpected FITS shape for {path}: {arr.shape}")
            
    # Downgrade if requested
    if quick_nside is not None and hp.get_nside(m) != quick_nside:
        # Power = -2 for temperature-like scalar field
        m = hp.ud_grade(m, nside_out=quick_nside, power=-2)
        
    return m

def save_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True))

def save_npy(arr: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path.as_posix(), arr)

def summary_line(msg: str) -> None:
    print(json.dumps({"msg": msg}))
