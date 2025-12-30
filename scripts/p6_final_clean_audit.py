import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import rotate
from scipy.fft import fft
import os

# CONFIG
DATA_PATH = 'data/raw/COSMOS2020_subset.csv'
GRID_SIZE = 50       # Resolution of the density grid
N_SHUFFLES = 20      # Number of null tests (keep low for speed, raise to 100 for paper)

def remove_density_dipole(grid):
    """Fits and subtracts a 2D plane (Gradient) from a density grid."""
    ny, nx = grid.shape
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    
    x_flat = X.ravel()
    y_flat = Y.ravel()
    z_flat = grid.ravel()
    
    # Fit plane: z = ax + by + c
    def plane(coords, a, b, c):
        x, y = coords
        return a*x + b*y + c
    
    # Fit only to valid data (nonzero) to avoid edge bias
    mask = z_flat > 0
    if np.sum(mask) < 10: return grid # Too sparse to fit
    
    popt, _ = curve_fit(plane, (x_flat[mask], y_flat[mask]), z_flat[mask])
    
    # Subtract plane from whole grid
    z_plane = plane((x_flat, y_flat), *popt).reshape(ny, nx)
    return grid - z_plane

def get_grid_variance(grid, angle):
    """Rotates the grid and measures the variance of the structure."""
    # reshape=False keeps the "Window" static (like the telescope detector)
    rot_grid = rotate(grid, angle, reshape=False, mode='constant', cval=0)
    
    # Measure variance of the signal inside the window
    valid_pixels = rot_grid[rot_grid != 0]
    if len(valid_pixels) == 0: return 0
    
    # Index of Dispersion (Clumpiness)
    return np.var(valid_pixels) / (np.mean(valid_pixels) + 1e-9)

def p6_clean_audit():
    print("--- STARTING P6 AUDIT (Dipole-Corrected) ---")
    
    if not os.path.exists(DATA_PATH):
        print(f"ERROR: {DATA_PATH} not found.")
        return

    # 1. Load Data
    df = pd.read_csv(DATA_PATH)
    print(f"Loaded {len(df)} objects.")
    
    # 2. Grid the Real Data
    print("Gridding and removing dipole from Real Data...")
    H_data, _, _ = np.histogram2d(df['ra'], df['dec'], bins=GRID_SIZE)
    H_data = H_data.T # Transpose for image coords
    H_clean = remove_density_dipole(H_data)
    
    # 3. Analyze Real Data
    angles = np.linspace(0, 180, 37)
    real_curve = [get_grid_variance(H_clean, a) for a in angles]
    
    # 4. Analyze Shuffles (The Null Hypothesis)
    print(f"Running {N_SHUFFLES} Shuffled Null Tests...")
    shuffle_curves = []
    
    for i in range(N_SHUFFLES):
        # Shuffle RA to destroy physical alignment but keep mask/density
        df_shuf = df.copy()
        df_shuf['ra'] = np.random.permutation(df['ra'].values)
        
        # Grid -> Remove Dipole -> Measure
        H_shuf, _, _ = np.histogram2d(df_shuf['ra'], df_shuf['dec'], bins=GRID_SIZE)
        H_shuf_clean = remove_density_dipole(H_shuf.T)
        
        curve = [get_grid_variance(H_shuf_clean, a) for a in angles]
        shuffle_curves.append(curve)
        
    # 5. Statistics
    baseline_mean = np.mean(shuffle_curves, axis=0)
    baseline_std = np.std(shuffle_curves, axis=0)
    
    # Excess Signal (Data - Null)
    excess = np.array(real_curve) - baseline_mean
    
    # Harmonic Analysis
    fft_vals = np.abs(fft(excess))
    mode2 = fft_vals[2] # Quadrupole (Physical Alignment)
    mode4 = fft_vals[4] # Grid Artifact (Square)
    
    print("\n--- RESULTS ---")
    print(f"Mode 2 (Signal):   {mode2:.2f}")
    print(f"Mode 4 (Artifact): {mode4:.2f}")
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(angles, excess, 'g-', linewidth=2, label='Clean Excess Signal')
    plt.fill_between(angles, -baseline_std, baseline_std, color='gray', alpha=0.3, label='1-sigma (Null)')
    plt.title(f"P6 Clean Audit: Signal vs Grid Artifact\nMode2={mode2:.1f} | Mode4={mode4:.1f}")
    plt.xlabel("Rotation Angle")
    plt.ylabel("Variance Excess")
    plt.legend()
    plt.savefig('p6_clean_result.png')
    print("Plot saved to p6_clean_result.png")

if __name__ == "__main__":
    p6_clean_audit()
