import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft

# --- Configuration ---
DATA_PATH = 'data/raw/COSMOS2020_subset.csv' # Point this to your actual file
GRID_SIZE = 50       # Binning resolution
N_SHUFFLES = 50      # Number of control realizations (Higher = smoother baseline)

def load_data(path):
    # Load data and ensure we are using the correct columns
    df = pd.read_csv(path)
    # Simple filter to remove edge garbage if necessary
    return df[['ra', 'dec']].dropna()

def get_rotational_variance(df, angle, grid_size=50):
    """
    Rotates the coordinates by 'angle' (degrees) and calculates 
    the variance of the counts in a 2D grid.
    """
    theta = np.radians(angle)
    ra_c = df['ra'].mean()
    dec_c = df['dec'].mean()
    
    # 1. Rotate
    x = (df['ra'] - ra_c) * np.cos(theta) - (df['dec'] - dec_c) * np.sin(theta)
    y = (df['ra'] - ra_c) * np.sin(theta) + (df['dec'] - dec_c) * np.cos(theta)
    
    # 2. Bin
    # Note: We don't set fixed ranges here, letting the grid float with the rotation
    # This is consistent for both data and shuffles.
    H, _, _ = np.histogram2d(x, y, bins=grid_size)
    
    # 3. Measure Variance (Clumpiness)
    # We only care about populated bins to avoid edge-zero dominance
    valid_bins = H[H > 0]
    if len(valid_bins) == 0: return 0
    
    # Index of Dispersion (Variance / Mean)
    return np.var(valid_bins) / np.mean(valid_bins)

def generate_shuffled_control(df):
    """
    Shuffles RA values to destroy physical alignment while preserving
    clustering statistics and mask interactions.
    """
    df_new = df.copy()
    df_new['ra'] = np.random.permutation(df['ra'].values)
    return df_new

# --- Execution ---

print(f"1. Loading Data from {DATA_PATH}...")
try:
    df_data = load_data(DATA_PATH)
except FileNotFoundError:
    # Fallback for testing without the file
    print("   [!] File not found. Generating Mock Data for demonstration.")
    from p6_stress_test import generate_mock_cosmos_data
    df_data = generate_mock_cosmos_data(n_points=5000)

angles = np.linspace(0, 180, 37) # Scan 0 to 180 degrees
print(f"2. Analyzing Real Data ({len(df_data)} objects)...")
signal_curve = [get_rotational_variance(df_data, a, GRID_SIZE) for a in angles]

print(f"3. Analyzing {N_SHUFFLES} Shuffled Controls (this may take a moment)...")
noise_curves = []
for i in range(N_SHUFFLES):
    df_shuf = generate_shuffled_control(df_data)
    curve = [get_rotational_variance(df_shuf, a, GRID_SIZE) for a in angles]
    noise_curves.append(curve)

# Calculate Mean Baseline and Error
baseline_curve = np.mean(noise_curves, axis=0)
baseline_std = np.std(noise_curves, axis=0)

# The "True" Signal is the difference
excess_signal = np.array(signal_curve) - baseline_curve
significance = excess_signal / baseline_std

# --- Metric Calculation ---
# Check for Mode 2 (Physical Alignment) strength in the clean signal
fft_vals = np.abs(fft(excess_signal))
mode2_strength = fft_vals[2]
mean_sig = np.mean(np.abs(significance))

print("\n--- RESULTS ---")
print(f"Mean Deviation: {mean_sig:.2f} sigma")
print(f"Dipole/Quadrupole Power: {mode2_strength:.2f}")

# --- Visualization ---
plt.figure(figsize=(10, 6))

# Top Panel: Raw Curves
plt.subplot(2, 1, 1)
plt.plot(angles, signal_curve, 'b-', label='Raw Data (W-Shape?)', linewidth=1.5)
plt.plot(angles, baseline_curve, 'k--', label='Grid Artifact Baseline', alpha=0.7)
plt.fill_between(angles, baseline_curve - baseline_std, baseline_curve + baseline_std, color='gray', alpha=0.2)
plt.title("P6 Corrected Analysis: Raw Variance vs. Grid Artifact")
plt.ylabel("Variance Index")
plt.legend()
plt.grid(True, alpha=0.3)

# Bottom Panel: The Clean Signal
plt.subplot(2, 1, 2)
plt.plot(angles, excess_signal, 'g-', linewidth=2, label='Excess (Data - Artifact)')
plt.axhline(0, color='k', linewidth=1)
# Add error bars (1-sigma from the shuffles)
plt.fill_between(angles, -baseline_std, baseline_std, color='green', alpha=0.1, label='1-Sigma Noise')
plt.title(f"The Clean Signal (Mean Deviation: {mean_sig:.2f}$\sigma$)")
plt.xlabel("Rotation Angle (Degrees)")
plt.ylabel("Residual Variance")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('p6_corrected_result.png')
plt.show()
