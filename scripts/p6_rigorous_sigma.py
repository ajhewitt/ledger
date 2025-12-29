import numpy as np
import pandas as pd

# Data from previous high-res audit
real_v = np.array([2.23, 3.48, 4.03, 4.45, 4.62, 4.65, 4.39, 4.05, 3.50, 2.23])
lcdm_v = np.array([1.12, 1.87, 2.09, 2.41, 2.56, 2.50, 2.37, 2.14, 1.87, 1.12])

# Empirical Error Estimation (based on N=200k sources and 100^2 bins)
# Standard error for Variance Ratio is approx sqrt(2/N_bins)
sig_v = np.sqrt(2 / 10000) # ~0.014

# Calculate Hump Height (Delta V)
dv_real = real_v[5] - real_v[0]  # 4.65 - 2.23 = 2.42
dv_lcdm = lcdm_v[5] - lcdm_v[0]  # 2.50 - 1.12 = 1.38

# Final Sigma (Difference of Deltas over pooled error)
# Note: In a full Jackknife, these errors would be calculated empirically.
# Using the theoretical floor for now:
pooled_err = np.sqrt(sig_v**2 + sig_v**2)
final_sigma = (dv_real - dv_lcdm) / pooled_err

print(f"Real Delta V: {dv_real:.2f}")
print(f"LCDM Delta V: {dv_lcdm:.2f}")
print(f"Calculated Sigma: {final_sigma:.2f} σ")
