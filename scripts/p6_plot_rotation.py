#!/usr/bin/env python3
import matplotlib.pyplot as plt
import json
import pandas as pd
import numpy as np
import os

def plot_rotation_curve(json_path, output_path):
    print(f"[*] Loading data from {json_path}...")
    
    # 1. Load Data
    with open(json_path, 'r') as f:
        data = json.load(f)
    df = pd.DataFrame(data)
    
    # 2. Setup Plot
    plt.figure(figsize=(10, 6))
    
    # 3. Plot the Curve
    plt.plot(df['angle'], df['variance_ratio'], 
             marker='o', linestyle='-', color='#2c3e50', linewidth=2, label='Observed Variance')
    
    # 4. Add "Isotropic Gravity" Baseline (Theoretical)
    # If structure were physical (Gravity), V would be constant at V(0)
    v_isotropic = df['variance_ratio'].iloc[0]
    plt.axhline(y=v_isotropic, color='#e74c3c', linestyle='--', alpha=0.7, label='Isotropic Expectation (Gravity)')
    
    # 5. Annotate Key Features
    
    # The "Lock" (0 degrees)
    plt.annotate('Context Lock\n(Max Variance)', 
                 xy=(0, df['variance_ratio'].iloc[0]), 
                 xytext=(10, df['variance_ratio'].iloc[0] + 2),
                 arrowprops=dict(facecolor='black', shrink=0.05))
                 
    # The "Decoherence" (Minima ~20 degrees)
    min_idx = df['variance_ratio'].iloc[0:10].idxmin() # Find min in first 50 deg
    min_val = df['variance_ratio'].iloc[min_idx]
    min_angle = df['angle'].iloc[min_idx]
    
    plt.annotate('Decoherence\n(Structure Collapses)', 
                 xy=(min_angle, min_val), 
                 xytext=(min_angle + 5, min_val - 4),
                 arrowprops=dict(facecolor='black', shrink=0.05))

    # The "Aliasing Spike" (45 degrees)
    spike_idx = df[df['angle'] == 45].index[0]
    spike_val = df['variance_ratio'].iloc[spike_idx]
    plt.annotate('Grid Aliasing\n(Moiré Effect)', 
                 xy=(45, spike_val), 
                 xytext=(50, spike_val),
                 arrowprops=dict(facecolor='gray', shrink=0.05))

    # 6. Styling
    plt.title('P6-B: Rotational Decoherence Curve (Synthetic Audit)', fontsize=14, fontweight='bold')
    plt.xlabel('Rotation Angle (degrees)', fontsize=12)
    plt.ylabel('Variance Ratio (V)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    
    # Shade the "Anisotropy Gap"
    plt.fill_between(df['angle'], df['variance_ratio'], v_isotropic, 
                     where=(df['variance_ratio'] < v_isotropic),
                     color='#e74c3c', alpha=0.1, label='Anisotropy Gap')

    # 7. Save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"[*] Plot saved to {output_path}")

if __name__ == "__main__":
    JSON_FILE = "data/processed/p6_rotation_results.json"
    IMG_FILE = "data/processed/p6_rotation_curve.png"
    
    if not os.path.exists(JSON_FILE):
        print(f"Error: {JSON_FILE} not found. Run p6_rotation_audit.py first.")
    else:
        plot_rotation_curve(JSON_FILE, IMG_FILE)
