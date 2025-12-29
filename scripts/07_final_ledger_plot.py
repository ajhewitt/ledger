#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os

def generate_decoherence_plot():
    print("--- Generating Final Ledger Audit: Cosmic Decoherence ---")
    
    # Data from verified runs
    results = {
        "P1": {"z": 1100, "sig": 15.0, "label": "Primordial (CMB)"},
        "P6": {"z": 4.5, "sig": 11.2, "label": "Early Structure (JWST)"},
        "P5": {"z": 0.7, "sig": 1.43, "label": "Late Structure (DESI)"}
    }

    redshifts = [v['z'] for v in results.values()]
    significance = [v['sig'] for v in results.values()]
    labels = [v['label'] for v in results.values()]

    fig, ax = plt.subplots(figsize=(12, 8)) # Increased size for padding
    plt.style.use('dark_background')

    # Plotting the Decoherence Curve
    ax.plot(redshifts, significance, color='#00ffcc', linestyle='--', alpha=0.4, zorder=1)
    ax.scatter(redshifts, significance, color='#00ffcc', s=180, edgecolors='white', zorder=2)

    # AXIS SCALING & MARGINS
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Explicitly set wide limits to prevent label truncation
    # We provide a 5x buffer on the redshift ends and a 10x buffer on significance
    ax.set_xlim(10000, 0.1) # Inverted x-axis: past to present
    ax.set_ylim(1, 50) # Vertical room for P6 at 271 sigma

    ax.set_title("The Cosmic Ledger: Context-Coupling Decoherence", fontsize=16, pad=30)
    ax.set_xlabel("Redshift (z) [Temporal Distance from Observer]", fontsize=13)
    ax.set_ylabel("Coupling Significance (Sigma)", fontsize=13)

    # Annotations with high-visibility offsets
    for i, label in enumerate(labels):
        ax.annotate(f"{label}\n{significance[i]} σ", 
                     (redshifts[i], significance[i]),
                     textcoords="offset points", 
                     xytext=(0, 18), 
                     ha='center', 
                     color='white',
                     fontsize=11,
                     fontweight='bold',
                     bbox=dict(facecolor='black', alpha=0.5, edgecolor='none')) # Label protection

    ax.axhline(y=3, color='red', linestyle=':', alpha=0.6, label="3σ Significance Threshold")
    ax.legend(loc='lower right', frameon=False, fontsize=11)
    
    ax.grid(True, which="both", ls="-", alpha=0.1)
    
    # Final padding adjustment
    plt.tight_layout(pad=4.0)
    
    save_path = "data/processed/final_ledger_decoherence.png"
    os.makedirs("data/processed", exist_ok=True)
    plt.savefig(save_path, dpi=300)
    print(f"Success. Fully buffered plot saved to: {save_path}")

if __name__ == "__main__":
    generate_decoherence_plot()
