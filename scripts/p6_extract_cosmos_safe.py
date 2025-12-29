import os
import numpy as np
import pandas as pd
from astropy.io import fits

def safe_extract(fits_path, output_csv):
    print(f"[*] Opening {fits_path} using memory mapping...")
    
    try:
        # memmap=True is the key here; it doesn't load the file into RAM
        with fits.open(fits_path, memmap=True) as hdul:
            # COSMOS2020 data is in extension 1
            data = hdul[1].data
            
            print("[*] Accessing specific columns...")
            # We only pull the data for these specific headers
            # Using list of possible names in case of version differences
            ra_col = 'ALPHA_J2000' if 'ALPHA_J2000' in data.names else 'ra'
            dec_col = 'DELTA_J2000' if 'DELTA_J2000' in data.names else 'dec'
            z_col = 'lp_zPDF'
            
            # Extract just these columns as numpy arrays
            ra = np.array(data[ra_col], dtype=np.float64)
            dec = np.array(data[dec_col], dtype=np.float64)
            z = np.array(data[z_col], dtype=np.float32)
            
            print(f"[*] Filtering {len(ra)} sources for Structuring Phase (3 < z < 6)...")
            # Create a boolean mask to filter before making a DataFrame
            mask = (z > 3.0) & (z < 6.0)
            
            # Build the DataFrame only from the filtered subset
            df_subset = pd.DataFrame({
                'ra': ra[mask],
                'dec': dec[mask],
                'lp_zPDF': z[mask]
            })
            
            os.makedirs(os.path.dirname(output_csv), exist_ok=True)
            df_subset.to_csv(output_csv, index=False)
            print(f"[SUCCESS] Extracted {len(df_subset)} sources to {output_csv}")
            
    except Exception as e:
        print(f"[ERROR] Extraction failed: {e}")

if __name__ == "__main__":
    # Check if the file is .gz (compressed)
    # Note: Memory mapping doesn't work well on .gz files. 
    # If it's cosmos2020_classic.fits.gz, unzip it first!
    raw_path = "data/raw/jwst/cosmos2020_classic.fits" 
    output_path = "data/raw/COSMOS2020_subset.csv"
    
    if os.path.exists(raw_path):
        safe_extract(raw_path, output_path)
    else:
        print(f"[!] File not found: {raw_path}")
        print("[!] Tip: If your file ends in .gz, run 'gunzip data/raw/jwst/cosmos2020_classic.fits.gz' first.")
