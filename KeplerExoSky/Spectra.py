import matplotlib.pyplot as plt
import numpy as np
from astroquery.gaia import Gaia
from astropy import units as u
from astropy.coordinates import SkyCoord
from gaiaxpy import calibrate
import gaiaxpy
from matplotlib.colors import LinearSegmentedColormap

Debug = True


def process_spectra(stars):
    if stars is None or len(stars) == 0:
        print("No stars data to process.")
        return None
    
    # Check if necessary columns exist
    required_columns = ['phot_bp_mean_mag', 'phot_rp_mean_mag', 'SOURCE_ID']
    for col in required_columns:
        if col not in stars.colnames:
            print(f"Error: Missing column {col} in stars data")
            return None
 
    try:
        coefficients = {
            'bp': np.array(stars['phot_bp_mean_mag']),
            'rp': np.array(stars['phot_rp_mean_mag'])
        }
        source_ids = stars['SOURCE_ID'] if 'SOURCE_ID' in stars.colnames else np.arange(len(stars))
        
        # Set the sampling parameter to ensure it's within the valid range
        sampling = np.geomspace(331, 1049, 64)  # Generate values from 330 to 1050 with step 10
        
        # Debugging: Print the values to check
        if Debug == True:
            print("Sampling:", sampling)
            print("Coefficients:", coefficients)
            print("Source IDs:", source_ids)
        
        # Remove the sampling argument from here
        xp_continuous, xp_samples = calibrate(coefficients, source_ids, sampling=sampling)
        return xp_continuous
    except KeyError as e:
        print(f"Error: Missing column in stars data: {e}")
        return None
    except Exception as e:
        print(f"Error processing spectra: {e}")
        return None
    