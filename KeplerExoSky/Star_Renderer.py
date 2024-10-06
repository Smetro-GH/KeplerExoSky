import matplotlib.pyplot as plt
import numpy as np
from astroquery.gaia import Gaia
from astropy import units as u
from astropy.coordinates import SkyCoord
from gaiaxpy import calibrate
import gaiaxpy
from matplotlib.colors import LinearSegmentedColormap
from Spectra import process_spectra


Debug = True

stars = 100

def query_gaia_dr3(limit):
    query = f"""
    SELECT TOP {limit}
        SOURCE_ID, ra, dec, bp_rp, phot_g_mean_mag,
        phot_bp_mean_mag, phot_rp_mean_mag
    FROM gaiadr3.gaia_source
    WHERE phot_g_mean_mag < 6
      AND bp_rp IS NOT NULL
      AND phot_bp_mean_mag IS NOT NULL
      AND phot_rp_mean_mag IS NOT NULL
    ORDER BY phot_g_mean_mag ASC
    """
    if Debug == True:
        try:
            job = Gaia.launch_job(query)
            results = job.get_results()
            if len(results) == 0:
                print("Query returned no results.")
                return None
            print("Available columns:", results.colnames)
            return results
        except Exception as e:
            print(f"Error executing Gaia query: {e}")
            return None
    

# Simple color temperature calculation
def calculate_color_temperature(spectra):
    # Simple color temperature calculation
    blue_flux = np.mean(spectra[:, 100:150], axis=1)  # ~400-450 nm
    red_flux = np.mean(spectra[:, 450:500], axis=1)   # ~750-800 nm
    return blue_flux / red_flux

def create_starmap(stars, spectra):
    fig, ax = plt.subplots(figsize=(15, 12), facecolor='black')
    ax.set_facecolor('black')
    
    # Create custom colormap
    colors = ['blue', 'white', 'yellow', 'orange', 'red']
    cmap = LinearSegmentedColormap.from_list('star_cmap', colors, N=100)
    
    coords = SkyCoord(ra=stars['ra']*u.degree, dec=stars['dec']*u.degree, frame='icrs')
    l = coords.galactic.l.degree
    b = coords.galactic.b.degree
    
    sizes = 50 * 10**(-stars['phot_g_mean_mag']/5)
    color_temp = calculate_color_temperature(spectra)
    
    scatter = ax.scatter(l, b, s=sizes, c=color_temp, cmap=cmap, 
                         norm=plt.Normalize(vmin=0.5, vmax=1.5),
                         edgecolors='none', alpha=0.8)
    
    ax.set_xlabel('Galactic Longitude', color='white')
    ax.set_ylabel('Galactic Latitude', color='white')
    ax.set_title(f'Gaia DR3 Star Map with Spectral Data ({len(stars)} Brightest Stars)', color='white')

    ax.set_xlim(360, 0)
    ax.set_ylim(-90, 90)

    cbar = plt.colorbar(scatter)
    cbar.set_label('Color Temperature (Blue/Red Flux Ratio)', color='white')
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

    plt.grid(color='gray', linestyle=':', alpha=0.3)
    plt.tight_layout()
    plt.show()

# Main execution
stars = query_gaia_dr3(stars)
if stars is not None:
    spectra = process_spectra(stars)
    if spectra is not None:
        create_starmap(stars, spectra)
        plt.show()