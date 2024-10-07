import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from astroquery.gaia import Gaia
from astropy import units as u
from astropy.coordinates import SkyCoord
import csv
import streamlit as st

# Configuración de la consulta a Gaia
def query_gaia_exoplanets(limit):
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
    job = Gaia.launch_job(query)
    return job.get_results()

# Calcular la temperatura de color
def calculate_color_temperature(bp_rp):
    # Relación aproximada de BP-RP con la temperatura (K)
    return 10400 * (1 / (0.92 * bp_rp + 1.7) + 1 / (0.92 * bp_rp + 0.62))

# Crear mapa estelar con colores reales

def read_exoplanet_data(file_path):
    exoplanets = []
    with open(file_path, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            exoplanets.append(row)
    return exoplanets

def kelvin_to_rgb(temperature):
    """Convert a temperature in Kelvin to an RGB color."""
    temperature = temperature / 100

    if temperature <= 66:
        r = 255
    else:
        r = 329.698727446 * ((temperature - 60) ** -0.1332047592)
    
    if temperature <= 66:
        g = 99.4708025861 * np.log(temperature) - 161.1195681661
    else:
        g = 288.1221695283 * ((temperature - 60) ** -0.0755148492)
    
    if temperature >= 66:
        b = 255
    elif temperature <= 19:
        b = 0
    else:
        b = 138.5177312231 * np.log(temperature - 10) - 305.0447927307

    return np.clip([r, g, b], 0, 255) / 255

def create_star_colormap():
    """Create a custom colormap for star temperatures."""
    temperatures = np.linspace(2000, 30000, 256)
    colors = [kelvin_to_rgb(temp) for temp in temperatures]
    return LinearSegmentedColormap.from_list("star_colors", colors)

def create_2d_starmap(stars, exoplanet):
    fig, ax = plt.subplots(figsize=(12, 8), facecolor='black')
    
    x = stars['ra']
    y = stars['dec']  # Invert the declination values

    sizes = 1 * 10**(-stars['phot_g_mean_mag']/5)  # Increased size factor
    temperatures = calculate_color_temperature(stars['bp_rp'])

    norm = plt.Normalize(vmin=2000, vmax=30000)  # Set fixed range for better color representation
    star_cmap = create_star_colormap()

    scatter = ax.scatter(x, y, s=sizes, c=temperatures, cmap=star_cmap, norm=norm, alpha=0.8)

    ax.set_facecolor('black')
    ax.set_xlim(0, 360)
    ax.set_ylim(90, -90)  # Reverse the y-axis limits
    
    ax.set_xlabel('Right Ascension (degrees)', color='white')
    ax.set_ylabel('Declination (degrees)', color='white')
    ax.tick_params(colors='white')

    ax.set_title(f'View from {exoplanet["pl_name"]}', color='red')

    cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('Temperature (K)', color='white')
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

    plt.tight_layout()

    # En lugar de plt.show(), usa st.pyplot
    st.pyplot(fig)



# Main execution
stars_limit = 100000
stars = query_gaia_exoplanets(stars_limit)

# Ruta al archivo CSV de exoplanetas
exoplanet_file_path = r'C:\Users\JG Rodrigo\Downloads\Codigo Fuente proyecto\gaiadata.csv'
exoplanets = read_exoplanet_data(exoplanet_file_path)

# Seleccionar un exoplaneta al azar
random_exoplanet = np.random.choice(exoplanets)

# Crear el mapa estelar desde la perspectiva del exoplaneta seleccionado
create_2d_starmap(stars, random_exoplanet)
