import streamlit as st
import backend as back
import matplotlib.pyplot as plt
import numpy as np

def main():
    # Configuración de la página
    st.set_page_config(
        page_title="KEPLER EXOSKY",
        page_icon="🌌",
        layout="wide"
    )


# Estilos CSS personalizados
st.markdown("""
<style>
    /* Degradado de fondo entre tres colores */
    .stApp {
        background: linear-gradient(to bottom, #4c6eb7, #ffffff, #000000);
        display: flex;
        flex-direction: column;
        align-items: center;
    }
    
    /* Estilos para el título EXOSKY */
    .titulo-exosky {
        font-size: 80px;
        color: #ffffff;
        font-family: 'Trebuchet MS', sans-serif;
        text-align: center;
        margin-bottom: 0;
    }
    
    /* Estilos para el subtítulo */
    .titulo-subtitulo {
        font-size: 30px;
        color: #ffffff;
        font-family: 'Trebuchet MS', sans-serif;
        text-align: center;
        margin-top: 0;
    }
    
    /* Estilos personalizados para el botón */
    .stButton > button {
        background-color: #4c6eb7;
        color: white;
        border: none;
        padding: 15px 32px;
        text-align: center;
        text-decoration: none;
        display: inline-block;
        font-size: 16px;
        margin: 20px 0;
        cursor: pointer;
        border-radius: 5px;
        transition: background-color 0.3s;
        font-weight: bold;
    }
    
    /* Estilo al pasar el mouse por encima del botón */
    .stButton > button:hover {
        background-color: #A9A9A9;
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
    }
    
    /* Estilo cuando se presiona el botón */
    .stButton > button:active {
        background-color: #3e8e41;
        transform: translateY(2px);
    }
    
    /* Centrar el botón */
    div.stButton {
        display: flex;
        justify-content: center;
    }

    /* Estilo para la imagen en la esquina superior derecha */
    .corner-image {
        position: absolute;
        top: 10px;
        right: -50px;  /* Cambiado a -50px para permitir ajuste manual */
        width: 100px;
        height: auto;
    }
</style>
""", unsafe_allow_html=True)

# Añadir imagen en la esquina superior derecha
st.markdown("""
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e5/NASA_logo.svg/1200px-NASA_logo.svg.png" alt="Corner Image" class="corner-image">
""", unsafe_allow_html=True)

# EXOSKY TITLE
st.markdown("<h1 class='titulo-exosky'>KEPLER EXOSKY</h1>", unsafe_allow_html=True)

# SUBTITLE SEARCH ENGINE
st.markdown("<h3 class='titulo-subtitulo'>SEARCH ENGINE</h3>", unsafe_allow_html=True)

# Estado del botón
if 'button_pressed' not in st.session_state:
    st.session_state.button_pressed = False

# Crear el texto del botón
button_text = "Travel to another ExoPlanet" if st.session_state.button_pressed else "Travel to a random ExoPlanet"

# Crear un botón de Streamlit (ahora centrado debido al CSS)
if st.button(button_text):
    st.session_state.button_pressed = True  # Cambiar el estado a presionado
    
    # Main execution
    stars_limit = 100000
    stars = back.query_gaia_exoplanets(stars_limit)
    
    # Ruta al archivo CSV de exoplanetas
    exoplanet_file_path = r'C:\Users\JG Rodrigo\Downloads\Codigo Fuente proyecto\gaiadata.csv'
    exoplanets = back.read_exoplanet_data(exoplanet_file_path)
    
    # Seleccionar un exoplaneta al azar
    random_exoplanet = np.random.choice(exoplanets)
    
    # Crear el mapa estelar desde la perspectiva del exoplaneta seleccionado
    back.create_2d_starmap(stars, random_exoplanet)