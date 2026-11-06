# =============================================================================
# Constants.jl
#
# Módulo Constants: constantes físicas y rutas comunes a todo SIRIUS, en las
# unidades del proyecto (km, s, rad). Los parámetros del programa que no son
# constantes físicas están en Parameters.jl.
#
# Constantes:
#   MU_EARTH, R_EARTH, J2_EARTH  : campo gravitatorio terrestre (EGM2008, coherentes con Propagator)
#   OMEGA_EARTH                  : velocidad de rotación terrestre
#   JD_J2000                     : fecha juliana de la época J2000
#   WGS84_A, WGS84_F, WGS84_E2   : elipsoide WGS84 (posición de las estaciones)
#   R_SUN, D_SUN                 : radio del Sol y distancia media Tierra-Sol (eclipses)
#   DATA_DIR                     : carpeta data/ (kernels SPICE y catálogo Hipparcos)
# =============================================================================

module Constants

export MU_EARTH, R_EARTH, J2_EARTH, OMEGA_EARTH, JD_J2000, DATA_DIR
export WGS84_A, WGS84_F, WGS84_E2, R_SUN, D_SUN

# Campo gravitatorio terrestre (EGM2008); J2_EARTH = -C20 de la tabla de Propagator.load_coefs_n4
const MU_EARTH = 398600.4415       # km³/s²
const R_EARTH  = 6378.1363         # km, radio de referencia del modelo de gravedad
const J2_EARTH = 1.0826261739e-3   # adimensional, sin normalizar

const OMEGA_EARTH = 7.2921150e-5   # rad/s (rotación terrestre)
const JD_J2000 = 2451545.0         # 2000-01-01 12h TT

# Carpeta data/ de la raíz del repositorio (kernels SPICE y catálogo Hipparcos)
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

# Elipsoide WGS84: semieje mayor (km), achatamiento y excentricidad al cuadrado
const WGS84_A = 6378.137
const WGS84_F = 1/298.257223563
const WGS84_E2 = WGS84_F*(2 - WGS84_F)

# Sol, para los conos de sombra (Visibility)
const R_SUN = 696340.0             # km
const D_SUN = 149597870.7          # km

end
