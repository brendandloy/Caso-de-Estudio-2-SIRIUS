module Constants

export μ, R, ω_earth, jd_J2000, DATA_DIR, WGS84_A, WGS84_f, WGS84_E2
export R_sun, D_sun

const μ  = 398600.4415      # km³/s²
const R  = 6378.1363          # km
const ω_earth = 7.2921150e-5 # rad/s (rotación terrestre)
const jd_J2000 = 2451545.0 # 2000-01-01 12h TT

# Carpeta data/ de la raíz del repositorio (kernels SPICE y catálogo Hipparcos)
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

const WGS84_A = 6378.137
const WGS84_f = 1/298.257223563
const WGS84_E2 = WGS84_f*(2 - WGS84_f)

# Photometry/environment constants in project units
const R_sun = 696340.0         # km
const D_sun = 149597870.7      # km

end
