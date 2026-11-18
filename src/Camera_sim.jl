# =============================================================================
# Camera_sim.jl
#
# Módulo CameraSim: simulador de la cámara de cada observador (CE3 de I. Ferrando,
# ampliado con la orientación de la cámara). En cada instante selecciona las estrellas
# (Hipparcos) y los RSOs dentro del FOV, los proyecta sobre el sensor, renderiza la
# imagen con PSF gaussiana y ruido, y guarda medidas LOS y posiciones en píxeles.
#
# Orientación de la cámara: llega como C_ci (ECI -> cámara) por instante, C_ci = C_cb · C_bi
# (ver Attitude/Mount.jl). Ejes cámara dextrógiros: x_c derecha de la imagen, y_c abajo,
# z_c = boresight.
#
# Intrínsecos (CameraIntrinsics, bloque camera_config): resolución, focales f_x, f_y (px),
# punto principal (u0, v0), skew α y distorsión de Brown-Conrady (k1, k2, k3 radiales y
# p1, p2 tangenciales; Brown, 1966 y 1971). Proyección pinhole: coordenadas normalizadas
# x' = x_c/z_c, y' = y_c/z_c, distorsión (x'', y'') y matriz K (Memoria Nagoya, ec. 2.5):
# u = f_x x'' + α y'' + u0, v = f_y y'' + v0. Píxeles
# desde 1 (el píxel k ocupa [k - ½, k + ½]). Una fuente se ve si su píxel cae dentro del
# sensor (FOV rectangular); antes se descartan las que están fuera de un cono que contiene
# todo el sensor (prefiltro rápido para el catálogo de estrellas).
#
# Exclusión solar: si el Sol está a menos de camera_config.sun_exclusion_deg del eje óptico,
# su luz entra en la óptica y satura la imagen; en ese instante no se genera imagen ni
# medidas. Cada instante queda en frame_status.csv con el ángulo Sol-boresight.
#
# Tierra: las estrellas que tapa la Tierra (esfera de radio R⊕ + altura de la atmósfera)
# no se pintan; los RSOs tapados ya los quita la fotometría (losECI). El disco de la Tierra
# queda negro, también de día. frame_status.csv guarda si la Tierra está en la imagen y qué
# fracción ocupa, y cuánta de ella está de día, medidas con una malla de rayos.
#
# Funciones principales: CameraIntrinsics, camera_intrinsics, distort, undistort,
# sensor_fov_deg, prefilter_cone_deg, sun_exclusion_deg, load_hipparcos, prefilter_box, filter_cone, radec_to_pixels!,
# prepare_render_dict, in_sensor, sources_in_view, gaussian_kernel, render_psf_gaussian!,
# ruido (add_*), build_debris_catalog, err_meas!, star_frame, camera_sim! (punto de entrada)
# y sus pasos: setup_star_field, frame_star_field!, start_frame!, paint_star_field!,
# sky_background, frame_conditions, earth_coverage (con is_earth_near_view, is_earth_on_border,
# pixel_ray_eci y grid_centers), hide_stars_behind_earth, debris_in_view, log_debris!,
# finish_frame!, frame_name y star_catalog_rows.
#
# Dependencias: Types, Constants (DATA_DIR), Parameters (FOV_PREFILTER_MARGIN,
# EARTH_COVERAGE_GRID_SIZE), JSONUtils (lectura de camera_config) y Visibility (ángulo
# Sol-boresight y geometría de la Tierra), que deben estar cargados antes.
# =============================================================================

module CameraSim

export camera_sim!, CameraIntrinsics, camera_intrinsics, sensor_fov_deg, sun_exclusion_deg, distort, undistort

using FITSIO
using DataFrames
using Images, ImageIO
using Random
using Statistics
using StaticArrays
using LinearAlgebra
using CSV
using Dates

using ..Types
using ..Constants: DATA_DIR     # carpeta data/ (catálogo Hipparcos)
using ..Parameters: FOV_PREFILTER_MARGIN, EARTH_COVERAGE_GRID_SIZE   # parámetros maestros de la cámara
using ..JSONUtils: todict, getnum, check_keys   # lectura y validación de camera_config
using ..Visibility: sun_boresight_angle_deg, earth_angular_radius_deg, earth_ray_hit, is_earth_occulted

# -----------------------------------------------------------------------------
# Intrínsecos de la cámara
# -----------------------------------------------------------------------------

"""
    CameraIntrinsics(resolution_x, resolution_y, fx_px, fy_px, u0_px, v0_px, skew_px[, k1, k2, k3, p1, p2])

Parámetros intrínsecos de la cámara pinhole (Memoria Nagoya, ec. 2.5): resolución (px),
focales efectivas en píxeles por radián, punto principal (px, desde 1), skew α (px) y
coeficientes de distorsión de Brown-Conrady, adimensionales (por defecto 0, sin distorsión).
Se construye y valida con `camera_intrinsics`.
"""
struct CameraIntrinsics
    resolution_x::Int
    resolution_y::Int
    fx_px::Float64
    fy_px::Float64
    u0_px::Float64
    v0_px::Float64
    skew_px::Float64
    k1::Float64
    k2::Float64
    k3::Float64
    p1::Float64
    p2::Float64
end

# Sin distorsión
CameraIntrinsics(nx, ny, fx, fy, u0, v0, skew) = CameraIntrinsics(nx, ny, fx, fy, u0, v0, skew, 0.0, 0.0, 0.0, 0.0, 0.0)

# Coeficientes del bloque "distortion" (todos opcionales)
const DISTORTION_KEYS = ("k1", "k2", "k3", "p1", "p2")

"""
    camera_intrinsics(camera_config) -> CameraIntrinsics

Lee y valida los intrínsecos del bloque `camera_config`. Obligatorios: `resolution_x`,
`resolution_y`, `fx_px`, `fy_px`. Opcionales: `u0_px`, `v0_px` (por defecto, o si valen
`null`, el centro del sensor, (N + 1)/2), `skew_px` (por defecto 0) y el bloque
`distortion` con `k1`, `k2`, `k3`, `p1`, `p2` (cada uno por defecto 0).

Por qué existe: los intrínsecos se validan una vez, al leer el JSON, y la cámara trabaja
con un tipo de campos fijos.

# Argumentos
- `camera_config`: bloque del JSON (JSON3 o `Dict`).

# Devuelve
- `CameraIntrinsics`.

# Excepciones
- `ErrorException`: si falta una clave obligatoria, la resolución no es un entero positivo,
  una focal no es positiva o el bloque `distortion` tiene una clave desconocida.

# Ejemplo
```julia
cam = camera_intrinsics(Dict("resolution_x" => 1024, "resolution_y" => 1024,
                             "fx_px" => 2132.63, "fy_px" => 2132.63))   # u0 = v0 = 512.5
```
"""
function camera_intrinsics(camera_config)
    b = todict(camera_config)
    who = "camera_config"
    nx, ny = (getnum(b, k, who) for k in ("resolution_x", "resolution_y"))
    (isinteger(nx) && isinteger(ny) && nx > 0 && ny > 0) ||
        error("$who: 'resolution_x' and 'resolution_y' must be positive integers.")
    fx, fy = getnum(b, "fx_px", who), getnum(b, "fy_px", who)
    (fx > 0 && fy > 0) || error("$who: 'fx_px' and 'fy_px' must be > 0 (got $fx, $fy).")

    # Punto principal: por defecto, el centro del sensor con píxeles desde 1
    u0 = get(b, "u0_px", nothing) === nothing ? (nx + 1) / 2 : getnum(b, "u0_px", who)
    v0 = get(b, "v0_px", nothing) === nothing ? (ny + 1) / 2 : getnum(b, "v0_px", who)
    skew = haskey(b, "skew_px") ? getnum(b, "skew_px", who) : 0.0
    k1, k2, k3, p1, p2 = distortion_coeffs(get(b, "distortion", Dict{String,Any}()), "$who.distortion")
    return CameraIntrinsics(Int(nx), Int(ny), fx, fy, u0, v0, skew, k1, k2, k3, p1, p2)
end

"""
    sun_exclusion_deg(camera_config) -> Float64

Lee y valida `sun_exclusion_deg` de `camera_config`: semiángulo (°) del cono alrededor del
eje óptico dentro del cual el Sol satura la imagen. Lo fija el diseño del parasol (*baffle*);
valores típicos de sensores estelares, 30–45°.

# Argumentos
- `camera_config`: bloque `camera_config` del JSON (Dict o JSON3.Object).

# Devuelve
- `Float64` en [0, 180).

# Excepciones
- `ErrorException`: si falta la clave, no es un número o está fuera de [0, 180).
"""
function sun_exclusion_deg(camera_config)
    who = "camera_config"
    θ = getnum(todict(camera_config), "sun_exclusion_deg", who)
    0 <= θ < 180 || error("$who: 'sun_exclusion_deg' must be in [0, 180) (got $θ).")
    return θ
end

# Coeficientes (k1, k2, k3, p1, p2) del bloque "distortion"; los que faltan valen 0
function distortion_coeffs(d, who::AbstractString)
    d isa AbstractDict || error("$who: must be a block {k1, k2, k3, p1, p2}.")
    check_keys(d, DISTORTION_KEYS; who = who)
    return Tuple(haskey(d, k) ? getnum(d, k, who) : 0.0 for k in DISTORTION_KEYS)
end

"""
    distort(cam, x, y) -> (x_d, y_d)

Aplica la distorsión de Brown-Conrady a unas coordenadas normalizadas (x, y) = (x_c/z_c, y_c/z_c)
(Brown, 1966 y 1971; modelo de OpenCV y de `cameraIntrinsics` de MATLAB):

    r² = x² + y²,   R = 1 + k1 r² + k2 r⁴ + k3 r⁶
    x_d = x R + 2 p1 x y + p2 (r² + 2x²)
    y_d = y R + p1 (r² + 2y²) + 2 p2 x y

Los términos radiales (k) estiran o encogen la imagen con la distancia al centro; los
tangenciales (p), la lente descentrada respecto al sensor.

# Argumentos
- `cam::CameraIntrinsics`; `x`, `y`: coordenadas normalizadas.

# Devuelve
- `(x_d, y_d)`, coordenadas normalizadas con distorsión.

# Excepciones
- Ninguna.

# Ejemplo
```julia
distort(CameraIntrinsics(1024, 1024, 1e3, 1e3, 512.5, 512.5, 0.0, 0.1, 0, 0, 0, 0), 0.2, 0.0)   # (0.2008, 0.0)
```
"""
function distort(cam::CameraIntrinsics, x::Real, y::Real)
    r2 = x^2 + y^2
    R = 1 + r2 * (cam.k1 + r2 * (cam.k2 + r2 * cam.k3))
    xd = x * R + 2cam.p1 * x * y + cam.p2 * (r2 + 2x^2)
    yd = y * R + cam.p1 * (r2 + 2y^2) + 2cam.p2 * x * y
    return xd, yd
end

"""
    undistort(cam, x_d, y_d; iters=50) -> (x, y)

Inversa de `distort` por punto fijo (como `undistortPoints` de OpenCV): x = (x_d − Δ_t(x, y))/R(x, y),
iterando desde (x_d, y_d). Converge mientras la distorsión sea moderada.

Por qué existe: para saber qué dirección cae en una esquina del sensor (`prefilter_cone_deg`).

# Argumentos
- `cam::CameraIntrinsics`; `x_d`, `y_d`: coordenadas normalizadas con distorsión.
- `iters::Int`: número de iteraciones.

# Devuelve
- `(x, y)`, coordenadas normalizadas sin distorsión.

# Excepciones
- Ninguna.
"""
function undistort(cam::CameraIntrinsics, xd::Real, yd::Real; iters::Int = 50)
    x, y = Float64(xd), Float64(yd)
    for _ in 1:iters
        r2 = x^2 + y^2
        R = 1 + r2 * (cam.k1 + r2 * (cam.k2 + r2 * cam.k3))
        dx = 2cam.p1 * x * y + cam.p2 * (r2 + 2x^2)
        dy = cam.p1 * (r2 + 2y^2) + 2cam.p2 * x * y
        x, y = (xd - dx) / R, (yd - dy) / R
    end
    return x, y
end

# Distorsión radial dentro de su dominio de validez: r_d = r R(r) creciente con r,
# d(r R)/dr = 1 + 3k1 r² + 5k2 r⁴ + 7k3 r⁶ > 0. Fuera, el modelo "dobla" la imagen
radial_distortion_ok(cam::CameraIntrinsics, r2::Real) = 1 + r2 * (3cam.k1 + r2 * (5cam.k2 + r2 * 7cam.k3)) > 0

# Coordenadas normalizadas (x', y') = (x_c/z_c, y_c/z_c) de un píxel (u, v): inversa de K y de la distorsión
function pixel_to_normalized(cam::CameraIntrinsics, u::Real, v::Real)
    yd = (v - cam.v0_px) / cam.fy_px
    xd = (u - cam.u0_px - cam.skew_px * yd) / cam.fx_px
    return undistort(cam, xd, yd)
end

"""
    sensor_fov_deg(cam) -> (fov_x, fov_y)

Campo de visión (°) horizontal y vertical del sensor: ángulo entre los bordes izquierdo y
derecho (superior e inferior) vistos desde la cámara, `atan((u0 - ½)/f_x) + atan((N + ½ - u0)/f_x)`.
Con el punto principal centrado es `2 atan(N/(2 f_x))`. No es una entrada: se informa.

# Argumentos
- `cam::CameraIntrinsics`.

# Devuelve
- `(fov_x, fov_y)` en grados (sin skew ni distorsión: es el FOV nominal de K).

# Excepciones
- Ninguna.

# Ejemplo
```julia
sensor_fov_deg(CameraIntrinsics(1024, 1024, 2132.63, 2132.63, 512.5, 512.5, 0.0))   # ≈ (27.0, 27.0)
```
"""
function sensor_fov_deg(cam::CameraIntrinsics)
    edges(c, n, f) = atand((c - 0.5) / f) + atand((n + 0.5 - c) / f)
    return edges(cam.u0_px, cam.resolution_x, cam.fx_px), edges(cam.v0_px, cam.resolution_y, cam.fy_px)
end

"""
    prefilter_cone_deg(cam) -> Float64

Apertura total (°) del cono de prefiltro alrededor del boresight: el doble del ángulo a la
esquina del sensor más alejada, multiplicado por `FOV_PREFILTER_MARGIN` y limitado a 179.8°.

La esquina se pasa a dirección deshaciendo K y la distorsión (`pixel_to_normalized`).

Por qué existe: filtrar el catálogo con un cono es rápido; el margen cubre los errores de
la inversa de la distorsión y evita perder fuentes cerca del borde.

# Argumentos
- `cam::CameraIntrinsics`.

# Devuelve
- `Float64` en grados.

# Excepciones
- Ninguna.
"""
function prefilter_cone_deg(cam::CameraIntrinsics)
    corners = ((u, v) for u in (0.5, cam.resolution_x + 0.5), v in (0.5, cam.resolution_y + 0.5))
    θ_corner = maximum(atand(hypot(pixel_to_normalized(cam, u, v)...)) for (u, v) in corners)
    return min(2 * FOV_PREFILTER_MARGIN * θ_corner, 179.8)
end

# Normalice RA to [0,360)
@inline function wrap360(x)
    return mod(x, 360.0)
end

@inline function mag_to_flux(m::Real, mag0::Real, flux0::Real)
    return flux0 * 10.0^(-0.4 * (m - mag0))
end

function load_hipparcos(path::AbstractString, mag_lim::Real)

    f = FITS(path)
    hdu = f[2] 

    # --- Coordenadas ---
    ra   = read(hdu, "RAdeg")   # grados
    dec  = read(hdu, "DEdeg")   # grados

    # --- Magnitud ---
    vmag = read(hdu, "Vmag")    # Johnson V
    
    close(f)

    # --- Filtro por magnitud ---
    mask = (vmag .<= mag_lim) .& isfinite.(Float64.(vmag))

    return DataFrame(
        ra  = ra[mask],
        dec = dec[mask],
        mag = vmag[mask]
    )
end


function vec_to_radec(v::AbstractVector{<:Real})
    x, y, z = v

    dec = rad2deg(asin(z))
    ra  = rad2deg(atan(y, x))   # atan2(y, x)

    return wrap360(ra), dec
end


function radec_to_vec(ra_deg::Real, dec_deg::Real)
    α = deg2rad(Float64(ra_deg))
    δ = deg2rad(Float64(dec_deg))
    cδ = cos(δ)

    return SVector{3,Float64}(
        cδ * cos(α),
        cδ * sin(α),
        sin(δ)
    )
end

function radec_box_square(ra0::Real, dec0::Real, fov_deg::Real)

    half = fov_deg / 2

    dec_min_raw = dec0 - half
    dec_max_raw = dec0 + half

    dec_min = max(dec_min_raw, -90.0)
    dec_max = min(dec_max_raw,  90.0)

    # If FOV contains pole
    is_pole_in_fov = (dec_max_raw >= 90.0) || (dec_min_raw <= -90.0)

    c = cos(deg2rad(dec0))
    dra = is_pole_in_fov ? 180.0 : (half / max(c, 1e-12))  # 180 deg makes RA full-range

    ra_min = wrap360(ra0 - dra)
    ra_max = wrap360(ra0 + dra)

    # If FOV contains RA = 0
    is_ra_wrapped = ra_min > ra_max 

    # If RA wraps, use full RA range to keep the output simple
    ra_min = (is_pole_in_fov || is_ra_wrapped) ? 0.0 : ra_min
    ra_max = (is_pole_in_fov || is_ra_wrapped) ? 360.0 : ra_max

    return (ra_min = ra_min, ra_max = ra_max, 
        dec_min = dec_min, dec_max = dec_max)
end


"""
    prefilter_box(stars, boresight_vec, cone_deg) -> DataFrame

Primer filtro, muy rápido: fuentes con declinación a menos de `cone_deg/2` de la del
boresight (banda en Dec; el filtro en RA está desactivado).

# Argumentos
- `stars::DataFrame`: fuentes con `ra`, `dec` (°).
- `boresight_vec`: boresight unitario en ECI.
- `cone_deg::Real`: apertura total del cono de prefiltro (°).

# Devuelve
- Las filas de `stars` dentro de la banda.

# Excepciones
- Ninguna.
"""
function prefilter_box(stars, boresight_vec, cone_deg::Real)

    # Convert camera pointing to RA/Dec center
    ra0, dec0 = vec_to_radec(boresight_vec)

    # Build bounding box in RA/Dec
    box = radec_box_square(ra0, dec0, cone_deg)

    # Apply declination bounds
    mask = (stars.dec .>= box.dec_min) .& (stars.dec .<= box.dec_max)

    # Apply RA bounds (may be full range if pole/wrap was detected)
    # mask .&= (stars.ra .>= box.ra_min) .& (stars.ra .<= box.ra_max)

    return stars[mask, :]
end

"""
    filter_cone(stars, boresight_vec, cone_deg) -> DataFrame

Fuentes a menos de `cone_deg/2` del boresight (separación angular).

# Argumentos
- `stars::DataFrame`: fuentes con `ra`, `dec` (°).
- `boresight_vec`: boresight unitario en ECI.
- `cone_deg::Real`: apertura total del cono (°).

# Devuelve
- Las filas de `stars` dentro del cono.

# Excepciones
- Ninguna.
"""
function filter_cone(stars, boresight_vec, cone_deg::Real)

    bx, by, bz = boresight_vec
    cos_fov = cos(deg2rad(cone_deg / 2))

    is_in_cone = BitVector(undef, nrow(stars))

    for i in 1:nrow(stars)
        # Convert RA/Dec to a unit vector (ICRS)
        α = deg2rad(stars.ra[i])
        δ = deg2rad(stars.dec[i])
        cδ = cos(δ)

        sx = cδ * cos(α)
        sy = cδ * sin(α)
        sz = sin(δ)

        # Keep star if angular separation <= radius
        is_in_cone[i] = (sx*bx + sy*by + sz*bz) >= cos_fov
    end

    return stars[is_in_cone, :]
end


# Píxel que se asigna a una fuente detrás de la cámara (z_c <= 0): queda fuera de la imagen
const OFF_IMAGE = typemin(Int) ÷ 2

"""
    radec_to_pixels!(x, y, ra, dec, C_ci, cam) -> Nothing

Proyecta direcciones (RA, Dec) sobre la imagen con el modelo pinhole y los intrínsecos de
la cámara: se pasa cada dirección a ejes cámara y se divide por z_c (coordenadas
normalizadas), se aplica la distorsión (`distort`) y después la matriz K (Memoria Nagoya,
ec. 2.5):

  ρ_c = C_ci · û_i,   x' = x_c/z_c,   y' = y_c/z_c,   (x'', y'') = distort(x', y'),
  u = f_x x'' + α y'' + u0,   v = f_y y'' + v0.

Ejes cámara dextrógiros (x_c derecha, y_c abajo, z_c = boresight); píxeles desde 1. Es la
proyección gnomónica (TAN) del CE3 de I. Ferrando (ec. 2.15–2.19), con focales y punto
principal propios.

# Argumentos
- `x::Vector{Int}`, `y::Vector{Int}`: salida, columna y fila (1 = arriba a la izquierda).
- `ra::Vector{Float64}`, `dec::Vector{Float64}`: direcciones en ECI/ICRS (°).
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `cam::CameraIntrinsics`: intrínsecos.

# Devuelve
- `nothing` (rellena `x` e `y`, redondeados al píxel). Las fuentes detrás de la cámara, o
  fuera del dominio de validez de la distorsión radial (donde r(1 + k1 r² + ...) deja de
  crecer y la imagen se doblaría), reciben `OFF_IMAGE`.

# Excepciones
- `AssertionError`: si los vectores no tienen la misma longitud.

# Ejemplo
```julia
x = zeros(Int, 1); y = zeros(Int, 1)
radec_to_pixels!(x, y, [0.0], [0.0], C_ci, cam)    # boresight -> (u0, v0) redondeado
```
"""
function radec_to_pixels!(x::Vector{Int}, y::Vector{Int},
                          ra::Vector{Float64}, dec::Vector{Float64},
                          C_ci::AbstractMatrix, cam::CameraIntrinsics)

    @assert length(x) == length(y) == length(ra) == length(dec)

    C = SMatrix{3,3,Float64,9}(C_ci)

    @inbounds for i in eachindex(ra)
        # Dirección en ejes cámara
        ρ = C * radec_to_vec(ra[i], dec[i])

        # Detrás del plano imagen (no ocurre si se ha filtrado con el cono)
        if ρ[3] <= 0.0
            x[i] = OFF_IMAGE
            y[i] = OFF_IMAGE
            continue
        end

        # Coordenadas normalizadas; fuera del dominio de la distorsión no hay imagen
        xn, yn = ρ[1] / ρ[3], ρ[2] / ρ[3]
        if !radial_distortion_ok(cam, xn^2 + yn^2)
            x[i] = OFF_IMAGE
            y[i] = OFF_IMAGE
            continue
        end

        # Distorsión y matriz K: x hacia la derecha y y hacia abajo
        xd, yd = distort(cam, xn, yn)
        x[i] = round(Int, cam.fx_px * xd + cam.skew_px * yd + cam.u0_px)
        y[i] = round(Int, cam.fy_px * yd + cam.v0_px)
    end

    return nothing
end


"""
    prepare_render_dict(stars, C_ci, cam; mag0=8.0, flux0=20_000.0) -> NamedTuple

Calcula los píxeles y el flujo relativo de una lista de fuentes (estrellas o RSOs)
para una orientación de la cámara.

# Argumentos
- `stars::DataFrame`: fuentes con columnas `ra`, `dec` (°) y `mag`.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `cam::CameraIntrinsics`: intrínsecos.
- `mag0`, `flux0`: magnitud y flujo de referencia de la escala de flujo.

# Devuelve
- `(x, y, flux)`: columnas, filas y flujos relativos.

# Excepciones
- Ninguna.
"""
function prepare_render_dict(stars::DataFrame, C_ci::AbstractMatrix, cam::CameraIntrinsics;
                             mag0::Real = 8.0, flux0::Real = 20_000.0)

    n = nrow(stars)

    # Force dense, type-stable arrays (fast in hot loops)
    ra  = Vector{Float64}(stars.ra)
    dec = Vector{Float64}(stars.dec)
    mag = Vector{Float64}(stars.mag)

    x = Vector{Int}(undef, n)
    y = Vector{Int}(undef, n)

    # Project to pixels with the camera orientation
    radec_to_pixels!(x, y, ra, dec, C_ci, cam)

    # Magnitude -> flux (relative scaling)
    flux = similar(x, Float64)
    @inbounds for i in eachindex(mag)
        flux[i] = mag_to_flux(mag[i], mag0, flux0)
    end

    return (x = x, y = y, flux = flux)
end


# Píxel (x, y) dentro del sensor (píxeles desde 1); OFF_IMAGE queda fuera
in_sensor(x::Integer, y::Integer, cam::CameraIntrinsics) = 1 <= x <= cam.resolution_x && 1 <= y <= cam.resolution_y

"""
    sources_in_view(sources, C_ci, config) -> (DataFrame, NamedTuple)

Fuentes (estrellas o RSOs) que caen en el sensor con una orientación de la cámara, y sus
píxeles y flujos: cono de prefiltro (`prefilter_cone_deg`), proyección con los intrínsecos
y selección de las que caen dentro del sensor (FOV rectangular).

Por qué existe: es la misma selección para estrellas y RSOs.

# Argumentos
- `sources::DataFrame`: fuentes con `ra`, `dec` (°) y `mag`.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `config`: `camera_config` (usa `intrinsics`, `mag0` y `flux0`).

# Devuelve
- `(sources_in, render)`: filas de `sources` en el sensor y su `(x, y, flux)`.

# Excepciones
- Ninguna.
"""
function sources_in_view(sources::DataFrame, C_ci::AbstractMatrix, config)
    cam = config.intrinsics
    b = boresight_eci(C_ci)
    cone = prefilter_cone_deg(cam)
    near = filter_cone(prefilter_box(sources, b, cone), b, cone)
    render = prepare_render_dict(near, C_ci, cam; mag0 = config.mag0, flux0 = config.flux0)
    keep = in_sensor.(render.x, render.y, Ref(cam))
    return near[keep, :], (x = render.x[keep], y = render.y[keep], flux = render.flux[keep])
end

"""
    gaussian_kernel(fwhm_x_pix, fwhm_y_pix; k_sigma=4.0) -> (K, rx, ry)

Núcleo gaussiano normalizado (suma 1) de la PSF, con anchura propia en cada eje:
`exp(-Δx²/2σ_x² - Δy²/2σ_y²)`, con σ = FWHM / (2√(2 ln 2)). Se corta a `k_sigma` σ.

Por qué existe: con f_x ≠ f_y (píxeles no cuadrados) una PSF circular en el cielo es
elíptica en el sensor; con f_x = f_y es circular.

# Argumentos
- `fwhm_x_pix`, `fwhm_y_pix`: FWHM en columnas y en filas (px).
- `k_sigma::Real`: radio del núcleo en σ.

# Devuelve
- `K::Matrix{Float64}` de (2ry + 1) × (2rx + 1) y los semianchos `rx`, `ry` (px).

# Excepciones
- Ninguna.
"""
function gaussian_kernel(fwhm_x_pix::Real, fwhm_y_pix::Real; k_sigma::Real = 4.0)

    σx = Float64(fwhm_x_pix) / 2.3548200450309493
    σy = Float64(fwhm_y_pix) / 2.3548200450309493
    rx, ry = ceil(Int, k_sigma * σx), ceil(Int, k_sigma * σy)

    K = Matrix{Float64}(undef, 2ry + 1, 2rx + 1)
    @inbounds for j in 1:(2ry + 1), i in 1:(2rx + 1)
        dx = Float64(i - (rx + 1))
        dy = Float64(j - (ry + 1))
        K[j, i] = exp(-dx^2 / (2σx^2) - dy^2 / (2σy^2))
    end

    # Normalize once so total kernel sum = 1
    K ./= sum(K)

    return K, rx, ry
end


"""
    render_psf_gaussian!(img, render, cam, fwhm_deg; k_sigma=4.0) -> Nothing

Suma a `img` la PSF gaussiana de cada fuente de `render`. La FWHM angular se pasa a
píxeles con la escala del centro de la imagen en cada eje: `f_x · fwhm` y `f_y · fwhm`
(fwhm en rad); con píxeles cuadrados la mancha es circular. El skew no se tiene en cuenta
en la PSF (α ≈ 0 en sensores reales).

# Argumentos
- `img::Matrix{Float64}`: imagen (ny × nx); se modifica.
- `render::NamedTuple`: `(x, y, flux)` de las fuentes (`sources_in_view`).
- `cam::CameraIntrinsics`: intrínsecos.
- `fwhm_deg::Real`: FWHM de la PSF (°).

# Devuelve
- `nothing`.

# Excepciones
- Ninguna.
"""
function render_psf_gaussian!(img::Matrix{Float64}, render::NamedTuple, cam::CameraIntrinsics,
                              fwhm_deg::Real; k_sigma::Real = 4.0)

    x = render.x
    y = render.y
    flux = render.flux

    ny, nx = size(img)

    fwhm_rad = deg2rad(Float64(fwhm_deg))
    K, rx, ry = gaussian_kernel(cam.fx_px * fwhm_rad, cam.fy_px * fwhm_rad; k_sigma = k_sigma)

    @inbounds for k in eachindex(x)

        ix0 = x[k]
        iy0 = y[k]

        # Check if kernel is in image
        if ix0 + rx < 1 || ix0 - rx > nx || iy0 + ry < 1 || iy0 - ry > ny
            continue
        end

        F = flux[k]

        for kj in 1:(2ry + 1)
            j = iy0 + (kj - (ry + 1))
            if j < 1 || j > ny
                continue
            end
            for ki in 1:(2rx + 1)
                i = ix0 + (ki - (rx + 1))
                if i < 1 || i > nx
                    continue
                end
                img[j, i] += F * K[kj, ki]
            end
        end
    end

    return nothing
end


function add_background!(img::Matrix{Float64}, level::Real;
                         grad_x::Real = 0.0, grad_y::Real = 0.0)

    ny, nx = size(img)
    b0 = Float64(level)
    gx = Float64(grad_x)
    gy = Float64(grad_y)

    @inbounds for j in 1:ny
        row_off = gy * (j - 1)
        for i in 1:nx
            img[j, i] += b0 + gx * (i - 1) + row_off
        end
    end

    return nothing
end


@inline function poisson_knuth(rng::AbstractRNG, λ::Float64)
    # Knuth algorithm; efficient for small λ
    L = exp(-λ)
    k = 0
    p = 1.0
    while p > L
        k += 1
        p *= rand(rng)
    end
    return k - 1
end


@inline function poisson_sample(rng::AbstractRNG, λ::Float64)
    # Hybrid: Knuth for small λ, Normal approx for large λ
    if λ <= 0.0
        return 0
    elseif λ < 30.0
        return poisson_knuth(rng, λ)
    else
        # Normal approximation; clamp to nonnegative integer
        x = λ + sqrt(λ) * randn(rng)
        return Int(floor(x < 0.0 ? 0.0 : x + 0.5))
    end
end


function add_poisson_noise!(img::Matrix{Float64}; rng::AbstractRNG = Random.default_rng())

    @inbounds for idx in eachindex(img)
        λ = img[idx]
        img[idx] = Float64(poisson_sample(rng, λ))
    end

    return nothing
end


function add_read_noise!(img::Matrix{Float64}, sigma::Real;
                         rng::AbstractRNG = Random.default_rng(),
                         clamp_to_zero::Bool = false)

    σ = Float64(sigma)

    @inbounds for idx in eachindex(img)
        v = img[idx] + σ * randn(rng)
        img[idx] = clamp_to_zero ? (v < 0.0 ? 0.0 : v) : v
    end

    return nothing
end


function build_debris_catalog(obj_names::Vector{Symbol}, obj_trajs::Vector{TrajectoryData},
    photometry_data::AbstractDict{Tuple{Symbol,Symbol},<:NamedTuple},
    obs_name::Symbol, r_obs::SVector{3,Float64}, i::Int)

    obj_name_col = Symbol[]
    ra_col  = Float64[]
    dec_col = Float64[]
    mag_col = Float64[]

    @inbounds for (obj_name, obj_traj) in zip(obj_names, obj_trajs)
        phot = photometry_data[(obs_name, obj_name)]
        phot.is_visible[i] || continue

        m_obj = phot.Vmag[i]
        isfinite(m_obj) || continue

        rel = obj_traj.r[i] - r_obs
        d = norm(rel)
        d > 0.0 || continue
        
        α, δ = vec_to_radec(rel / d)

        push!(obj_name_col, obj_name)
        push!(ra_col, α)
        push!(dec_col, δ)
        push!(mag_col, m_obj)
    end

    return DataFrame(obj_name = obj_name_col, ra = ra_col, dec = dec_col, mag = mag_col)
end


function perturb_los(u::SVector{3,Float64}, σ_rad::Real)
    σ = Float64(σ_rad)

    if σ <= 0.0
        return u
    end

    # Build an orthonormal basis in the plane normal to the LOS.
    ref = abs(u[3]) < 0.9 ? SVector{3,Float64}(0.0, 0.0, 1.0) : SVector{3,Float64}(1.0, 0.0, 0.0)

    e1 = normalize(cross(u, ref))
    e2 = normalize(cross(u, e1))

    du = σ * randn() * e1 + σ * randn() * e2

    return normalize(u + du)
end


function err_meas!(obs_measurements::Dict{Symbol,Vector{LOSMeasurement}}, debris::DataFrame, t::Float64, los_sigma_arcsec::Real, has_meas_error::Bool)

    σ_los = deg2rad(Float64(los_sigma_arcsec) / 3600.0)

    @inbounds for k in 1:nrow(debris)

        obj_name = debris.obj_name[k]

        u_true = radec_to_vec(debris.ra[k], debris.dec[k])
        u_meas = has_meas_error ? perturb_los(u_true, σ_los) : SVector{3,Float64}(NaN, NaN, NaN)

        obj_measurements = get!(obs_measurements, obj_name) do
            LOSMeasurement[]
        end

        push!(obj_measurements,
            LOSMeasurement(
                t,
                u_true,
                u_meas,
                Float64(debris.mag[k])
            )
        )
    end

    return nothing
end


function prepare_display!(img::Matrix{Float64}, percu::Real, percl::Real)

    v = vec(img)
    vmax = quantile(v, percu/100)
    # Lower level defined from zero, not from the image minimum
    vmin = (percl / 100) * vmax

    scale = max(vmax - vmin, eps(Float64))

    img .= clamp.((img .- vmin) ./ scale, 0, 1.0)

    return nothing
end


function save_frame(img::Matrix{Float64}, image_name::AbstractString, outdir::AbstractString, config)
    img_out = copy(img)
    prepare_display!(img_out, config.upper_percentile, config.lower_percentage)

    save(joinpath(outdir, image_name), colorview(Gray, img_out))

    return nothing
end


##################### MAIN ###############################

# Boresight de la cámara en ECI: fila 3 de C_ci (T6: las filas de una DCM son los ejes de llegada)
boresight_eci(C_ci::AbstractMatrix) = SVector{3,Float64}(C_ci[3, 1], C_ci[3, 2], C_ci[3, 3])

"""
    star_frame(stars_all, C_ci, config) -> NamedTuple

Estrellas que ve la cámara con una orientación dada: las que caen en el sensor y sus
píxeles (`sources_in_view`).

Por qué existe: con actitud constante en ECI se llama una sola vez; si la actitud
cambia, una vez por imagen (las estrellas se mueven en la imagen).

# Argumentos
- `stars_all::DataFrame`: catálogo Hipparcos ya filtrado por magnitud.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `config`: `camera_config` (usa `intrinsics`, `mag0` y `flux0`).

# Devuelve
- `(stars, render)`: estrellas en el sensor y sus `(x, y, flux)`.

# Excepciones
- Ninguna.
"""
function star_frame(stars_all::DataFrame, C_ci::AbstractMatrix, config)
    stars, render = sources_in_view(stars_all, C_ci, config)
    return (stars = stars, render = render)
end

"""
    camera_sim!(data, los_measurements, photometry_data, C_ci, sun_dir, obs_name, obj_names, output_path, config)

Simula la cámara de un observador: en cada instante comprueba la exclusión solar
y, si el Sol no deslumbra la cámara, selecciona las estrellas y los RSOs visibles
dentro del FOV, los proyecta con la orientación de la cámara, genera la imagen (si
`image_generation_flag`) y guarda medidas LOS y posiciones en píxeles. En un
instante deslumbrado no hay imagen ni medidas. Las estrellas tapadas por la Tierra
no se pintan.

La orientación de la cámara llega como C_ci (ECI -> cámara) por instante,
calculada como C_cb · C_bi (montaje × actitud del satélite). Si es la misma en
todos los instantes, el campo estelar se calcula una sola vez, como antes; si no,
se recalcula en cada imagen.

# Argumentos
- `data::Dict{Symbol,TrajectoryData}`: trayectorias propagadas (ECI, km).
- `los_measurements`: diccionario de salida con las medidas LOS de cada observador.
- `photometry_data`: magnitud y visibilidad de cada par (observador, objetivo).
- `C_ci::AbstractVector`: DCM ECI -> cámara de este observador, una por instante.
- `sun_dir::AbstractVector`: dirección unitaria del Sol en ECI, una por instante.
- `obs_name::Symbol`: observador.
- `obj_names::Vector{Symbol}`: objetivos.
- `output_path`: carpeta de resultados (se crea `output_path/obs_name/`).
- `config`: `camera_config` con los flags y `sun_exclusion_deg`.

# Devuelve
- `nothing`. Escribe `debris_positions.csv`, `frame_status.csv` (en cada instante: ángulo
  Sol-boresight, si la cámara está deslumbrada, si la Tierra está en la imagen, qué fracción
  ocupa y qué fracción de la imagen es Tierra de día), los PNG (si se piden) y
  `observed_star_catalog.csv` (si `star_csv_flag`; con actitud variable lleva además
  la columna `image_name`).

# Excepciones
- `ErrorException`: si `C_ci` o `sun_dir` no tienen un valor por instante de la trayectoria.
"""
function camera_sim!(data::Dict{Symbol,TrajectoryData}, los_measurements::Dict{Symbol, Dict{Symbol, Vector{LOSMeasurement}}}, photometry_data::Dict{Tuple{Symbol,Symbol},NamedTuple},
    C_ci::AbstractVector, sun_dir::AbstractVector, obs_name::Symbol, obj_names::Vector{Symbol}, output_path, config)

    outdir = joinpath(output_path, String(obs_name))
    mkpath(outdir)

    # --- Trayectorias y medidas de este observador ---
    obs_traj = data[obs_name]
    obj_trajs = [data[name] for name in obj_names]
    length(C_ci) == length(obs_traj.r) ||
        error("camera_sim!: $(length(C_ci)) camera orientations for $(length(obs_traj.r)) time steps of '$obs_name'.")
    length(sun_dir) == length(obs_traj.r) ||
        error("camera_sim!: $(length(sun_dir)) Sun directions for $(length(obs_traj.r)) time steps of '$obs_name'.")
    obs_los = los_measurements[obs_name] = Dict{Symbol, Vector{LOSMeasurement}}()

    # --- Estrellas e imágenes, y tablas de salida ---
    sky = setup_star_field(C_ci, outdir, config)
    debris_log = DataFrame(observer_name = String[], debris_name = String[], image_name = String[],
                           x_pix = Int[], y_pix = Int[], Vmag = String[])
    star_catalog = DataFrame(image_name = String[], x = Int[], y = Int[], Vmag = Float64[])
    frame_status = DataFrame(image_name = String[], t_s = Float64[], sun_angle_deg = Float64[], is_sun_blinded = Bool[],
                             is_earth_in_view = Bool[], earth_fraction = Float64[], earth_daylight_fraction = Float64[])

    # --- Bucle temporal: Sol y Tierra -> estrellas -> RSOs (medidas, píxeles, PSF) -> ruido y PNG ---
    for i in eachindex(obs_traj.r)
        C = C_ci[i]                                # orientación de la cámara en este instante
        image_name = frame_name(obs_traj.t[i])
        status = frame_conditions(C, obs_traj.r[i], sun_dir[i], config)
        push!(frame_status, (image_name = image_name, t_s = Float64(obs_traj.t[i]), status...))
        status.is_sun_blinded && continue          # cámara deslumbrada: sin imagen ni medidas

        star_field = frame_star_field!(star_catalog, sky, C, obs_traj.r[i], image_name, config)
        config.image_generation_flag && start_frame!(sky, star_field, config)

        debris, debris_fov, render_debris = debris_in_view(obj_names, obj_trajs, photometry_data, obs_name, obs_traj.r[i], i, C, config)
        debris_meas = config.ignore_fov_for_measurements_flag ? debris : debris_fov
        nrow(debris_meas) > 0 &&
            err_meas!(obs_los, debris_meas, obs_traj.t[i], config.los_sigma_arcsec, config.navigation_error_flag)
        if nrow(debris_fov) > 0
            log_debris!(debris_log, obs_name, debris_fov, render_debris, image_name)
            config.image_generation_flag && render_psf_gaussian!(sky.img, render_debris, config.intrinsics, config.fwhm_deg)
        end

        config.image_generation_flag && finish_frame!(sky.img, image_name, outdir, config)
    end

    # Con actitud variable, un único CSV con las estrellas de cada imagen
    !sky.is_static && config.star_csv_flag && CSV.write(joinpath(outdir, "observed_star_catalog.csv"), star_catalog)

    # Se escriben siempre (aunque no haya PNG): qué RSOs vio la cámara, cuándo y dónde, y en qué instantes observó
    CSV.write(joinpath(outdir, "debris_positions.csv"), debris_log; delim=';')
    CSV.write(joinpath(outdir, "frame_status.csv"), frame_status)
    return nothing
end

"""
    frame_conditions(C_ci, r_obs, u_sun, config) -> NamedTuple

Condiciones de observación de un instante, las de `frame_status.csv`: ángulo entre el Sol y
el eje óptico, si la cámara está deslumbrada (ángulo menor que `config.sun_exclusion_deg`)
y la Tierra en la imagen (`earth_coverage`).

Por qué existe: separa del bucle de `camera_sim!` la decisión de observar o no y lo que se
registra de cada instante.

# Argumentos
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `r_obs::SVector{3,Float64}`: posición del observador en ECI (km).
- `u_sun::SVector{3,Float64}`: dirección unitaria del Sol en ECI.
- `config`: `camera_config` (usa `sun_exclusion_deg` e `intrinsics`).

# Devuelve
- `(sun_angle_deg, is_sun_blinded, is_earth_in_view, earth_fraction, earth_daylight_fraction)`.

# Excepciones
- Ninguna.
"""
function frame_conditions(C_ci::AbstractMatrix, r_obs::SVector{3,Float64}, u_sun::SVector{3,Float64}, config)
    θ_sun = sun_boresight_angle_deg(C_ci, u_sun)
    sun = (sun_angle_deg = round(θ_sun, digits = 3), is_sun_blinded = θ_sun < config.sun_exclusion_deg)
    return merge(sun, earth_coverage(C_ci, r_obs, u_sun, config.intrinsics))
end

"""
    earth_coverage(C_ci, r_obs, u_sun, cam) -> NamedTuple

Tierra en la imagen de un instante:
- `is_earth_in_view`: si la Tierra toca la imagen. Un disco sin huecos que toca el sensor o
  cruza su borde o cabe dentro, así que basta con mirar los píxeles del borde
  (`is_earth_on_border`) y la malla; solo se escaparía una Tierra más pequeña que una celda
  de la malla y lejos del borde.
- `earth_fraction`: fracción de la imagen que es Tierra, con una malla de n × n rayos
  (n = `EARTH_COVERAGE_GRID_SIZE`): rayos que chocan con la Tierra / n².
- `earth_daylight_fraction`: fracción de la imagen que es Tierra de día: rayos cuyo punto
  de choque P cumple P·û☉ > 0 (sin crepúsculo) / n².

Si el disco de la Tierra no toca el cono de prefiltro (`is_earth_near_view`), no se lanza
ningún rayo y todo vale `false` o 0.

Por qué existe: la Tierra de día satura la cámara y está pendiente decidir qué hacer con
esas imágenes; con estas columnas se pueden filtrar después sin repetir la simulación.

# Argumentos
- `C_ci::AbstractMatrix`: DCM ECI -> cámara; `r_obs`: posición del observador (ECI, km).
- `u_sun::SVector{3,Float64}`: dirección unitaria del Sol en ECI.
- `cam::CameraIntrinsics`.

# Devuelve
- `(is_earth_in_view, earth_fraction, earth_daylight_fraction)`.

# Excepciones
- Ninguna.
"""
function earth_coverage(C_ci::AbstractMatrix, r_obs::SVector{3,Float64}, u_sun::SVector{3,Float64}, cam::CameraIntrinsics)
    is_earth_near_view(C_ci, r_obs, cam) ||
        return (is_earth_in_view = false, earth_fraction = 0.0, earth_daylight_fraction = 0.0)
    n = EARTH_COVERAGE_GRID_SIZE
    hits = daylight = 0
    for v in grid_centers(cam.resolution_y, n), u in grid_centers(cam.resolution_x, n)
        P = earth_ray_hit(r_obs, pixel_ray_eci(cam, C_ci, u, v))
        P === nothing && continue
        hits += 1
        daylight += dot(P, u_sun) > 0
    end
    is_in_view = hits > 0 || is_earth_on_border(C_ci, r_obs, cam)
    return (is_earth_in_view = is_in_view, earth_fraction = hits / n^2, earth_daylight_fraction = daylight / n^2)
end

# ¿Toca el disco de la Tierra el cono de prefiltro? Ángulo eje óptico-nadir γ menos el
# semiángulo de la Tierra ρ⊕, frente al semiángulo del cono (prefiltro rápido)
function is_earth_near_view(C_ci::AbstractMatrix, r_obs::SVector{3,Float64}, cam::CameraIntrinsics)
    nadir = -r_obs / norm(r_obs)
    γ = acosd(clamp(dot(boresight_eci(C_ci), nadir), -1.0, 1.0))
    return γ - earth_angular_radius_deg(r_obs) <= prefilter_cone_deg(cam) / 2
end

# ¿Choca con la Tierra el rayo de algún píxel del borde del sensor?
function is_earth_on_border(C_ci::AbstractMatrix, r_obs::SVector{3,Float64}, cam::CameraIntrinsics)
    nx, ny = cam.resolution_x, cam.resolution_y
    border = Iterators.flatten((((u, v) for u in 1:nx, v in (1, ny)), ((u, v) for u in (1, nx), v in 2:ny-1)))
    return any(px -> earth_ray_hit(r_obs, pixel_ray_eci(cam, C_ci, px...)) !== nothing, border)
end

# Dirección unitaria en ECI del rayo que llega al píxel (u, v): deshace K y la distorsión
# (coordenadas normalizadas) y pasa de ejes cámara a ECI con C_ciᵀ
function pixel_ray_eci(cam::CameraIntrinsics, C_ci::AbstractMatrix, u::Real, v::Real)
    x, y = pixel_to_normalized(cam, u, v)
    return normalize(SVector{3,Float64}(C_ci' * SVector(x, y, 1.0)))
end

# Centros de n celdas iguales a lo largo de N píxeles (el sensor va de ½ a N + ½)
grid_centers(N::Integer, n::Integer) = range(0.5 + N / (2n), N + 0.5 - N / (2n); length = n)

"""
    hide_stars_behind_earth(field, r_obs, C_ci, cam) -> NamedTuple

Campo estelar sin las estrellas que tapa la Tierra vista desde `r_obs` (`is_earth_occulted`).
Si la Tierra está lejos de la imagen (`is_earth_near_view`) o no tapa ninguna, devuelve el
mismo `field`, sin copiarlo.

Por qué existe: las estrellas están en el infinito, así que la Tierra tapa siempre la misma
parte del cielo vista desde un punto; la fotometría ya quita los RSOs tapados, pero las
estrellas no pasan por ella. Con actitud constante, devolver el mismo objeto permite seguir
usando la imagen de fondo precalculada en los instantes sin Tierra.

# Argumentos
- `field`: campo estelar (`star_frame`), o `nothing`.
- `r_obs::SVector{3,Float64}`: posición del observador en ECI (km).
- `C_ci::AbstractMatrix`: DCM ECI -> cámara; `cam::CameraIntrinsics`.

# Devuelve
- `(stars, render)` con las estrellas visibles, o `field` sin cambios.

# Excepciones
- Ninguna.
"""
function hide_stars_behind_earth(field, r_obs::SVector{3,Float64}, C_ci::AbstractMatrix, cam::CameraIntrinsics)
    (field === nothing || !is_earth_near_view(C_ci, r_obs, cam)) && return field
    is_visible = [!is_earth_occulted(r_obs, radec_to_vec(ra, dec)) for (ra, dec) in zip(field.stars.ra, field.stars.dec)]
    all(is_visible) && return field
    r = field.render
    return (stars = field.stars[is_visible, :], render = (x = r.x[is_visible], y = r.y[is_visible], flux = r.flux[is_visible]))
end

# Nombre del PNG de un instante t (s): "120.0.png"
frame_name(t) = replace(string(t), ":" => "", " " => "_") * ".png"

# Filas de observed_star_catalog.csv para un campo estelar: píxel (x, y) y magnitud de cada estrella
star_catalog_rows(star_field) = DataFrame(x = star_field.render.x, y = star_field.render.y, Vmag = star_field.stars.mag)

"""
    setup_star_field(C_ci, outdir, config) -> NamedTuple

Prepara lo que la cámara necesita de las estrellas antes del bucle temporal: carga
el catálogo Hipparcos y, si la actitud es constante en ECI, calcula una sola vez el
campo estelar, su imagen de fondo y `observed_star_catalog.csv`.

Por qué existe: con actitud constante las estrellas no se mueven en la imagen y no
hace falta recalcularlas en cada instante; con actitud variable se recalculan en
`frame_star_field!`.

# Argumentos
- `C_ci::AbstractVector`: DCM ECI -> cámara en cada instante.
- `outdir::AbstractString`: carpeta del observador.
- `config`: `camera_config` con los flags.

# Devuelve
- `(stars_all, is_static, needs_stars, field, img_static, img)`: catálogo, si la
  actitud es constante, si hacen falta estrellas (PNG o CSV), campo estático (o
  `nothing`), imagen de fondo e imagen de trabajo (`nothing` si no se generan PNG).

# Excepciones
- Las de lectura del catálogo FITS (`load_hipparcos`).
"""
function setup_star_field(C_ci::AbstractVector, outdir::AbstractString, config)
    stars_all = load_hipparcos(joinpath(DATA_DIR, "I_239_vmag10_reduced.fits"), config.Vmag_lim)
    is_static = all(C -> C == C_ci[1], C_ci)                   # actitud constante en ECI
    needs_stars = config.image_generation_flag || config.star_csv_flag

    field = is_static && needs_stars ? star_frame(stars_all, C_ci[1], config) : nothing
    field !== nothing && config.star_csv_flag &&
        CSV.write(joinpath(outdir, "observed_star_catalog.csv"), star_catalog_rows(field))

    img_static = config.image_generation_flag ? sky_background(field, config) : nothing
    img = config.image_generation_flag ? similar(img_static) : nothing    # imagen de trabajo, reutilizada
    return (stars_all = stars_all, is_static = is_static, needs_stars = needs_stars,
            field = field, img_static = img_static, img = img)
end

"""
    frame_star_field!(star_catalog, sky, C_ci, r_obs, image_name, config) -> Union{NamedTuple, Nothing}

Campo estelar de un instante: el estático de `setup_star_field` o, con actitud
variable, uno nuevo con la orientación de la cámara, sin las estrellas que tapa la
Tierra (`hide_stars_behind_earth`). Con actitud variable, sus estrellas se añaden a
`star_catalog` (si `star_csv_flag`).

Por qué existe: separa del bucle de `camera_sim!` la diferencia entre actitud
constante y variable.

# Argumentos
- `star_catalog::DataFrame`: estrellas de cada imagen (actitud variable); se modifica.
- `sky`: salida de `setup_star_field`.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `r_obs::SVector{3,Float64}`: posición del observador en ECI (km).
- `image_name::AbstractString`: imagen del instante.
- `config`: `camera_config`.

# Devuelve
- El campo estelar (`star_frame`), o `nothing` si no hacen falta estrellas. Con actitud
  constante y sin estrellas tapadas, es el mismo objeto que `sky.field`.

# Excepciones
- Ninguna.
"""
function frame_star_field!(star_catalog::DataFrame, sky, C_ci::AbstractMatrix, r_obs::SVector{3,Float64},
                           image_name::AbstractString, config)
    sky.needs_stars || return nothing
    cam = config.intrinsics
    sky.is_static && return hide_stars_behind_earth(sky.field, r_obs, C_ci, cam)
    field = hide_stars_behind_earth(star_frame(sky.stars_all, C_ci, config), r_obs, C_ci, cam)
    if config.star_csv_flag
        rows = star_catalog_rows(field)
        insertcols!(rows, 1, :image_name => fill(image_name, nrow(rows)))   # primera columna: imagen
        append!(star_catalog, rows)
    end
    return field
end

"""
    start_frame!(sky, star_field, config) -> Nothing

Deja en la imagen de trabajo `sky.img` las estrellas y el fondo del instante: copia
la imagen de fondo si el campo es el estático sin cambios, o la pinta de nuevo (actitud
variable, o la Tierra tapa estrellas del campo estático).

Por qué existe: es el primer paso de cada imagen en `camera_sim!`, antes de los RSOs.

# Argumentos
- `sky`: salida de `setup_star_field` (con `img` e `img_static`).
- `star_field`: campo estelar del instante (`frame_star_field!`).
- `config`: `camera_config`.

# Devuelve
- `nothing`.

# Excepciones
- Ninguna.
"""
function start_frame!(sky, star_field, config)
    if sky.is_static && star_field === sky.field     # campo estático sin estrellas tapadas
        copyto!(sky.img, sky.img_static)
    else
        fill!(sky.img, 0.0)
        paint_star_field!(sky.img, star_field, config)
    end
    return nothing
end

"""
    paint_star_field!(img, star_field, config) -> Nothing

Pinta en `img` las estrellas de `star_field` con la PSF gaussiana y suma el fondo
de cielo uniforme.

Por qué existe: es lo mismo con actitud constante (una vez) y variable (en cada imagen).

# Argumentos
- `img::Matrix{Float64}`: imagen (ny × nx); se modifica.
- `star_field`: salida de `star_frame`.
- `config`: `camera_config` (usa `intrinsics`, `fwhm_deg`, `background_level`).

# Devuelve
- `nothing`.

# Excepciones
- Ninguna.
"""
function paint_star_field!(img::Matrix{Float64}, star_field, config)
    render_psf_gaussian!(img, star_field.render, config.intrinsics, config.fwhm_deg)
    add_background!(img, config.background_level)
    return nothing
end

"""
    sky_background(star_field, config) -> Matrix{Float64}

Imagen de partida de cada instante con actitud constante en ECI: estrellas y fondo
de cielo, que no cambian entre imágenes. Con actitud variable (`star_field = nothing`)
devuelve una imagen vacía, que solo sirve de molde para la imagen de trabajo.

Por qué existe: con actitud constante el campo estelar se pinta una sola vez y en
cada instante se copia, en lugar de recalcularlo.

# Argumentos
- `star_field`: salida de `star_frame`, o `nothing` si la actitud cambia.
- `config`: `camera_config` (resolución, `intrinsics`, `fwhm_deg`, `background_level`).

# Devuelve
- `Matrix{Float64}` de tamaño `resolution_y × resolution_x`.

# Excepciones
- Ninguna.
"""
function sky_background(star_field, config)
    img = zeros(Float64, config.resolution_y, config.resolution_x)
    star_field === nothing || paint_star_field!(img, star_field, config)
    return img
end

"""
    debris_in_view(obj_names, obj_trajs, photometry_data, obs_name, r_obs, i, C_ci, config) -> (DataFrame, DataFrame, NamedTuple)

RSOs visibles desde el observador en el instante `i` (iluminados y en LOS, según la
fotometría) y, de ellos, los que caen en el sensor de la cámara
(`sources_in_view`), con sus píxeles.

Por qué existe: separa la selección de RSOs del bucle de `camera_sim!`.

# Argumentos
- `obj_names`, `obj_trajs`: objetivos y sus trayectorias.
- `photometry_data`: magnitud y visibilidad de cada par (observador, objetivo).
- `obs_name::Symbol`, `r_obs::SVector{3,Float64}`: observador y su posición (ECI, km).
- `i::Int`: índice del instante.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `config`: `camera_config` (usa `intrinsics`, `mag0` y `flux0`).

# Devuelve
- `(debris, debris_fov, render_debris)`: tablas con `obj_name`, `ra`, `dec` (°) y `mag`,
  y los `(x, y, flux)` de los que caen en el sensor.

# Excepciones
- Ninguna.
"""
function debris_in_view(obj_names, obj_trajs, photometry_data, obs_name::Symbol, r_obs::SVector{3,Float64},
                        i::Int, C_ci::AbstractMatrix, config)
    debris = build_debris_catalog(obj_names, obj_trajs, photometry_data, obs_name, r_obs, i)
    nrow(debris) == 0 && return debris, debris, (x = Int[], y = Int[], flux = Float64[])

    debris_fov, render_debris = sources_in_view(debris, C_ci, config)
    return debris, debris_fov, render_debris
end

"""
    log_debris!(debris_log, obs_name, debris_fov, render_debris, image_name) -> Nothing

Añade a `debris_log` una fila por RSO en el FOV: observador, objetivo, imagen,
píxel (x, y) y magnitud (texto con 3 decimales y coma decimal).

Por qué existe: es el contenido de `debris_positions.csv`; separa su formato del bucle.

# Argumentos
- `debris_log::DataFrame`: tabla de salida; se modifica.
- `obs_name::Symbol`: observador.
- `debris_fov::DataFrame`: RSOs en el FOV (`debris_in_view`).
- `render_debris`: sus píxeles (`prepare_render_dict`).
- `image_name::AbstractString`: imagen (instante) a la que corresponden.

# Devuelve
- `nothing`.

# Excepciones
- Ninguna.
"""
function log_debris!(debris_log::DataFrame, obs_name::Symbol, debris_fov::DataFrame, render_debris,
                     image_name::AbstractString)
    for j in 1:nrow(debris_fov)
        push!(debris_log, (
            observer_name = String(obs_name),
            debris_name   = String(debris_fov.obj_name[j]),
            image_name    = image_name,
            x_pix         = render_debris.x[j],
            y_pix         = render_debris.y[j],
            Vmag          = replace(string(round(debris_fov.mag[j], digits = 3)), "." => ",")
        ))
    end
    return nothing
end

"""
    finish_frame!(img, image_name, outdir, config) -> Nothing

Termina la imagen de un instante: le añade ruido de Poisson (fotones) y de lectura
del sensor y la guarda como PNG en `outdir`.

Por qué existe: es el último paso de cada imagen en `camera_sim!`, después de pintar
estrellas y RSOs.

# Argumentos
- `img::Matrix{Float64}`: imagen sin ruido; se modifica.
- `image_name::AbstractString`: nombre del PNG (instante).
- `outdir::AbstractString`: carpeta del observador.
- `config`: `camera_config` (usa `sigma_read_noise` y los percentiles de `save_frame`).

# Devuelve
- `nothing`.

# Excepciones
- Las de escritura del PNG (carpeta sin permisos, disco lleno).
"""
function finish_frame!(img::Matrix{Float64}, image_name::AbstractString, outdir::AbstractString, config)
    add_poisson_noise!(img)
    add_read_noise!(img, config.sigma_read_noise)
    save_frame(img, image_name, outdir, config)
    return nothing
end

end # module CameraSim