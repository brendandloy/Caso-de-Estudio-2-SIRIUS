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
# z_c = boresight. Proyección pinhole (= gnomónica TAN) con f = ½·mín(nx, ny)/tan(FOV/2).
# El FOV es el cono de semiángulo FOV/2 alrededor del boresight (como en el CE3).
#
# Funciones principales: load_hipparcos, prefilter_box, filter_cone, radec_to_pixels!,
# prepare_render_dict, render_psf_gaussian!, ruido (add_*), build_debris_catalog, err_meas!,
# star_frame, camera_sim! (punto de entrada) y sus pasos: setup_star_field, frame_star_field!,
# start_frame!, paint_star_field!, sky_background, debris_in_view, log_debris!, finish_frame!,
# frame_name y star_catalog_rows.
# =============================================================================

module CameraSim

export camera_sim!

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


function prefilter_box(stars, boresight_vec, fov_deg::Real)

    # Convert camera pointing to RA/Dec center
    ra0, dec0 = vec_to_radec(boresight_vec)

    # Build bounding box in RA/Dec
    box = radec_box_square(ra0, dec0, fov_deg)

    # Apply declination bounds
    mask = (stars.dec .>= box.dec_min) .& (stars.dec .<= box.dec_max)

    # Apply RA bounds (may be full range if pole/wrap was detected)
    # mask .&= (stars.ra .>= box.ra_min) .& (stars.ra .<= box.ra_max)

    return stars[mask, :]
end

function filter_cone(stars, boresight_vec, fov_deg::Real)

    bx, by, bz = boresight_vec
    cos_fov = cos(deg2rad(fov_deg / 2))

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
    radec_to_pixels!(x, y, ra, dec, C_ci, fov_deg, nx, ny) -> Nothing

Proyecta direcciones (RA, Dec) sobre la imagen con el modelo pinhole: se pasa
cada dirección a ejes cámara con la orientación de la cámara y se divide por z_c.

Modelo (ejes cámara dextrógiros: x_c derecha, y_c abajo, z_c = boresight):
  ρ_c = C_ci · û_i,   x = c_x + f · x_c/z_c,   y = c_y + f · y_c/z_c,
con f = ½·mín(nx, ny) / tan(FOV/2) (focal en píxeles, como en el CE3 de I. Ferrando)
y (c_x, c_y) el centro de la imagen. Es la proyección gnomónica (TAN) sobre el plano
tangente en el boresight, pero con el giro de la imagen que da la actitud.

# Argumentos
- `x::Vector{Int}`, `y::Vector{Int}`: salida, columna y fila (1 = arriba a la izquierda).
- `ra::Vector{Float64}`, `dec::Vector{Float64}`: direcciones en ECI/ICRS (°).
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `fov_deg::Real`: campo de visión (°) sobre el lado menor de la imagen.
- `nx::Int`, `ny::Int`: resolución (píxeles).

# Devuelve
- `nothing` (rellena `x` e `y`). Las fuentes detrás de la cámara reciben `OFF_IMAGE`.

# Excepciones
- `AssertionError`: si los vectores no tienen la misma longitud.

# Ejemplo
```julia
x = zeros(Int, 1); y = zeros(Int, 1)
radec_to_pixels!(x, y, [0.0], [0.0], C_ci, 27, 1024, 1024)    # boresight -> centro (512 o 513)
```
"""
function radec_to_pixels!(x::Vector{Int}, y::Vector{Int},
                          ra::Vector{Float64}, dec::Vector{Float64},
                          C_ci::AbstractMatrix,
                          fov_deg::Real, nx::Int, ny::Int)

    @assert length(x) == length(y) == length(ra) == length(dec)

    C = SMatrix{3,3,Float64,9}(C_ci)

    # Cone half-angle in radians
    θ = deg2rad(Float64(fov_deg) / 2.0)

    # Focal length in pixels:
    # Ensures the cone circle fits inside the pixel square
    half_min = 0.5 * min(Float64(nx), Float64(ny))
    f = half_min / tan(θ)

    # Image center in 1-based pixel coords
    cx = (Float64(nx) + 1.0) / 2.0
    cy = (Float64(ny) + 1.0) / 2.0

    @inbounds for i in eachindex(ra)
        # Dirección en ejes cámara
        ρ = C * radec_to_vec(ra[i], dec[i])

        # Detrás del plano imagen (no ocurre si se ha filtrado con el cono)
        if ρ[3] <= 0.0
            x[i] = OFF_IMAGE
            y[i] = OFF_IMAGE
            continue
        end

        # Pinhole: x hacia la derecha y y hacia abajo, como los índices de la imagen
        x[i] = round(Int, cx + f * ρ[1] / ρ[3])
        y[i] = round(Int, cy + f * ρ[2] / ρ[3])
    end

    return nothing
end


"""
    prepare_render_dict(stars, C_ci, fov_deg, nx, ny; mag0=8.0, flux0=20_000.0) -> NamedTuple

Calcula los píxeles y el flujo relativo de una lista de fuentes (estrellas o RSOs)
para una orientación de la cámara.

# Argumentos
- `stars::DataFrame`: fuentes con columnas `ra`, `dec` (°) y `mag`.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `fov_deg::Real`, `nx::Int`, `ny::Int`: campo de visión (°) y resolución.
- `mag0`, `flux0`: magnitud y flujo de referencia de la escala de flujo.

# Devuelve
- `(x, y, flux)`: columnas, filas y flujos relativos.

# Excepciones
- Ninguna.
"""
function prepare_render_dict(stars::DataFrame, C_ci::AbstractMatrix, fov_deg::Real,
                           nx::Int, ny::Int; mag0::Real = 8.0, flux0::Real = 20_000.0)

    n = nrow(stars)

    # Force dense, type-stable arrays (fast in hot loops)
    ra  = Vector{Float64}(stars.ra)
    dec = Vector{Float64}(stars.dec)
    mag = Vector{Float64}(stars.mag)

    x = Vector{Int}(undef, n)
    y = Vector{Int}(undef, n)

    # Project to pixels with the camera orientation
    radec_to_pixels!(x, y, ra, dec, C_ci, fov_deg, nx, ny)

    # Magnitude -> flux (relative scaling)
    flux = similar(x, Float64)
    @inbounds for i in eachindex(mag)
        flux[i] = mag_to_flux(mag[i], mag0, flux0)
    end

    return (x = x, y = y, flux = flux)
end


function gaussian_kernel(fwhm_pix::Real; k_sigma::Real = 4.0)

    σ = Float64(fwhm_pix) / 2.3548200450309493
    rmax = ceil(Int, k_sigma * σ)
    inv2σ2 = 1.0 / (2σ^2)

    sz = 2*rmax + 1
    K = Matrix{Float64}(undef, sz, sz)

    @inbounds for j in 1:sz, i in 1:sz
        dx = Float64(i - (rmax + 1))
        dy = Float64(j - (rmax + 1))
        K[j, i] = exp(-(dx*dx + dy*dy) * inv2σ2)
    end

    # Normalize once so total kernel sum = 1
    K ./= sum(K)

    return K, rmax
end


function render_psf_gaussian!(img::Matrix{Float64}, render::NamedTuple, fov_deg::Real,
                              fwhm_deg::Real;
                              k_sigma::Real = 4.0)

    x = render.x
    y = render.y
    flux = render.flux

    ny, nx = size(img)

    fwhm_pix = fwhm_deg / (fov_deg / nx)
    K, rmax = gaussian_kernel(fwhm_pix; k_sigma=k_sigma)
    sz = 2*rmax + 1

    @inbounds for k in eachindex(x)

        ix0 = x[k]
        iy0 = y[k]

        # Check if kernel is in image
        if ix0 + rmax < 1 || ix0 - rmax > nx || iy0 + rmax < 1 || iy0 - rmax > ny
            continue
        end

        F = flux[k]

        for kj in 1:sz
            j = iy0 + (kj - (rmax + 1))
            if j < 1 || j > ny
                continue
            end
            for ki in 1:sz
                i = ix0 + (ki - (rmax + 1))
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

Estrellas que ve la cámara con una orientación dada: filtro por el cono del
FOV alrededor del boresight y proyección a píxeles.

Por qué existe: con actitud constante en ECI se llama una sola vez; si la actitud
cambia, una vez por imagen (las estrellas se mueven en la imagen).

# Argumentos
- `stars_all::DataFrame`: catálogo Hipparcos ya filtrado por magnitud.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `config`: `camera_config` (usa `fov_deg`, resoluciones, `mag0` y `flux0`).

# Devuelve
- `(stars, render)`: estrellas dentro del FOV y su `prepare_render_dict`.

# Excepciones
- Ninguna.
"""
function star_frame(stars_all::DataFrame, C_ci::AbstractMatrix, config)
    fov = config.fov_deg
    b = boresight_eci(C_ci)
    stars = prefilter_box(stars_all, b, fov)
    stars = filter_cone(stars, b, fov)
    render = prepare_render_dict(stars, C_ci, fov, config.resolution_x, config.resolution_y;
        mag0 = config.mag0, flux0 = config.flux0)
    return (stars = stars, render = render)
end

"""
    camera_sim!(data, los_measurements, photometry_data, C_ci, obs_name, obj_names, output_path, config)

Simula la cámara de un observador: en cada instante selecciona las estrellas y
los RSOs visibles dentro del FOV, los proyecta con la orientación de la cámara,
genera la imagen (si `image_generation_flag`) y guarda medidas LOS y posiciones
en píxeles.

La orientación de la cámara llega como C_ci (ECI -> cámara) por instante,
calculada como C_cb · C_bi (montaje × actitud del satélite). Si es la misma en
todos los instantes, el campo estelar se calcula una sola vez, como antes; si no,
se recalcula en cada imagen.

# Argumentos
- `data::Dict{Symbol,TrajectoryData}`: trayectorias propagadas (ECI, km).
- `los_measurements`: diccionario de salida con las medidas LOS de cada observador.
- `photometry_data`: magnitud y visibilidad de cada par (observador, objetivo).
- `C_ci::AbstractVector`: DCM ECI -> cámara de este observador, una por instante.
- `obs_name::Symbol`: observador.
- `obj_names::Vector{Symbol}`: objetivos.
- `output_path`: carpeta de resultados (se crea `output_path/obs_name/`).
- `config`: `camera_config` con los flags.

# Devuelve
- `nothing`. Escribe `debris_positions.csv`, los PNG (si se piden) y
  `observed_star_catalog.csv` (si `star_csv_flag`; con actitud variable lleva además
  la columna `image_name`).

# Excepciones
- `ErrorException`: si `C_ci` no tiene un valor por instante de la trayectoria.
"""
function camera_sim!(data::Dict{Symbol,TrajectoryData}, los_measurements::Dict{Symbol, Dict{Symbol, Vector{LOSMeasurement}}}, photometry_data::Dict{Tuple{Symbol,Symbol},NamedTuple},
    C_ci::AbstractVector, obs_name::Symbol, obj_names::Vector{Symbol}, output_path, config)

    outdir = joinpath(output_path, String(obs_name))
    mkpath(outdir)

    # --- Trayectorias y medidas de este observador ---
    obs_traj = data[obs_name]
    obj_trajs = [data[name] for name in obj_names]
    length(C_ci) == length(obs_traj.r) ||
        error("camera_sim!: $(length(C_ci)) camera orientations for $(length(obs_traj.r)) time steps of '$obs_name'.")
    obs_los = los_measurements[obs_name] = Dict{Symbol, Vector{LOSMeasurement}}()

    # --- Estrellas e imágenes, y tablas de salida ---
    sky = setup_star_field(C_ci, outdir, config)
    debris_log = DataFrame(observer_name = String[], debris_name = String[], image_name = String[],
                           x_pix = Int[], y_pix = Int[], Vmag = String[])
    star_catalog = DataFrame(image_name = String[], x = Int[], y = Int[], Vmag = Float64[])

    # --- Bucle temporal: estrellas -> RSOs (medidas, píxeles, PSF) -> ruido y PNG ---
    for i in eachindex(obs_traj.r)
        C = C_ci[i]                                # orientación de la cámara en este instante
        image_name = frame_name(obs_traj.t[i])

        star_field = frame_star_field!(star_catalog, sky, C, image_name, config)
        config.image_generation_flag && start_frame!(sky, star_field, config)

        debris, debris_fov = debris_in_view(obj_names, obj_trajs, photometry_data, obs_name, obs_traj.r[i], i, C, config)
        debris_meas = config.ignore_fov_for_measurements_flag ? debris : debris_fov
        nrow(debris_meas) > 0 &&
            err_meas!(obs_los, debris_meas, obs_traj.t[i], config.los_sigma_arcsec, config.navigation_error_flag)
        if nrow(debris_fov) > 0
            render_debris = prepare_render_dict(debris_fov, C, config.fov_deg, config.resolution_x, config.resolution_y;
                                                mag0 = config.mag0, flux0 = config.flux0)
            log_debris!(debris_log, obs_name, debris_fov, render_debris, image_name)
            config.image_generation_flag && render_psf_gaussian!(sky.img, render_debris, config.fov_deg, config.fwhm_deg)
        end

        config.image_generation_flag && finish_frame!(sky.img, image_name, outdir, config)
    end

    # Con actitud variable, un único CSV con las estrellas de cada imagen
    !sky.is_static && config.star_csv_flag && CSV.write(joinpath(outdir, "observed_star_catalog.csv"), star_catalog)

    # Se escribe siempre (aunque no haya PNG): indica qué RSOs vio la cámara, cuándo y dónde
    CSV.write(joinpath(outdir, "debris_positions.csv"), debris_log; delim=';')
    return nothing
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
    frame_star_field!(star_catalog, sky, C_ci, image_name, config) -> Union{NamedTuple, Nothing}

Campo estelar de un instante: el estático de `setup_star_field` o, con actitud
variable, uno nuevo con la orientación de la cámara, cuyas estrellas se añaden a
`star_catalog` (si `star_csv_flag`).

Por qué existe: separa del bucle de `camera_sim!` la diferencia entre actitud
constante y variable.

# Argumentos
- `star_catalog::DataFrame`: estrellas de cada imagen (actitud variable); se modifica.
- `sky`: salida de `setup_star_field`.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `image_name::AbstractString`: imagen del instante.
- `config`: `camera_config`.

# Devuelve
- El campo estelar (`star_frame`), o `nothing` si no hacen falta estrellas.

# Excepciones
- Ninguna.
"""
function frame_star_field!(star_catalog::DataFrame, sky, C_ci::AbstractMatrix, image_name::AbstractString, config)
    (sky.is_static || !sky.needs_stars) && return sky.field
    field = star_frame(sky.stars_all, C_ci, config)
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
la imagen de fondo con actitud constante, o la pinta de nuevo con actitud variable.

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
    if sky.is_static
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
- `config`: `camera_config` (usa `fov_deg`, `fwhm_deg`, `background_level`).

# Devuelve
- `nothing`.

# Excepciones
- Ninguna.
"""
function paint_star_field!(img::Matrix{Float64}, star_field, config)
    render_psf_gaussian!(img, star_field.render, config.fov_deg, config.fwhm_deg)
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
- `config`: `camera_config` (resolución, `fov_deg`, `fwhm_deg`, `background_level`).

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
    debris_in_view(obj_names, obj_trajs, photometry_data, obs_name, r_obs, i, C_ci, config) -> (DataFrame, DataFrame)

RSOs visibles desde el observador en el instante `i` (iluminados, en LOS y sin
deslumbramiento, según la fotometría) y, de ellos, los que caen en el FOV de la cámara.

Por qué existe: separa la selección de RSOs del bucle de `camera_sim!`.

# Argumentos
- `obj_names`, `obj_trajs`: objetivos y sus trayectorias.
- `photometry_data`: magnitud y visibilidad de cada par (observador, objetivo).
- `obs_name::Symbol`, `r_obs::SVector{3,Float64}`: observador y su posición (ECI, km).
- `i::Int`: índice del instante.
- `C_ci::AbstractMatrix`: DCM ECI -> cámara en este instante.
- `config`: `camera_config` (usa `fov_deg`).

# Devuelve
- `(debris, debris_fov)`: tablas con `obj_name`, `ra`, `dec` (°) y `mag`.

# Excepciones
- Ninguna.
"""
function debris_in_view(obj_names, obj_trajs, photometry_data, obs_name::Symbol, r_obs::SVector{3,Float64},
                        i::Int, C_ci::AbstractMatrix, config)
    debris = build_debris_catalog(obj_names, obj_trajs, photometry_data, obs_name, r_obs, i)
    nrow(debris) == 0 && return debris, debris

    b = boresight_eci(C_ci)
    debris_fov = filter_cone(prefilter_box(debris, b, config.fov_deg), b, config.fov_deg)
    return debris, debris_fov
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