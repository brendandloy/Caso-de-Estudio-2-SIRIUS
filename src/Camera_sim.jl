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
# star_frame, camera_sim! (punto de entrada).
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
    pole_cross = (dec_max_raw >= 90.0) || (dec_min_raw <= -90.0)

    c = cos(deg2rad(dec0))
    dra = pole_cross ? 180.0 : (half / max(c, 1e-12))  # 180 deg makes RA full-range

    ra_min = wrap360(ra0 - dra)
    ra_max = wrap360(ra0 + dra)

    # If FOV contains RA = 0
    ra_wrap = ra_min > ra_max 

    # If RA wraps, use full RA range to keep the output simple
    ra_min = (pole_cross || ra_wrap) ? 0.0 : ra_min
    ra_max = (pole_cross || ra_wrap) ? 360.0 : ra_max

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

    keep = BitVector(undef, nrow(stars))

    for i in 1:nrow(stars)
        # Convert RA/Dec to a unit vector (ICRS)
        α = deg2rad(stars.ra[i])
        δ = deg2rad(stars.dec[i])
        cδ = cos(δ)

        sx = cδ * cos(α)
        sy = cδ * sin(α)
        sz = sin(δ)

        # Keep star if angular separation <= radius
        keep[i] = (sx*bx + sy*by + sz*bz) >= cos_fov
    end

    return stars[keep, :]
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
        phot.visible[i] || continue

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


function err_meas!(obs_measurements::Dict{Symbol,Vector{LOSMeasurement}}, debris::DataFrame, t::Float64, los_sigma_arcsec::Real, meas_flag::Bool)

    σ_los = deg2rad(Float64(los_sigma_arcsec) / 3600.0)

    @inbounds for k in 1:nrow(debris)

        obj_name = debris.obj_name[k]

        u_true = radec_to_vec(debris.ra[k], debris.dec[k])
        u_meas = meas_flag ? perturb_los(u_true, σ_los) : SVector{3,Float64}(NaN, NaN, NaN)

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

    path = joinpath(DATA_DIR, "I_239_vmag10_reduced.fits")
    outdir = joinpath(output_path, String(obs_name))
    mkpath(outdir)

    fov = config.fov_deg
    nx, ny = config.resolution_x, config.resolution_y

    debris_log = DataFrame(
        observer_name = String[],
        debris_name   = String[],
        image_name    = String[],
        x_pix         = Int[],
        y_pix         = Int[],
        Vmag          = String[]
    )

    # --- Extract trajectories once ---
    obs_traj = data[obs_name]
    obj_trajs = [data[name] for name in obj_names]
    length(C_ci) == length(obs_traj.r) ||
        error("camera_sim!: $(length(C_ci)) camera orientations for $(length(obs_traj.r)) time steps of '$obs_name'.")

    los_measurements[obs_name] = Dict{Symbol, Vector{LOSMeasurement}}()

    # --- Catálogo estelar (una vez) y campo estelar ---
    stars_all = load_hipparcos(path, config.Vmag_lim)
    static = all(C -> C == C_ci[1], C_ci)          # actitud constante en ECI: mismas estrellas siempre
    need_stars = config.image_generation_flag || config.star_csv_flag

    # Estrellas en la imagen: columnas x, y, Vmag (+ image_name si la actitud cambia)
    star_catalog = DataFrame(image_name = String[], x = Int[], y = Int[], Vmag = Float64[])

    if static && need_stars
        sf = star_frame(stars_all, C_ci[1], config)
        if config.star_csv_flag
            CSV.write(joinpath(outdir, "observed_star_catalog.csv"),
                      DataFrame(x = sf.render.x, y = sf.render.y, Vmag = sf.stars.mag))
        end
    end

    if config.image_generation_flag
        # Static image: stars + deterministic background (solo si la actitud es constante)
        img_static = zeros(Float64, ny, nx)
        if static
            render_psf_gaussian!(img_static, sf.render, fov, config.fwhm_deg)
            add_background!(img_static, config.background_level)
        end

        # Working image reused at each frame
        img = similar(img_static)
    end

    # --- Time loop ---
    @inbounds for i in eachindex(obs_traj.r)

        C = C_ci[i]                                # orientación de la cámara en este instante
        b = boresight_eci(C)
        image_name = replace(string(obs_traj.t[i]), ":" => "", " " => "_") * ".png"

        if config.image_generation_flag
            if static
                copyto!(img, img_static)
            else
                # Actitud variable: las estrellas se mueven en la imagen y se pintan de nuevo
                fill!(img, 0.0)
                sf = star_frame(stars_all, C, config)
                render_psf_gaussian!(img, sf.render, fov, config.fwhm_deg)
                add_background!(img, config.background_level)
            end
        end

        if !static && config.star_csv_flag
            config.image_generation_flag || (sf = star_frame(stars_all, C, config))
            append!(star_catalog, DataFrame(image_name = fill(image_name, nrow(sf.stars)),
                                            x = sf.render.x, y = sf.render.y, Vmag = sf.stars.mag))
        end

        r_obs = obs_traj.r[i]

        debris = build_debris_catalog(obj_names, obj_trajs, photometry_data, obs_name, r_obs, i)

        if nrow(debris) > 0
            debris_fov = prefilter_box(debris, b, fov)
            debris_fov = filter_cone(debris_fov, b, fov)

            debris_meas = config.ignore_fov_for_measurements_flag ? debris : debris_fov

            if nrow(debris_meas) > 0
                err_meas!(los_measurements[obs_name], debris_meas, obs_traj.t[i], config.los_sigma_arcsec, config.navigation_error_flag)
            end

            # Posiciones en píxeles de los RSOs en el FOV: se registran siempre (debris_positions.csv);
            # image_name identifica el instante aunque no se guarde el PNG
            if nrow(debris_fov) > 0

                render_debris = prepare_render_dict(debris_fov, C, fov, nx, ny;
                    mag0 = config.mag0, flux0 = config.flux0)
                    
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

                # Solo se pintan en la imagen si se generan los PNG
                if config.image_generation_flag
                    render_psf_gaussian!(img, render_debris, fov, config.fwhm_deg)
                end
            end
        end

        if config.image_generation_flag
            add_poisson_noise!(img)
            add_read_noise!(img, config.sigma_read_noise)
            save_frame(img, image_name, outdir, config)
        end

    end

    # Con actitud variable, un único CSV con las estrellas de cada imagen
    if !static && config.star_csv_flag
        CSV.write(joinpath(outdir, "observed_star_catalog.csv"), star_catalog)
    end

    # Se escribe siempre (aunque no haya PNG): indica qué RSOs vio la cámara, cuándo y dónde
    CSV.write(joinpath(outdir, "debris_positions.csv"), debris_log; delim=';')
    return nothing
    
end

end # module CameraSim