# =============================================================================
# Relative_utils.jl
#
# Helpers de movimiento relativo para las disposiciones basadas en las
# ecuaciones de Clohessy-Wiltshire (CW). Se incluyen desde Scenario.jl, por lo
# que comparten su espacio de nombres (ven elements_to_rv, cartesian_to_kepler...).
#
# Sistema local (LVLH / Hill) del chief:
#   x : radial (hacia afuera)    y : a lo largo de la traza    z : normal al plano orbital
#
# CW solo se usa para ELEGIR las condiciones iniciales de cada miembro; después el
# simulador propaga cada satélite con su dinámica completa. Hipótesis de CW: chief
# casi circular, distancias relativas << SMA, gravedad kepleriana.
#
# Constantes:
#   CW_MAX_ECC, CW_MAX_RATIO
#
# Funciones:
#   lvlh_to_elements : estado relativo LVLH (posición y velocidad) -> elementos del miembro.
#   cw_members       : lista de estados relativos -> lista de miembros (con avisos de validez).
# =============================================================================

# Excentricidad del chief por encima de la cual CW deja de ser una buena aproximación
const CW_MAX_ECC = 0.01

# Cociente distancia relativa / SMA por encima del cual la linealización pierde precisión
const CW_MAX_RATIO = 1e-3

"""
    lvlh_to_elements(chief, rho, rho_dot, ctx) -> Dict{String,Any}

Convierte un estado relativo en el sistema LVLH del chief en los elementos
keplerianos del satélite correspondiente.

Pasos: base LVLH `[x̂ ŷ ẑ]` con `x̂ = r/|r|`, `ẑ = h/|h|`, `ŷ = ẑ × x̂`; posición
inercial `r_c + Q·ρ`; velocidad inercial `v_c + Q·ρ̇ + ω × (Q·ρ)`, donde
`ω = h/|r|²` es la velocidad angular del sistema LVLH. El término `ω × ρ` es
necesario porque ρ̇ se mide en ejes que giran con el chief.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del chief (de `build_orbit`).
- `rho::AbstractVector`: posición relativa `[x, y, z]` (km) en LVLH.
- `rho_dot::AbstractVector`: velocidad relativa `[ẋ, ẏ, ż]` (km/s) medida en LVLH.
- `ctx`: contexto del Builder (usa `ctx.mu`).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA` (km y grados).

# Excepciones
- `ErrorException`: si el estado resultante no es una órbita ligada.

# Ejemplo
```julia
chief = Dict("SMA" => 7000.0, "ECC" => 0.0, "INC" => 90.0, "RAAN" => 0.0, "AOP" => 0.0, "TA" => 0.0)
lvlh_to_elements(chief, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], ctx)      # 1 km por delante: AOP + TA ≈ 0.0082 grados, SMA ≈ 7000.0003 km
```
"""
function lvlh_to_elements(chief::AbstractDict, rho::AbstractVector, rho_dot::AbstractVector, ctx)
    r_c, v_c = elements_to_rv(chief, ctx.mu)
    h = cross(r_c, v_c)

    # Base LVLH expresada en el sistema inercial (columnas x̂, ŷ, ẑ)
    ex = r_c / norm(r_c)
    ez = h / norm(h)
    ey = cross(ez, ex)
    Q = hcat(ex, ey, ez)

    ω = h / dot(r_c, r_c)                       # velocidad angular del sistema LVLH (inercial)
    ρ_I = Q * rho
    ρdot_I = Q * rho_dot + cross(ω, ρ_I)

    return cartesian_to_kepler(r_c + ρ_I, v_c + ρdot_I, ctx.mu)
end

"""
    cw_members(chief, states, layout, ctx; who) -> Vector{Dict{String,Any}}

Convierte una lista de estados relativos LVLH en miembros de la formación y
avisa si las hipótesis de CW no se cumplen. Si `layout` contiene
`"include_chief": true`, el chief se añade como primer miembro.

Por qué existe: las tres disposiciones CW (cw_state, cw_ellipse, pco) solo
difieren en cómo generan los estados; el resto es común y vive aquí.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del chief.
- `states::AbstractVector`: lista de vectores `[x, y, z, ẋ, ẏ, ż]` (km y km/s).
- `layout::AbstractDict`: bloque `layout` (se lee `include_chief`).
- `ctx`: contexto del Builder.
- `who::AbstractString`: nombre de la disposición, para los mensajes.

# Devuelve
- `Vector` de diccionarios de elementos (el chief primero si `include_chief`).

# Excepciones
- `ErrorException`: si `include_chief` no es `true` o `false`.
- `ErrorException`: si algún estado no produce una órbita ligada.

# Ejemplo
```julia
states = [[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]]
cw_members(chief, states, Dict("include_chief" => true), ctx; who = "cw_state")   # [chief, deputy]
```
"""
function cw_members(chief::AbstractDict, states::AbstractVector, layout::AbstractDict, ctx;
                    who::AbstractString)
    include_chief = get(layout, "include_chief", false)
    include_chief isa Bool || error("'$who': 'include_chief' must be true or false.")

    # Avisos de validez de la linealización
    chief["ECC"] > CW_MAX_ECC &&
        @warn "'$who': chief ECC = $(chief["ECC"]) > $CW_MAX_ECC; CW assumes a near-circular chief."
    ratio = maximum(norm(s[1:3]) for s in states) / chief["SMA"]
    ratio > CW_MAX_RATIO &&
        @warn "'$who': relative distance / SMA = $(round(ratio, sigdigits=2)) > $CW_MAX_RATIO; CW loses accuracy."

    members = Dict{String,Any}[]
    include_chief && push!(members, Dict{String,Any}(chief))

    for s in states
        push!(members, lvlh_to_elements(chief, s[1:3], s[4:6], ctx))
    end
    return members
end
