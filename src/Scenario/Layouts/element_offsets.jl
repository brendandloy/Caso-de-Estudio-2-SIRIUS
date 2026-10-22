# =============================================================================
# element_offsets.jl
#
# Plugin de disposición: cada miembro es el jefe más un desplazamiento en sus
# elementos orbitales. Registra la disposición "element_offsets".
#
#   "layout": {"type": "element_offsets", "offsets": [
#       {},                                   # el jefe tal cual
#       {"dM": 0.05},                         # desfase en anomalía media (deg)
#       {"dINC": 0.01, "dRAAN": 0.01}         # órbita ligeramente inclinada
#   ]}
#
# Claves de cada offset (todas opcionales, 0 si faltan):
#   dSMA (km), dECC, dINC, dRAAN, dAOP, dM (deg)
# Se lista un offset por miembro; "{}" deja el miembro en la órbita del jefe.
#
# Física (perturbaciones pequeñas respecto al jefe):
#   dM         -> separación a lo largo de la traza
#   dINC/dRAAN -> separación fuera del plano (cruces de plano cada medio periodo)
#   dECC       -> oscilación radial/along-track con el periodo orbital
#   dSMA       -> deriva secular: periodos distintos hacen que se separen con el tiempo
# El desplazamiento se aplica en anomalía media; la TA final se recalcula con la
# ECC ya modificada.
#
# Funciones: element_offsets_layout
# =============================================================================

# Claves admitidas en cada offset
const OFFSET_KEYS = ("dSMA", "dECC", "dINC", "dRAAN", "dAOP", "dM")

"""
    element_offsets_layout(chief, layout, ctx) -> Vector{Dict{String,Any}}

Genera un miembro por cada entrada de `offsets`, sumando los desplazamientos a
los elementos del jefe.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del jefe (de `build_orbit`).
- `layout::AbstractDict`: bloque `layout` con la lista `offsets`.
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` de diccionarios de elementos, uno por offset y en el mismo orden.

# Excepciones
- `ErrorException`: si falta `offsets`, está vacía o aparece una clave no admitida.
- `ErrorException`: si la excentricidad resultante sale de [0, 1).

# Ejemplo
```julia
chief = Dict("SMA" => 7000.0, "ECC" => 0.0, "INC" => 90.0, "RAAN" => 0.0, "AOP" => 0.0, "TA" => 0.0)
layout = Dict("type" => "element_offsets", "offsets" => [Dict(), Dict("dM" => 0.1, "dINC" => 0.02)])
element_offsets_layout(chief, layout, ctx)   # [jefe, miembro con TA = 0.1 e INC = 90.02]
```
"""
function element_offsets_layout(chief::AbstractDict, layout::AbstractDict, ctx)
    check_keys(layout, ("type", "offsets"); who = "element_offsets")

    haskey(layout, "offsets") || error("'element_offsets': missing required parameter 'offsets'.")
    offsets = layout["offsets"]
    isempty(offsets) && error("'element_offsets': 'offsets' is empty.")

    M0 = true_to_mean_anomaly_deg(chief["TA"], chief["ECC"])     # se desplaza en M, no en TA

    members = Dict{String,Any}[]
    for (k, o) in enumerate(offsets)
        check_keys(o, OFFSET_KEYS; who = "element_offsets (offset $k)")

        el = Dict{String,Any}(chief)                             # copia: no se modifica el jefe
        el["SMA"]  = chief["SMA"]  + getopt(o, "dSMA", 0.0)
        el["ECC"]  = chief["ECC"]  + getopt(o, "dECC", 0.0)
        el["INC"]  = chief["INC"]  + getopt(o, "dINC", 0.0)
        el["RAAN"] = mod(chief["RAAN"] + getopt(o, "dRAAN", 0.0), 360.0)
        el["AOP"]  = mod(chief["AOP"]  + getopt(o, "dAOP", 0.0), 360.0)

        0 <= el["ECC"] < 1 || error("'element_offsets' (offset $k): resulting ECC = $(el["ECC"]) is outside [0, 1).")

        # TA nueva con la ECC ya modificada; sin dM ni dECC el miembro conserva la TA del jefe
        el["TA"] = mean_to_true_anomaly_deg(M0 + getopt(o, "dM", 0.0), el["ECC"])
        push!(members, el)
    end
    return members
end

register_layout!("element_offsets", element_offsets_layout)
