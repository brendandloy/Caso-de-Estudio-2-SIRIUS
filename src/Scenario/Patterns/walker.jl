# =============================================================================
# walker.jl
#
# Patrón de constelación Walker T/P/F. Registra el patrón "walker".
#
#   "sat": {
#     "pattern":  "walker",
#     "orbit":    {"type": "LEO", "altitude_km": 800, "INC": 55, "RAAN": 0, "TA": 0},
#     "walker":   {"total": 12, "planes": 3, "phasing": 1, "spread": "delta"},
#     "pointing": {"type": "fixed", "dir": [0, 1, 0]}
#   }
#
# - "orbit": órbita del satélite 0 del plano 0 (cualquier tipo registrado). Aporta
#   SMA, ECC, INC, AOP, el RAAN del plano 0 y la posición del primer satélite.
# - "walker": T satélites en P planos (S = T/P por plano) con factor de fase F.
#     spread "delta" (por defecto): planos repartidos en 360° de RAAN.
#     spread "star": planos repartidos en 180° (constelaciones polares tipo Iridium).
#
# Fórmulas (p = 0..P-1 plano, s = 0..S-1 satélite dentro del plano):
#   RAAN(p)   = RAAN_0 + spread·p/P
#   M(p, s)   = M_0 + 360·s/S + 360·F·p/T          (anomalía media; vale también si ECC > 0)
#
# Orden y nombres (los pone el Builder): plano por plano, es decir, el miembro k
# (k = 1..T) es el satélite s = (k-1) mod S del plano p = (k-1) ÷ S.
#
# Funciones: walker_pattern
# =============================================================================

"""
    walker_pattern(group, ctx) -> Vector{Dict{String,Any}}

Genera los `T` satélites de una constelación Walker `T/P/F` a partir de la órbita
del satélite de referencia.

# Argumentos
- `group::AbstractDict`: grupo de la spec con los bloques `orbit` y `walker`
  (`total`, `planes`, `phasing` y, opcionalmente, `spread` = "delta" o "star").
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` de `T` diccionarios de elementos, ordenados plano por plano.

# Excepciones
- `ErrorException`: si falta el bloque `orbit` o `walker`, o aparece una clave no admitida.
- `ErrorException`: si `total`, `planes` o `phasing` no son enteros.
- `ErrorException`: si `total` no es múltiplo de `planes`, o `phasing` no está en 0..planes-1.
- `ErrorException`: si `spread` no es "delta" ni "star".

# Ejemplo
```julia
group = Dict(
    "orbit"  => Dict("type" => "LEO", "altitude_km" => 800, "INC" => 55, "RAAN" => 0, "TA" => 0),
    "walker" => Dict("total" => 6, "planes" => 3, "phasing" => 1))
walker_pattern(group, ctx)      # RAAN 0, 120, 240 (2 por plano); TA 0 y 180, 60 y 240, 120 y 300
```
"""
function walker_pattern(group::AbstractDict, ctx)
    haskey(group, "orbit")  || error("Pattern 'walker' needs an 'orbit' block.")
    haskey(group, "walker") || error("Pattern 'walker' needs a 'walker' block.")

    w = group["walker"]
    check_keys(w, ("total", "planes", "phasing", "spread"); who = "walker")

    T = getint(w, "total")
    P = getint(w, "planes")
    F = getint(w, "phasing")

    P >= 1     || error("'walker': 'planes' must be >= 1 (got $P).")
    T >= P     || error("'walker': 'total' ($T) must be >= 'planes' ($P).")
    T % P == 0 || error("'walker': 'total' ($T) must be a multiple of 'planes' ($P).")
    0 <= F < P || error("'walker': 'phasing' must be in 0..planes-1 (got $F with planes = $P).")

    # Reparto de planos en RAAN: 360° (delta) o 180° (star)
    spread_name = haskey(w, "spread") ? String(w["spread"]) : "delta"
    spread = spread_name == "delta" ? 360.0 :
             spread_name == "star"  ? 180.0 :
             error("'walker': 'spread' must be \"delta\" or \"star\" (got \"$spread_name\").")

    S = T ÷ P                                                   # satélites por plano
    base = build_orbit(group["orbit"], ctx)
    M0 = true_to_mean_anomaly_deg(base["TA"], base["ECC"])

    members = Dict{String,Any}[]
    for p in 0:P-1, s in 0:S-1
        el = Dict{String,Any}(base)                             # copia: no se modifica la órbita base
        el["RAAN"] = mod(base["RAAN"] + spread * p / P, 360.0)

        # El satélite de referencia conserva su TA exacta; el resto se desplaza en anomalía media
        if !(p == 0 && s == 0)
            M = M0 + 360.0 * s / S + 360.0 * F * p / T
            el["TA"] = mean_to_true_anomaly_deg(M, base["ECC"])
        end
        push!(members, el)
    end
    return members
end

register_pattern!("walker", walker_pattern)
