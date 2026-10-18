# =============================================================================
# fixed.jl
#
# Apuntado inercial fijo (ECI), como el obs_dir actual de tus JSON. Registra
# el apuntado "fixed".
#
#   "pointing": {"type": "fixed", "dir": [0, 1, 0]}
#
# El vector se normaliza aquí (load_config también lo hace, pero así el JSON queda limpio).
#
# Funciones: fixed_pointing
# =============================================================================

"""
    fixed_pointing(spec, member, ctx) -> Vector{Float64}

Devuelve la dirección de apuntado fija indicada en la spec, normalizada.

Por qué recibe `member` y `ctx` sin usarlos: todos los plugins de apuntado
comparten la misma firma. Otros apuntados (p. ej. hacia un objetivo) sí
necesitan la posición del observador y la lista de objetivos.

# Argumentos
- `spec::AbstractDict`: bloque `pointing`; debe contener `dir` = [x, y, z].
- `member::AbstractDict`: elementos orbitales del observador (no se usan).
- `ctx`: contexto del Builder (no se usa).

# Devuelve
- `Vector{Float64}` unitario de 3 componentes.

# Excepciones
- `ErrorException`: si falta `dir`, no tiene 3 componentes o su norma es cero.

# Ejemplo
```julia
fixed_pointing(Dict("type" => "fixed", "dir" => [0, 2, 0]), member, ctx)     # [0.0, 1.0, 0.0]
```
"""
function fixed_pointing(spec::AbstractDict, member::AbstractDict, ctx)
    haskey(spec, "dir") || error("Pointing 'fixed' needs 'dir' = [x, y, z].")
    v = Float64.(spec["dir"])
    length(v) == 3 || error("Pointing 'fixed': 'dir' must have 3 components.")
    n = norm(v)
    n > 0 || error("Pointing 'fixed': 'dir' has zero norm.")
    return v ./ n
end

register_pointing!("fixed", fixed_pointing)