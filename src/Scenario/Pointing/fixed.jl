# Dirección de apuntado inercial fija (ECI), como el obs_dir actual de tus JSON.
#
#   "pointing": {"type": "fixed", "dir": [0, 1, 0]}
#
# Se normaliza aquí (load_config también normaliza, pero así el JSON queda limpio).

function fixed_pointing(spec::AbstractDict, member::AbstractDict, ctx)
    haskey(spec, "dir") || error("Pointing 'fixed' needs 'dir' = [x, y, z].")
    v = Float64.(spec["dir"])
    length(v) == 3 || error("Pointing 'fixed': 'dir' must have 3 components.")
    n = norm(v)
    n > 0 || error("Pointing 'fixed': 'dir' has zero norm.")
    return v ./ n
end

register_pointing!("fixed", fixed_pointing)
