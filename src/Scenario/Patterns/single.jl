# Patrón más simple: un único objeto con una órbita.
#
#   "cam1": {"pattern": "single", "orbit": {"type": "leo", "altitude_km": 621.86, "INC": 90}, ...}
#
# Si se omite "pattern", el Builder usa "single" por defecto.

function single_pattern(group::AbstractDict, ctx)
    haskey(group, "orbit") || error("Pattern 'single' needs an 'orbit' block.")
    return [build_orbit(group["orbit"], ctx)]
end

register_pattern!("single", single_pattern)
