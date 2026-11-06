# =============================================================================
# Mount.jl
#
# Montaje de la cámara en el satélite observador (parámetros extrínsecos de la
# cámara respecto a los ejes cuerpo) y orientación de la cámara respecto a ECI.
# Se incluye desde Attitude.jl.
#
#   "camera_mount": {"boresight": [0, 0, 1], "up": [0, -1, 0]}
#
# - "boresight": dirección hacia la que mira la cámara, en ejes cuerpo del satélite (eje z_c).
# - "up": dirección del satélite (ejes cuerpo) que sale hacia arriba en la imagen; fija el
#   giro de la cámara alrededor del boresight. No puede ser paralela a "boresight".
# - Sin bloque se usa DEFAULT_CAMERA_MOUNT (C_cb = I) y se avisa.
#
# Ejes cámara (dextrógiros, como en el modelo pinhole): z_c = boresight, x_c = derecha de la
# imagen (columnas crecientes), y_c = abajo en la imagen (filas crecientes).
#
# Modelo (T6: filas de una DCM = ejes del sistema de llegada; composición C_ca = C_cb C_ba):
#   z_c = b̂
#   y_c = -(u - (u·z_c) z_c) / |u - (u·z_c) z_c|      ("arriba" proyectado en el plano del sensor)
#   x_c = y_c × z_c
#   C_cb = [x_cᵀ; y_cᵀ; z_cᵀ]   (cuerpo -> cámara)          C_ci = C_cb · C_bi = C_cb · C_ibᵀ
#
# Funciones: camera_mount_block, camera_mount_dcm, parse_camera_mount, camera_dcm_history
# =============================================================================

# Montaje por defecto: boresight = +z cuerpo y arriba de la imagen = -y cuerpo, es decir, C_cb = I
const DEFAULT_CAMERA_MOUNT = Dict{String,Any}("boresight" => [0, 0, 1], "up" => [0, -1, 0])

# Claves admitidas en el bloque "camera_mount"
const CAMERA_MOUNT_KEYS = ("boresight", "up")

"""
    camera_mount_block(entry, name) -> Dict{String,Any}

Devuelve el bloque `camera_mount` de un observador (del JSON del simulador o de
un grupo de la spec). Si no existe, devuelve una copia de `DEFAULT_CAMERA_MOUNT`
e imprime un aviso para que quede claro qué montaje se usa.

# Argumentos
- `entry::AbstractDict`: entrada del observador (o grupo de observadores de la spec).
- `name::AbstractString`: nombre del observador, para los mensajes.

# Devuelve
- `Dict{String,Any}`: copia del bloque `camera_mount` (sin validar sus parámetros).

# Excepciones
- `ErrorException`: si `camera_mount` existe pero no es un bloque `{...}`.

# Ejemplo
```julia
camera_mount_block(Dict("SMA" => 7000), "cam1")
# aviso + Dict("boresight" => [0, 0, 1], "up" => [0, -1, 0])
```
"""
function camera_mount_block(entry, name::AbstractString)
    e = todict(entry)
    if !haskey(e, "camera_mount")
        println("Observer '$name': no 'camera_mount' block; assuming boresight = body +z ",
                "and image up = body -y (camera axes = body axes).")
        return deepcopy(DEFAULT_CAMERA_MOUNT)
    end
    e["camera_mount"] isa Dict{String,Any} || error("Observer '$name': 'camera_mount' must be a block {...}.")
    return e["camera_mount"]
end

"""
    camera_mount_dcm(boresight, up) -> SMatrix{3,3,Float64,9}

DCM C_cb (cuerpo -> cámara) a partir del boresight y de la dirección "arriba" de
la imagen, ambos en ejes cuerpo. "Arriba" se proyecta sobre el plano
perpendicular al boresight, así que basta con que no sea paralelo a él.

# Argumentos
- `boresight::AbstractVector`: eje de visión en ejes cuerpo (no hace falta que sea unitario).
- `up::AbstractVector`: dirección del cuerpo que sale arriba en la imagen.

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_cb, ortogonal y con det = +1 (filas = ejes cámara en ejes cuerpo).

# Excepciones
- `ErrorException`: si algún vector tiene norma cero o son paralelos.

# Ejemplo
```julia
camera_mount_dcm([0, 0, 1], [0, -1, 0])     # identidad
camera_mount_dcm([1, 0, 0], [0, 0, 1])      # cámara mirando a +x cuerpo, con +z cuerpo arriba
```
"""
function camera_mount_dcm(boresight::AbstractVector, up::AbstractVector)
    b = Vec3(boresight)
    u = Vec3(up)
    (norm(b) > 0 && norm(u) > 0) || error("camera_mount: 'boresight' and 'up' must have non-zero norm.")

    z_c = b / norm(b)                          # eje de visión
    u_perp = u - dot(u, z_c) * z_c             # "arriba" proyectado en el plano del sensor
    norm(u_perp) > 1e-8 * norm(u) || error("camera_mount: 'up' is parallel to 'boresight'.")
    y_c = -u_perp / norm(u_perp)               # y_c apunta hacia abajo en la imagen
    x_c = cross(y_c, z_c)                      # derecha de la imagen; x_c × y_c = z_c

    return Mat3(vcat(x_c', y_c', z_c'))        # filas = ejes cámara en ejes cuerpo
end

"""
    parse_camera_mount(block, name) -> SMatrix{3,3,Float64,9}

Lee y valida un bloque `camera_mount` y devuelve la DCM C_cb.

# Argumentos
- `block::AbstractDict`: bloque `camera_mount` (de `camera_mount_block`).
- `name::AbstractString`: nombre del observador, para los mensajes de error.

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_cb (cuerpo -> cámara).

# Excepciones
- `ErrorException`: si hay claves no admitidas, falta `boresight` o `up`, no tienen 3
  componentes, tienen norma cero o son paralelos.

# Ejemplo
```julia
parse_camera_mount(Dict{String,Any}("boresight" => [1, 0, 0], "up" => [0, 0, 1]), "cam1")
```
"""
function parse_camera_mount(block, name::AbstractString)
    who = "Observer '$name' (camera_mount)"
    b = todict(block)
    b isa Dict{String,Any} || error("$who: must be a block {...}.")
    check_keys(b, CAMERA_MOUNT_KEYS; who = who)

    bs = getunit(b, "boresight", who)
    up = getunit(b, "up", who)
    # Se repite el error de camera_mount_dcm con el nombre del observador delante
    return try
        camera_mount_dcm(bs, up)
    catch err
        err isa ErrorException ? error("$who: " * err.msg) : rethrow()
    end
end

"""
    camera_dcm_history(C_cb, att) -> Vector{SMatrix{3,3,Float64,9}}

Orientación de la cámara respecto a ECI en cada instante: C_ci = C_cb · C_bi,
con C_bi = DCM_b2iᵀ (T6: composición de DCM y inversa = transpuesta).

# Argumentos
- `C_cb::AbstractMatrix`: montaje de la cámara (cuerpo -> cámara).
- `att::AttitudeData`: actitud del satélite observador (DCM_b2i = C_ib por instante).

# Devuelve
- `Vector{SMatrix{3,3,Float64,9}}`: C_ci (ECI -> cámara) por instante. Su fila 3 es el
  boresight en ECI.

# Excepciones
- Ninguna.

# Ejemplo
```julia
C_ci = camera_dcm_history(C_cb, compute_attitude(spec, traj, sun_dir))
C_ci[1][3, :]        # boresight en ECI en el primer instante
```
"""
function camera_dcm_history(C_cb::AbstractMatrix, att::AttitudeData)
    C = Mat3(C_cb)
    return [C * D' for D in att.DCM_b2i]
end
