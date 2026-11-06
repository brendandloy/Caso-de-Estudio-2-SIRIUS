# =============================================================================
# test_photometry.jl
#
# Pruebas del módulo Photometry con una placa de una faceta y casos calculados a
# mano: flujo de Phong con la placa de frente al Sol y al observador, ley 1/d²
# (duplicar la distancia suma 5·log10(2) = 1,505 mag), giro de la normal con la
# DCM, faceta de espaldas y objetivo en la sombra de la Tierra.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_photometry.jl
#
# Carga Types, Constants, Parameters, JSON_utils, Transformations, Geometry,
# Visibility y Photometry. Sol en +x ECI; distancias en km y áreas en km².
#
# Funciones: plate, magnitude
# =============================================================================

using Test
using LinearAlgebra
using StaticArrays

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Parameters)      || include(joinpath(ROOT, "src", "Parameters.jl"))
isdefined(Main, :JSONUtils)       || include(joinpath(ROOT, "src", "JSON_utils.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))
isdefined(Main, :Geometry)        || include(joinpath(ROOT, "src", "Geometry", "Geometry.jl"))
isdefined(Main, :Visibility)      || include(joinpath(ROOT, "src", "Visibility.jl"))
isdefined(Main, :Photometry)      || include(joinpath(ROOT, "src", "Photometry.jl"))

using .Types, .Geometry, .Photometry

const SUN_X = SVector(1.0, 0.0, 0.0)                 # dirección del Sol: +x ECI
const PHOT_CFG = (phi_sun = 455.0, glare_angle_deg = 15.0)
const C_SPEC, C_DIFF, AREA_KM2 = 0.3, 0.5, 1e-6      # placa de 1 m² (1e-6 km²)

# Placa de una faceta con normal `normal` (ejes cuerpo)
plate(normal) = DebrisModel(:plate, Dict{String,Float64}(), 1.0, 0, reshape(collect(normal), 3, 1),
                            [AREA_KM2], AREA_KM2, 0.2, C_SPEC, C_DIFF, 1.2)

"""
    magnitude(debris, d; r_obj, dcm) -> NamedTuple

Resultado de `compute_apparent_magnitude` en un instante, con el observador a una
distancia `d` (km) del objetivo en la dirección +x (hacia el Sol).

Por qué existe: todas las pruebas usan la misma geometría y solo cambian la
distancia, la posición del objetivo o su actitud.

# Argumentos
- `debris::DebrisModel`: objetivo.
- `d::Real`: distancia observador-objetivo (km).
- `r_obj`: posición del objetivo en ECI (km); por defecto, en GEO del lado del Sol.
- `dcm`: DCM cuerpo -> ECI del objetivo; por defecto, la identidad.

# Devuelve
- `NamedTuple` de `compute_apparent_magnitude` (vectores de un elemento).

# Excepciones
- Ninguna.
"""
function magnitude(debris, d; r_obj = SVector(42164.0, 0, 0), dcm = SMatrix{3,3,Float64,9}(I))
    v0 = [zero(SVector{3,Float64})]                       # velocidades: no intervienen
    obs = TrajectoryData([r_obj + SVector(d, 0, 0)], v0, [0.0])
    obj = TrajectoryData([r_obj], v0, [0.0])
    return compute_apparent_magnitude(obs, obj, [SUN_X], debris, PHOT_CFG, [dcm])
end

@testset "Photometry" begin

    @testset "Placa de frente al Sol y al observador" begin
        # cos θ_in = cos θ_v = cos θ_s = 1 -> F = Φ (C_diff + C_spec) A / d²
        d = 1000.0
        res = magnitude(plate((1.0, 0, 0)), d)
        F = PHOT_CFG.phi_sun * (C_DIFF + C_SPEC) * AREA_KM2 / d^2
        @test res.visible[1]
        @test res.Vmag[1] ≈ -12.6988 - 2.5 * log10(F)

        # Ley 1/d²: al duplicar la distancia, la magnitud aumenta 5·log10(2) ≈ 1,505
        @test magnitude(plate((1.0, 0, 0)), 2d).Vmag[1] - res.Vmag[1] ≈ 5 * log10(2)
    end

    @testset "Actitud, faceta de espaldas y sombra" begin
        # Normal +z en ejes cuerpo girada 90° alrededor de y: apunta a +x ECI, mismo resultado
        Ry = @SMatrix [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0]
        @test magnitude(plate((0, 0, 1.0)), 1000.0; dcm = Ry).Vmag[1] ≈ magnitude(plate((1.0, 0, 0)), 1000.0).Vmag[1]

        # Faceta de espaldas al Sol y al observador: no se ve
        res = magnitude(plate((-1.0, 0, 0)), 1000.0)
        @test !res.visible[1] && isinf(res.Vmag[1])

        # Objetivo detrás de la Tierra (en umbra): no se ve
        res = magnitude(plate((1.0, 0, 0)), 1000.0; r_obj = SVector(-42164.0, 0, 0))
        @test res.in_eclipse[1] && isinf(res.Vmag[1])
    end
end
