# =============================================================================
# test_camera.jl
#
# Pruebas del módulo CameraSim con casos calculados a mano: conversión RA/Dec ->
# vector y vuelta, proyección pinhole (centro, 100 px a la derecha y abajo, detrás
# de la cámara, cámara girada), PSF gaussiana (suma 1, mitad del máximo a FWHM/2,
# conserva el flujo), escala de magnitudes y filtro del cono del FOV.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_camera.jl
#
# Carga Types, Constants y Camera_sim (FITSIO, Images... tardan en cargar la
# primera vez). No lee el catálogo Hipparcos ni escribe imágenes.
#
# Funciones: ninguna (solo pruebas)
# =============================================================================

using Test
using LinearAlgebra
using StaticArrays
using DataFrames

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)     || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants) || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :CameraSim) || include(joinpath(ROOT, "src", "Camera_sim.jl"))

using .Types, .Constants
const C = CameraSim                    # acceso corto a las funciones internas

@testset "Camera" begin

    @testset "RA/Dec <-> vector" begin
        for (ra, dec) in ((0.0, 0.0), (123.4, -56.7), (359.0, 89.0))
            @test collect(C.vec_to_radec(C.radec_to_vec(ra, dec))) ≈ [ra, dec]
        end
    end

    @testset "Proyección pinhole" begin
        # Imagen de 1001 × 1001: centro en el píxel 501; una dirección con tan α = 100 / f cae 100 px más allá
        nx = ny = 1001
        fov = 20.0
        f = 0.5 * nx / tand(fov / 2)                     # focal en píxeles
        t = 100 / f
        dirs = [SVector(0.0, 0, 1), SVector(t, 0, 1), SVector(0, t, 1), SVector(0.0, 0, -1)]
        radec = [C.vec_to_radec(normalize(d)) for d in dirs]
        x = zeros(Int, 4); y = zeros(Int, 4)

        # Cámara alineada con ECI: x_c = x, y_c = y (abajo en la imagen), boresight = z
        C.radec_to_pixels!(x, y, first.(radec), last.(radec), SMatrix{3,3,Float64,9}(I), fov, nx, ny)
        @test (x[1], y[1]) == (501, 501)                 # boresight
        @test (x[2], y[2]) == (601, 501)                 # +x_c: a la derecha
        @test (x[3], y[3]) == (501, 601)                 # +y_c: hacia abajo
        @test x[4] == C.OFF_IMAGE                        # detrás de la cámara

        # Cámara girada con el boresight hacia +x ECI (filas de C_ci = ejes cámara en ECI)
        C_ci = @SMatrix [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
        @test C.boresight_eci(C_ci) == SVector(1.0, 0, 0)
        xb = zeros(Int, 1); yb = zeros(Int, 1)
        C.radec_to_pixels!(xb, yb, [0.0], [0.0], C_ci, fov, nx, ny)
        @test (xb[1], yb[1]) == (501, 501)
    end

    @testset "PSF gaussiana y flujo" begin
        # Núcleo normalizado; a FWHM/2 = 2 px del centro vale la mitad del máximo
        K, rmax = C.gaussian_kernel(4.0)
        @test sum(K) ≈ 1
        @test K[rmax + 1, rmax + 3] / K[rmax + 1, rmax + 1] ≈ 0.5

        # Una fuente en el centro de la imagen conserva todo su flujo (FWHM = 3 px)
        img = zeros(101, 101)
        C.render_psf_gaussian!(img, (x = [51], y = [51], flux = [1000.0]), 10.0, 3 * 10.0 / 101)
        @test sum(img) ≈ 1000.0

        # 5 magnitudes más débil = 100 veces menos flujo
        @test C.mag_to_flux(13.0, 8.0, 20000.0) ≈ 200.0
    end

    @testset "Cono del FOV" begin
        # Boresight +x y FOV de 20°: una estrella a 0° entra y otra a 15° (> 10°) no
        stars = DataFrame(ra = [0.0, 0.0], dec = [0.0, 15.0], mag = [5.0, 5.0])
        kept = C.filter_cone(stars, SVector(1.0, 0, 0), 20.0)
        @test nrow(kept) == 1 && kept.dec[1] == 0.0
    end
end
