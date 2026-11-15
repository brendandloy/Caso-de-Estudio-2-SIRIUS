# =============================================================================
# test_camera.jl
#
# Pruebas del módulo CameraSim con casos calculados a mano: conversión RA/Dec ->
# vector y vuelta, intrínsecos (lectura, punto principal por defecto, errores, FOV del
# sensor), proyección pinhole con K (centro, 100 px a la derecha y abajo, detrás de la
# cámara, cámara girada, punto principal desplazado, skew), selección por sensor (FOV
# rectangular: esquina dentro, fuera del sensor no), PSF gaussiana (suma 1, mitad del
# máximo a FWHM/2, elíptica con f_x ≠ f_y, conserva el flujo), escala de magnitudes y
# filtro del cono.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_camera.jl
#
# Carga Types, Constants, Parameters, JSON_utils y Camera_sim (FITSIO, Images... tardan
# en cargar la primera vez). No lee el catálogo Hipparcos ni escribe imágenes.
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
isdefined(Main, :Parameters) || include(joinpath(ROOT, "src", "Parameters.jl"))
isdefined(Main, :JSONUtils) || include(joinpath(ROOT, "src", "JSON_utils.jl"))
isdefined(Main, :CameraSim) || include(joinpath(ROOT, "src", "Camera_sim.jl"))

using .Types, .Constants
const C = CameraSim                    # acceso corto a las funciones internas

@testset "Camera" begin

    @testset "RA/Dec <-> vector" begin
        for (ra, dec) in ((0.0, 0.0), (123.4, -56.7), (359.0, 89.0))
            @test collect(C.vec_to_radec(C.radec_to_vec(ra, dec))) ≈ [ra, dec]
        end
    end

    @testset "Intrínsecos" begin
        b = Dict("resolution_x" => 1024, "resolution_y" => 768, "fx_px" => 2000.0, "fy_px" => 2000.0)
        cam = C.camera_intrinsics(b)
        @test (cam.u0_px, cam.v0_px, cam.skew_px) == (512.5, 384.5, 0.0)     # centro y skew por defecto
        @test C.camera_intrinsics(merge(b, Dict("u0_px" => nothing))).u0_px == 512.5   # null = centro
        @test C.camera_intrinsics(merge(b, Dict("u0_px" => 500.0))).u0_px == 500.0
        @test_throws ErrorException C.camera_intrinsics(Dict("resolution_x" => 1024, "resolution_y" => 768, "fx_px" => 2000.0))
        @test_throws ErrorException C.camera_intrinsics(merge(b, Dict("fy_px" => -1.0)))
        @test_throws ErrorException C.camera_intrinsics(merge(b, Dict("resolution_x" => 10.5)))

        # FOV del sensor: con el punto principal centrado, 2 atan(N/2f); 27° con f = 512/tan(13.5°)
        f27 = 512 / tand(13.5)
        fx, fy = C.sensor_fov_deg(C.CameraIntrinsics(1024, 1024, f27, f27, 512.5, 512.5, 0.0))
        @test fx ≈ 27.0 && fy ≈ 27.0
        @test C.sensor_fov_deg(cam)[2] ≈ 2atand(384 / 2000)                   # vertical, 768 px
        # Cono de prefiltro: esquina del sensor (√2 · 512 px de un sensor cuadrado) por el margen
        cam27 = C.CameraIntrinsics(1024, 1024, f27, f27, 512.5, 512.5, 0.0)
        @test C.prefilter_cone_deg(cam27) ≈ 2 * C.FOV_PREFILTER_MARGIN * atand(sqrt(2) * 512 / f27)
    end

    @testset "Proyección pinhole" begin
        # Imagen de 1001 × 1001: centro en el píxel 501; una dirección con tan α = 100 / f cae 100 px más allá
        nx = ny = 1001
        f = 0.5 * nx / tand(10.0)                        # focal en píxeles (FOV de 20°)
        cam = C.CameraIntrinsics(nx, ny, f, f, 501.0, 501.0, 0.0)
        t = 100 / f
        dirs = [SVector(0.0, 0, 1), SVector(t, 0, 1), SVector(0, t, 1), SVector(0.0, 0, -1)]
        radec = [C.vec_to_radec(normalize(d)) for d in dirs]
        x = zeros(Int, 4); y = zeros(Int, 4)

        # Cámara alineada con ECI: x_c = x, y_c = y (abajo en la imagen), boresight = z
        C.radec_to_pixels!(x, y, first.(radec), last.(radec), SMatrix{3,3,Float64,9}(I), cam)
        @test (x[1], y[1]) == (501, 501)                 # boresight
        @test (x[2], y[2]) == (601, 501)                 # +x_c: a la derecha
        @test (x[3], y[3]) == (501, 601)                 # +y_c: hacia abajo
        @test x[4] == C.OFF_IMAGE                        # detrás de la cámara

        # Punto principal desplazado (+30, -20 px): toda la imagen se desplaza lo mismo
        cam_pp = C.CameraIntrinsics(nx, ny, f, f, 531.0, 481.0, 0.0)
        C.radec_to_pixels!(x, y, first.(radec), last.(radec), SMatrix{3,3,Float64,9}(I), cam_pp)
        @test (x[1], y[1]) == (531, 481) && (x[2], y[2]) == (631, 481)

        # Focales distintas y skew: u = f_x x' + α y' + u0, v = f_y y' + v0
        cam_k = C.CameraIntrinsics(nx, ny, f, 2f, 501.0, 501.0, 50.0)
        C.radec_to_pixels!(x, y, first.(radec), last.(radec), SMatrix{3,3,Float64,9}(I), cam_k)
        @test (x[3], y[3]) == (501 + round(Int, 50 * t), 701)    # y' = t: 200 px abajo y α t a la derecha

        # Cámara girada con el boresight hacia +x ECI (filas de C_ci = ejes cámara en ECI)
        C_ci = @SMatrix [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
        @test C.boresight_eci(C_ci) == SVector(1.0, 0, 0)
        xb = zeros(Int, 1); yb = zeros(Int, 1)
        C.radec_to_pixels!(xb, yb, [0.0], [0.0], C_ci, cam)
        @test (xb[1], yb[1]) == (501, 501)
    end

    @testset "Selección por sensor (FOV rectangular)" begin
        # Sensor 1001 × 1001 con FOV de 20° (semiángulo 10°). Una estrella cerca de la esquina,
        # a 12.6° del boresight, entra (fuera del antiguo cono de 10°); otra a 11° en x, no
        nx = ny = 1001
        f = 0.5 * nx / tand(10.0)
        cam = C.CameraIntrinsics(nx, ny, f, f, 501.0, 501.0, 0.0)
        corner = normalize(SVector(450 / f, 450 / f, 1.0))       # 450 px en x e y: dentro del sensor
        outside = normalize(SVector(tand(11.0), 0.0, 1.0))       # 11° en x: fuera del sensor
        radec = [C.vec_to_radec(corner), C.vec_to_radec(outside)]
        stars = DataFrame(ra = first.(radec), dec = last.(radec), mag = [5.0, 5.0])
        config = (intrinsics = cam, mag0 = 8.0, flux0 = 20000.0)
        kept, render = C.sources_in_view(stars, SMatrix{3,3,Float64,9}(I), config)
        @test acosd(corner[3]) > 10                               # fuera del antiguo cono
        @test nrow(kept) == 1 && (render.x[1], render.y[1]) == (951, 951)
        @test C.in_sensor(1, 1001, cam) && !C.in_sensor(0, 500, cam) && !C.in_sensor(C.OFF_IMAGE, 500, cam)
    end

    @testset "PSF gaussiana y flujo" begin
        # Núcleo normalizado; a FWHM/2 = 2 px del centro vale la mitad del máximo
        K, rx, ry = C.gaussian_kernel(4.0, 4.0)
        @test sum(K) ≈ 1 && rx == ry
        @test K[ry + 1, rx + 3] / K[ry + 1, rx + 1] ≈ 0.5

        # Elíptica (FWHM 4 px en x y 8 px en y): mitad del máximo a 2 px en x y a 4 px en y
        K, rx, ry = C.gaussian_kernel(4.0, 8.0)
        @test size(K) == (2ry + 1, 2rx + 1) && ry > rx
        @test K[ry + 1, rx + 3] / K[ry + 1, rx + 1] ≈ 0.5
        @test K[ry + 5, rx + 1] / K[ry + 1, rx + 1] ≈ 0.5

        # Una fuente en el centro de la imagen conserva todo su flujo (FWHM = f · fwhm = 3 px)
        img = zeros(101, 101)
        cam = C.CameraIntrinsics(101, 101, 1000.0, 1000.0, 51.0, 51.0, 0.0)
        C.render_psf_gaussian!(img, (x = [51], y = [51], flux = [1000.0]), cam, rad2deg(3 / 1000))
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
