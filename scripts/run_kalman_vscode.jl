# run_kalman_vscode.jl

if !isdefined(Main, :run_kalman)
    include(joinpath(@__DIR__, "Kalman_RSO.jl"))
end

run_kalman()