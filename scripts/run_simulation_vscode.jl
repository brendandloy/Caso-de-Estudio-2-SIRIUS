# run_simulation_vscode.jl

if !isdefined(Main, :run_simulation)
    include(joinpath(@__DIR__, "main.jl"))
end

run_simulation()