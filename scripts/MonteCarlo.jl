using JSON3
using CSV
using DataFrames
using Printf
using Statistics
using GLMakie
using Random

# =============================================================================
# JSON UTILITIES
# =============================================================================

function load_json_dict(path::AbstractString)
    isfile(path) || error("Config file not found: $path")
    return JSON3.read(read(path, String), Dict{String,Any})
end

function save_json(path::AbstractString, data)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, JSON3.write(data))
    end
    return nothing
end

function deep_merge!(base::Dict{String,Any}, overrides::Dict{String,Any})
    for (key, value) in overrides
        if value isa Dict{String,Any} &&
           haskey(base, key) &&
           base[key] isa Dict{String,Any}
            deep_merge!(base[key], value)
        else
            base[key] = deepcopy(value)
        end
    end
    return base
end

function resolve_path(path::AbstractString, base_dir::AbstractString)
    return isabspath(path) ? normpath(path) : normpath(joinpath(base_dir, path))
end

function format_decimal_2(value::Real)
    return replace(@sprintf("%.2f", Float64(value)), "." => ",")
end

# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

function create_next_mc_directory(output_root::AbstractString)
    mkpath(output_root)

    n = 1
    while isdir(joinpath(output_root, "MC-$n"))
        n += 1
    end

    mc_dir = joinpath(output_root, "MC-$n")
    mkpath(mc_dir)

    return mc_dir, n
end


# =============================================================================
# UKF SUMMARY EXTRACTION
# =============================================================================

function find_metric_row(df::DataFrame, metric_name::AbstractString)
    idx = findfirst(==(metric_name), String.(df.metric))
    return idx
end

function get_value(df::DataFrame, row_idx, column::Symbol)
    row_idx === nothing && return missing
    column ∈ propertynames(df) || return missing

    value = df[row_idx, column]
    return ismissing(value) ? missing : value
end

function append_summary_row!(
    summary_df::DataFrame;
    run_id::Int,
    seed::Int,
    case_tag::String,
    observers::Vector{String},
    sim_status::String,
    ukf_status::String,
    summary_csv::AbstractString
)
    position_idx = nothing
    velocity_idx = nothing
    conv_df = DataFrame()

    if ukf_status == "ok" && isfile(summary_csv)
        conv_df = CSV.read(summary_csv, DataFrame; delim=';')
        if :metric ∈ propertynames(conv_df)
            position_idx = find_metric_row(conv_df, "position")
            velocity_idx = find_metric_row(conv_df, "velocity")
        end
    end

    mean_dop = get_value(conv_df, position_idx, :mean_dop_km)

    push!(summary_df, (
        run = run_id,
        seed = seed,
        case_tag = case_tag,
        n_observers = length(observers),
        observers = join(observers, ","),
        sim_status = sim_status,
        ukf_status = ukf_status,

        position_converged = get_value(conv_df, position_idx, :converged_within_5_min),
        position_convergence_time_s = get_value(conv_df, position_idx, :convergence_time_s),
        position_post_mean_m = get_value(conv_df, position_idx, :post_convergence_mean),
        position_post_std_m = get_value(conv_df, position_idx, :post_convergence_std),
        position_post_rms_m = get_value(conv_df, position_idx, :post_convergence_rms),

        velocity_converged = get_value(conv_df, velocity_idx, :converged_within_5_min),
        velocity_convergence_time_s = get_value(conv_df, velocity_idx, :convergence_time_s),
        velocity_post_mean_m_s = get_value(conv_df, velocity_idx, :post_convergence_mean),
        velocity_post_std_m_s = get_value(conv_df, velocity_idx, :post_convergence_std),
        velocity_post_rms_m_s = get_value(conv_df, velocity_idx, :post_convergence_rms),

        mean_dop_km = mean_dop
    ))

    return nothing
end

function empty_summary_dataframe()
    return DataFrame(
        run = Int[],
        seed = Int[],
        case_tag = String[],
        n_observers = Int[],
        observers = String[],
        sim_status = String[],
        ukf_status = String[],

        position_converged = Union{Missing,Bool}[],
        position_convergence_time_s = Union{Missing,Float64}[],
        position_post_mean_m = Union{Missing,Float64}[],
        position_post_std_m = Union{Missing,Float64}[],
        position_post_rms_m = Union{Missing,Float64}[],

        velocity_converged = Union{Missing,Bool}[],
        velocity_convergence_time_s = Union{Missing,Float64}[],
        velocity_post_mean_m_s = Union{Missing,Float64}[],
        velocity_post_std_m_s = Union{Missing,Float64}[],
        velocity_post_rms_m_s = Union{Missing,Float64}[],

        mean_dop_km = Union{Missing,Float64}[]
    )
end


# =============================================================================
# MONTE CARLO STATISTICS
# =============================================================================

function finite_values(column)
    values = Float64[]

    for value in column
        if !ismissing(value)
            x = Float64(value)
            isfinite(x) && push!(values, x)
        end
    end

    return values
end

function safe_mean(values::Vector{Float64})
    return isempty(values) ? missing : mean(values)
end

function safe_std(values::Vector{Float64})
    if isempty(values)
        return missing
    elseif length(values) == 1
        return 0.0
    end

    return std(values)
end

function convergence_rate(column)
    values = Bool[]

    for value in column
        !ismissing(value) && push!(values, Bool(value))
    end

    return isempty(values) ? missing : 100.0 * count(identity, values) / length(values)
end

function build_mc_statistics(summary_df::DataFrame)
    stats_df = DataFrame(
        n_observers = Int[],
        n_cases = Int[],
        n_ukf_ok = Int[],

        position_rms_mean_m = Union{Missing,Float64}[],
        position_rms_std_m = Union{Missing,Float64}[],
        velocity_rms_mean_m_s = Union{Missing,Float64}[],
        velocity_rms_std_m_s = Union{Missing,Float64}[],

        position_tconv_mean_s = Union{Missing,Float64}[],
        position_tconv_std_s = Union{Missing,Float64}[],
        velocity_tconv_mean_s = Union{Missing,Float64}[],
        velocity_tconv_std_s = Union{Missing,Float64}[],

        position_convergence_rate_pct = Union{Missing,Float64}[],
        velocity_convergence_rate_pct = Union{Missing,Float64}[]
    )

    for group in groupby(summary_df, :n_observers)
        n_obs = Int(first(group.n_observers))
        n_cases = nrow(group)
        n_ukf_ok = count(==("ok"), group.ukf_status)

        position_rms = finite_values(group.position_post_rms_m)
        velocity_rms = finite_values(group.velocity_post_rms_m_s)

        position_tconv = finite_values(group.position_convergence_time_s)
        velocity_tconv = finite_values(group.velocity_convergence_time_s)

        push!(stats_df, (
            n_observers = n_obs,
            n_cases = n_cases,
            n_ukf_ok = n_ukf_ok,

            position_rms_mean_m = safe_mean(position_rms),
            position_rms_std_m = safe_std(position_rms),
            velocity_rms_mean_m_s = safe_mean(velocity_rms),
            velocity_rms_std_m_s = safe_std(velocity_rms),

            position_tconv_mean_s = safe_mean(position_tconv),
            position_tconv_std_s = safe_std(position_tconv),
            velocity_tconv_mean_s = safe_mean(velocity_tconv),
            velocity_tconv_std_s = safe_std(velocity_tconv),

            position_convergence_rate_pct = convergence_rate(group.position_converged),
            velocity_convergence_rate_pct = convergence_rate(group.velocity_converged)
        ))
    end

    sort!(stats_df, :n_observers)

    return stats_df
end

function extract_metric_values(stats_df::DataFrame, metric_col::Symbol)
    x = Int[]
    y = Float64[]

    for row in eachrow(stats_df)
        value = row[metric_col]

        if ismissing(value)
            continue
        end

        value_float = Float64(value)

        if isfinite(value_float)
            push!(x, Int(row.n_observers))
            push!(y, value_float)
        end
    end

    return x, y
end

function add_metric_series!(ax, stats_df::DataFrame, metric_col::Symbol)
    x, y = extract_metric_values(stats_df, metric_col)

    isempty(x) && return nothing

    lines!(ax, x, y)
    scatter!(ax, x, y)

    return nothing
end

function save_observer_effectiveness_plot(stats_df::DataFrame, output_dir::AbstractString)
    n_obs = Int.(stats_df.n_observers)

    fig = Figure(size=(1200, 900))

    ax1 = Axis(
        fig[1, 1],
        title="RMS medio de posición tras la convergencia",
        xlabel="Número de observadores",
        ylabel="RMS medio [m]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)   
    )

    ax2 = Axis(
        fig[1, 2],
        title="RMS medio de velocidad tras la convergencia",
        xlabel="Número de observadores",
        ylabel="RMS medio [m/s]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    ax3 = Axis(
        fig[2, 1],
        title="Tiempo medio de convergencia en posición",
        xlabel="Número de observadores",
        ylabel="Tiempo medio [s]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    ax4 = Axis(
        fig[2, 2],
        title="Tiempo medio de convergencia en velocidad",
        xlabel="Número de observadores",
        ylabel="Tiempo medio [s]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    add_metric_series!(
        ax1,
        stats_df,
        :position_rms_mean_m
    )

    add_metric_series!(
        ax2,
        stats_df,
        :velocity_rms_mean_m_s
    )

    add_metric_series!(
        ax3,
        stats_df,
        :position_tconv_mean_s
    )

    add_metric_series!(
        ax4,
        stats_df,
        :velocity_tconv_mean_s
    )

    for ax in (ax1, ax2, ax3, ax4)
        ax.xticks = n_obs
    end

    save(joinpath(output_dir, "mc_observer_effectiveness.png"), fig)

    return nothing
end

function save_observer_dispersion_plot(stats_df::DataFrame, output_dir::AbstractString)
    n_obs = Int.(stats_df.n_observers)

    fig = Figure(size=(1200, 900))

    ax1 = Axis(
        fig[1, 1],
        title="Dispersión del RMS del error de posición",
        xlabel="Número de observadores",
        ylabel="Desviación estándar [m]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    ax2 = Axis(
        fig[1, 2],
        title="Dispersión del RMS del error de velocidad",
        xlabel="Número de observadores",
        ylabel="Desviación estándar [m/s]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    ax3 = Axis(
        fig[2, 1],
        title="Dispersión del tiempo de convergencia de posición",
        xlabel="Número de observadores",
        ylabel="Desviación estándar [s]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    ax4 = Axis(
        fig[2, 2],
        title="Dispersión del tiempo de convergencia de velocidad",
        xlabel="Número de observadores",
        ylabel="Desviación estándar [s]",
        xtickformat = values -> string.(Int.(round.(values))),
        ytickformat = values -> format_decimal_2.(values)
    )

    add_metric_series!(
        ax1,
        stats_df,
        :position_rms_std_m
    )

    add_metric_series!(
        ax2,
        stats_df,
        :velocity_rms_std_m_s
    )

    add_metric_series!(
        ax3,
        stats_df,
        :position_tconv_std_s
    )

    add_metric_series!(
        ax4,
        stats_df,
        :velocity_tconv_std_s
    )

    for ax in (ax1, ax2, ax3, ax4)
        ax.xticks = n_obs
    end

    save(joinpath(output_dir, "mc_observer_dispersion.png"), fig)

    return nothing
end

function save_convergence_rate_plot(stats_df::DataFrame, output_dir::AbstractString)
    x = Int.(stats_df.n_observers)

    position_rate = Float64[]
    velocity_rate = Float64[]
    x_position = Int[]
    x_velocity = Int[]

    for row in eachrow(stats_df)
        if !ismissing(row.position_convergence_rate_pct)
            push!(x_position, Int(row.n_observers))
            push!(position_rate, Float64(row.position_convergence_rate_pct))
        end

        if !ismissing(row.velocity_convergence_rate_pct)
            push!(x_velocity, Int(row.n_observers))
            push!(velocity_rate, Float64(row.velocity_convergence_rate_pct))
        end
    end

    fig = Figure(size=(800, 500))
    ax = Axis(
        fig[1, 1],
        title="Tasa de convergencia del UKF",
        xlabel="Número de observadores",
        ylabel="Ejecuciones convergidas [%]",
        limits=(nothing, (0.0, 100.0))
    )

    if !isempty(x_position)
        lines!(ax, x_position, position_rate, label="Posición")
        scatter!(ax, x_position, position_rate)
    end

    if !isempty(x_velocity)
        lines!(ax, x_velocity, velocity_rate, label="Velocidad")
        scatter!(ax, x_velocity, velocity_rate)
    end

    ax.xticks = x
    axislegend(ax)

    save(joinpath(output_dir, "mc_convergence_rate.png"), fig)

    return nothing
end

function save_mc_statistics(summary_df::DataFrame, output_dir::AbstractString; make_plots::Bool=true)
    stats_df = build_mc_statistics(summary_df)

    CSV.write(joinpath(output_dir, "mc_statistics.csv"), stats_df; delim=';')

    if make_plots
        save_observer_effectiveness_plot(stats_df, output_dir)
        save_observer_dispersion_plot(stats_df, output_dir)
        save_convergence_rate_plot(stats_df, output_dir)
    end

    return stats_df
end

# =============================================================================
# VALIDATION
# =============================================================================

function validate_observer_cases(sim_cfg::Dict{String,Any}, observer_cases)
    haskey(sim_cfg, "observers") || error("Simulator base config does not contain 'observers'.")

    available = Set(String.(collect(keys(sim_cfg["observers"]))))

    for case in observer_cases
        observers = String.(case["observers"])
        missing_observers = setdiff(Set(observers), available)

        isempty(missing_observers) || error(
            "Observer case contains observers not present in simulator config: " *
            join(sort(collect(missing_observers)), ", ")
        )
    end

    return nothing
end

# =============================================================================
# SIMULATOR AND KALMAN LOADING
# =============================================================================

mc_boot_config_path = isempty(ARGS) ?
    joinpath(@__DIR__, "..", "config", "mc_config.json") :
    abspath(ARGS[1])

mc_boot_cfg = load_json_dict(mc_boot_config_path)
mc_boot_dir = dirname(mc_boot_config_path)

sim_script_boot = resolve_path(
    String(mc_boot_cfg["sim_script"]),
    mc_boot_dir
)

ukf_script_boot = resolve_path(
    String(mc_boot_cfg["ukf_script"]),
    mc_boot_dir
)

if !isdefined(Main, :run_simulation)
    include(sim_script_boot)
end

if !isdefined(Main, :run_kalman)
    include(ukf_script_boot)
end

# =============================================================================
# MAIN
# =============================================================================

function main()
    mc_config_path = isempty(ARGS) ?
        joinpath(@__DIR__, "..", "config", "mc_config.json") :
        abspath(ARGS[1])

    mc_cfg = load_json_dict(mc_config_path)
    mc_cfg_dir = dirname(mc_config_path)

    n_runs = Int(mc_cfg["n_runs"])
    n_runs > 0 || error("'n_runs' must be greater than zero.")

    base_seed = Int(mc_cfg["base_seed"])

    sim_base_config_path = resolve_path(String(mc_cfg["sim_base_config"]), mc_cfg_dir)
    ukf_base_config_path = resolve_path(String(mc_cfg["ukf_base_config"]), mc_cfg_dir)
    output_root = resolve_path(String(mc_cfg["output_root"]), mc_cfg_dir)
    make_plots = Bool(get(mc_cfg, "make_plots", true))

    sim_base_cfg = load_json_dict(sim_base_config_path)
    ukf_base_cfg = load_json_dict(ukf_base_config_path)

    observer_cases = mc_cfg["observer_cases"]
    isempty(observer_cases) && error("'observer_cases' cannot be empty.")

    validate_observer_cases(sim_base_cfg, observer_cases)

    sim_overrides = haskey(mc_cfg, "sim_overrides") ?
        Dict{String,Any}(mc_cfg["sim_overrides"]) :
        Dict{String,Any}()

    ukf_overrides = haskey(mc_cfg, "ukf_overrides") ?
        Dict{String,Any}(mc_cfg["ukf_overrides"]) :
        Dict{String,Any}()

    mc_dir, mc_number = create_next_mc_directory(output_root)

    cp(mc_config_path, joinpath(mc_dir, "mc_config.json"); force=true)
    cp(sim_base_config_path, joinpath(mc_dir, "sim_base_config.json"); force=true)
    cp(ukf_base_config_path, joinpath(mc_dir, "kalman_base_config.json"); force=true)


    summary_df = empty_summary_dataframe()
    summary_path = joinpath(mc_dir, "mc_summary.csv")

    max_observers = length(sim_base_cfg["observers"])

    println("============================================================")
    println("Monte Carlo batch: MC-$mc_number")
    println("Runs: $n_runs")
    println("Maximum simulated observers: $max_observers")
    println("Output: $mc_dir")
    println("============================================================")

    for run_id in 1:n_runs
        run_seed = base_seed + run_id - 1
        run_tag = @sprintf("run-%03d", run_id)

        println("\n[$run_tag] Starting simulator...")

        sim_output = joinpath(
            mc_dir,
            @sprintf("sim_%s_obsmax-%02d", run_tag, max_observers)
        )
        mkpath(sim_output)

        sim_cfg = deepcopy(sim_base_cfg)
        deep_merge!(sim_cfg, sim_overrides)

        sim_cfg_path = joinpath(sim_output, "sim_config.json")
        save_json(sim_cfg_path, sim_cfg)

        sim_status = "ok"

        try
            Random.seed!(run_seed)

            run_simulation(
                sim_cfg_path,
                sim_output
            )
        catch err
            sim_status = "failed"
            @warn "Simulator failed for $run_tag." exception=(err, catch_backtrace())
        end

        for case in observer_cases
            case_tag = String(case["tag"])
            observers = String.(case["observers"])
            n_obs = length(observers)

            ukf_output = joinpath(
                mc_dir,
                @sprintf("ukf_%s_%s", run_tag, case_tag)
            )
            mkpath(ukf_output)

            ukf_cfg = deepcopy(ukf_base_cfg)
            deep_merge!(ukf_cfg, ukf_overrides)
            ukf_cfg["filter_observers"] = observers

            # Generate Kalman plots only for the first Monte Carlo run.
            ukf_cfg["plot_results"] = run_id == 1

            ukf_cfg_path = joinpath(ukf_output, "kalman_config.json")
            save_json(ukf_cfg_path, ukf_cfg)

            summary_csv = joinpath(ukf_output, "summary.csv")

            ukf_status = "not_run"

            if sim_status == "ok"
                println(
                    @sprintf(
                        "[%s | %s] Running UKF with %d observers...",
                        run_tag,
                        case_tag,
                        n_obs
                    )
                )

                try
                    run_kalman(
                        sim_output,
                        ukf_cfg_path,
                        ukf_output
                    )

                    ukf_status = "ok"
                catch err
                    ukf_status = "failed"
                    @warn "UKF failed for $run_tag / $case_tag." exception=(err, catch_backtrace())
                end
            end

            append_summary_row!(
                summary_df;
                run_id=run_id,
                seed=run_seed,
                case_tag=case_tag,
                observers=observers,
                sim_status=sim_status,
                ukf_status=ukf_status,
                summary_csv=summary_csv
            )

            # The summary is rewritten after every case so partial results survive interruptions.
            CSV.write(summary_path, summary_df; delim=';')
        end
    end

    stats_df = save_mc_statistics(summary_df, mc_dir; make_plots=make_plots)

    println("\n============================================================")
    println("Monte Carlo batch completed.")
    println("Run summary: $summary_path")
    println("Statistics: ", joinpath(mc_dir, "mc_statistics.csv"))

    if make_plots
        println("Observer effectiveness plot: ", joinpath(mc_dir, "mc_observer_effectiveness.png"))
        println("Observer dispersion plot: ", joinpath(mc_dir, "mc_observer_dispersion.png"))
        println("Convergence rate plot: ", joinpath(mc_dir, "mc_convergence_rate.png"))
    end

    println("============================================================")

    return nothing
end

main()
