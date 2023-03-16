using CSV, StatsPlots, Plots, DataFrames

function plot_runtimes(_expnr_c, _expnr_o, _file_name)
    _c_results = load_tsp_results(_expnr_c)
    _o_results = load_tsp_results(_expnr_o)

    data_points = Vector{Float64}[]
    labels = String[]

    for algid in unique(_c_results.algid)
        _rel_results_c = filter(row -> row.algid == algid, _c_results)
        _rel_results_o = filter(row -> row.algid == algid, _o_results)

        _runtimes_c = map(eval ∘ Meta.parse, _rel_results_c.runtimes)
        _runtimes_o = map(eval ∘ Meta.parse, _rel_results_o.runtimes)

        _runtimes_c = reduce(vcat,_runtimes_c)
        _runtimes_o = reduce(vcat,_runtimes_o)

        push!(labels, algid, algid)
        push!(data_points, _runtimes_c, _runtimes_o)
    end

    _p = violin([labels[1]], data_points[1], alpha = .5, color = :red, label = "")

    xlabel!("convex solver")
    ylabel!("time to solve (s)")

    for i in 2:length(labels)
        _color = i % 2 == 0 ? :blue : :red
        _label = ""
        if i == length(labels)
            _label = "objective"
        elseif i == length(labels)-1
            _label = "constraint"
        end

        violin!(_p, [labels[i]], data_points[i], alpha = .5, color = _color, label = _label)
    end

    #return _p
    savefig(_p, "figs/$_file_name")
end

function scale_runtime_plots(_runtimes::DataFrame, systemid::AbstractString, algids = unique(_runtimes.algid), printscale::Bool=false)

    data_points = Vector{Float64}[]
    labels = String[]

    _runtimes = filter(row -> row.systemid == systemid, _runtimes)

    for scale in unique(_runtimes.scale)
        _alg_scales = filter(row -> row.scale == scale, _runtimes)
        for algid in algids
            _alg_scale_runtimes = filter(row -> row.algid == algid, _alg_scales)

            _alg_runtimes = map(eval ∘ Meta.parse, _alg_scale_runtimes.runtimes)

            _alg_all_runtimes = reduce(vcat,_alg_runtimes)

            label = printscale ? "$algid ($scale)" : "$algid"
            push!(labels, label)
            push!(data_points, _alg_all_runtimes)
        end
    end

    violin([labels[1]], data_points[1], color = :red, leg = false)

    xlabel!("convex solver")
    ylabel!("time to solve (s)")

    for i in 2:length(labels)
        _color = i % 2 == 0 ? :blue : :red
        violin!([labels[i]], data_points[i], color = _color, leg = false)
    end
end

function scale_eval_plots(_exp_data::DataFrame, systemid::String, algids = unique(_exp_data.algid), printscale::Bool=false)

    data_points = Vector{Float64}[]
    labels = String[]

    _evals = filter(row -> row.systemid == systemid, _exp_data)

    for scale in unique(_exp_data.scale)
        _alg_scales = filter(row -> row.scale == scale, _evals)
        for algid in algids
            _alg_scale_evaluations = filter(row -> row.algid == algid, _alg_scales)

            _alg_all_evals = reduce(vcat,_alg_scale_evaluations.evaluations)

            label = printscale ? "$algid ($scale)" : "$algid"
            push!(labels, label)
            push!(data_points, _alg_all_evals)
        end
    end

    plotfunc = std(data_points[1]) == 0 ? boxplot : violin # violin does not plot if std is zero, boxplot does
    plotfunc([labels[1]], data_points[1], color = :red, leg = false)

    xlabel!("convex solver")
    ylabel!("number of evaluations")

    for i in 2:length(labels)
        _color = i % 2 == 0 ? :blue : :red
        plotfunc = std(data_points[i]) == 0 ? boxplot! : violin!
        plotfunc([labels[i]], data_points[i], color = _color, leg = false)
    end
end

function scale_eval_plots(_exp_data::DataFrame, systemid::String, to_file::String, algids = unique(_exp_data.algid))
    scale_plots(_exp_data, systemid, to_file, algids, scale_eval_plots)
end

function scale_runtime_plots(_exp_data::DataFrame, systemid::AbstractString, to_file::String, algids = unique(_exp_data.algid))
    scale_plots(_exp_data, systemid, to_file, algids, scale_runtime_plots)
end

function scale_plots(_exp_data::DataFrame, systemid::String, to_file::String, algids, plotting_func::Function)
    plotting_func(_exp_data, systemid, algids)
    savefig(to_file)
end
