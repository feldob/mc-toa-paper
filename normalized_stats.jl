include("io.jl")

mapto(ts::AbstractVector{<:TspExperimentResultSummary}, f::Function) = map(f, ts)[:,1]

function minmax_criterion(_s::Vector{DataFrame}, criterion::Symbol)
    _min = map(s -> round(minimum(s[:,criterion]), digits = 2), _s)
    _max = map(s -> round(maximum(s[:,criterion]), digits = 2), _s)
    return _min, _max
end

function full_stats(exp_nr::Union{String, Integer}=1, results=load_tsp_results(exp_nr))
    _stats_dir = joinpath("experiment_$exp_nr", "results", "stats")
    _full_stats_file = joinpath(_stats_dir, "fullstats.csv")
    if isfile(_full_stats_file)
        return CSV.read(_full_stats_file, DataFrame)
    end

    _results = load_experimental_summary(exp_nr, results)

    _stats = map(stats, _results)
    _minsize = map(s -> minimum(s[:,:size]), _stats)

    _min_exercises, _max_exercises = minmax_criterion(_stats, :exercises)
    _min_failrate, _max_failrate = minmax_criterion(_stats, :failrate)

    _maxexpersize = map(s -> maximum(s[:,:exercises]), _stats)
    _maxfrpersize = map(s -> maximum(s[:,:failrate]), _stats)

    df = DataFrame(systemid = mapto(_results, systemid),
                    instanceid = mapto(_results, instanceid),
                    scale = mapto(_results, scale),
                    numtests = mapto(_results, numtests),
                    numfuncs = mapto(_results, numfuncs),
                    max_failrate = _max_failrate,
                    max_exercises = _max_exercises,
                    min_failrate = _min_failrate,
                    min_exercises = _min_exercises,
                    algid = mapto(_results, algid),
                    runid = mapto(_results, runid),
                    iterations = mapto(_results, iterations),
                    evaluations = mapto(_results, evaluations),
                    minsize = _minsize,
                    max_fr_per_size = _maxfrpersize,
                    max_ex_per_size = _maxexpersize)

    for ct in categories(_results[1])
        _min_cr, _max_cr = minmax_criterion(_stats, ct)
        df[!, "min_"*string(ct)] = _min_cr
        df[!, "max_"*string(ct)] = _max_cr
    end
    if !isdir(_stats_dir)
        mkpath(_stats_dir)
    end
    CSV.write(_full_stats_file, df)
    return df
end

function dist_to_opt_size!(df::DataFrame, sdf::SubDataFrame)
    minsize_col = sdf[:, :minsize]
    df.dist_to_opt_size = minsize_col .- minimum(minsize_col)
    df.abs_size_reduction = round.((minsize_col ./ minimum(minsize_col) .- 1) .* 100, RoundUp, digits = 2)
    return df
end

function dist_to_best_fr_per_size!(df::DataFrame, sdf::SubDataFrame)
    fr_col = sdf[:, :max_fr_per_size]
    fr_per_size_alltests = sdf[1, :max_failrate] /  sdf[1, :numtests]
    rev_dist = (fr_col .- fr_per_size_alltests) ./ (maximum(fr_col) - fr_per_size_alltests)
    df.dist_to_best_fr_per_size = 1 .- rev_dist
    df.best_fr_per_size_abs_impr = fr_col ./ fr_per_size_alltests .- 1
    return df
end

function dist_to_best_ex_per_size!(df::DataFrame, sdf::SubDataFrame)
    ex_col = sdf[:, :max_ex_per_size]
    ex_per_size_alltests = sdf[1, :max_exercises] /  sdf[1, :numtests]
    rev_dist = (ex_col .- ex_per_size_alltests) ./ (maximum(ex_col) - ex_per_size_alltests)
    df.dist_to_best_ex_per_size = 1 .- rev_dist
    df.best_ex_per_size_abs_impr = ex_col ./ ex_per_size_alltests .- 1
    return df
end

function dist_to_opt(sdf::SubDataFrame)::DataFrame
    df = DataFrame(sdf)
    dist_to_opt_size!(df, sdf)
    dist_to_best_ex_per_size!(df, sdf)
    dist_to_best_fr_per_size!(df, sdf)
    return df
end

function add_normalized_stats(df::DataFrame)
    a_per_inst = DataFrames.groupby(df, :instanceid)
    return reduce(append!, [ dist_to_opt(sdf) for sdf in a_per_inst ])
end
