include("bbo/experiment.jl")

using DelimitedFiles, DataFrames, CSV, CodecZlib, BenchmarkTools

const RESULTS_FILE = "results.csv"

const ACTIVATION_MATRIX_COMP = "activationmatrix.csv.gz"

const HIST_ENDING_COMP = "_hist.gz"

const HIST_ENDING_CSV = "_hist.csv"

const PF_ENDING_COMP = "_pf.gz"

const PF_ENDING_CSV = "_pf.csv"

const CAT_CSV_ENDING = "_cat.csv"

const FAIL_RATES = "failrates.vector"

# OBS minsize here assumes 100% coverage, i.e. all functions are exercised
function empty_results(r::TspExperimentResult)::DataFrame
    return DataFrame(algid = String[],
              systemid = String[],
              instanceid = Int64[],
              scale = Real[], # both scales in actual size and "scale" supported
              runid = Int64[],
              iterations = Int64[],
              evaluations = Int64[],
              runtimes = String[])
end

# OBS size is only right as we guarantee entires in both final column and row (each test covers at least one function and each function is covered)
save_binary_to(am::Matrix{T}, file_name::String) where {T <: Integer} = save_binary_to(sparse(am), file_name)
function save_binary_to(am::SparseMatrixCSC{T}, file_name::String) where {T <: Integer}
    I, J, _ = findnz(am)
    df = DataFrame([:I => I, :J => J])

    open(GzipCompressorStream, file_name, "w") do f
        CSV.write(f, df)
    end
end

function load_binary_matrix_from(file_name::String, datatype::Type=UInt8)
    df = open(GzipDecompressorStream, file_name, "r") do f
        DataFrame(CSV.File(f))
    end
    return sparse(df[:, :I], df[:, :J], ones(datatype, nrow(df)))
end

function save_to(what, file_name::String)
    open(file_name, "w") do f
        writedlm(f, what)
    end
end

function save_binary_vectors_to(vs::AbstractVector{<:AbstractVector{<:Integer}}, file_name::String)
    open(GzipCompressorStream, file_name, "w") do f
        writedlm(f, vs)
    end
end

function load_binary_vectors_from(file_name::String, datatype::Type{T}=UInt8)::Vector{SparseVector{T}} where {T <: Integer}
    pf_matrix = open(GzipDecompressorStream, file_name, "r") do f
        readdlm(f, datatype)
    end

    return [ sparsevec(pf_matrix[i,:]) for i in 1:size(pf_matrix,1) ]
end

function load_from(file_name::String, datatype::Type=UInt8)
    return open(file_name, "r") do f
        readdlm(f, datatype)
    end
end

function load_categories(dir::AbstractString)
    dfs = Dict{Symbol, DataFrame}()
    _categories = filter(f -> endswith(f, CAT_CSV_ENDING), readdir(dir))
    foreach(c -> dfs[Symbol(split(c, "_")[1])] = DataFrame(CSV.File(joinpath(dir,c))), _categories)
    return dfs
end

function load_case(dir::AbstractString, cached::Bool = false)::AbstractCase
    am = load_binary_matrix_from(joinpath(dir, ACTIVATION_MATRIX_COMP))
    fr = load_from(joinpath(dir, FAIL_RATES), Float64)[:,1]
    _categories = load_categories(dir)

    case = Case(am, _categories, fr)

    return cached ? CachedCase(case) : case
end

function save(c::AbstractCase, dir::String)
    mkdir(dir)
    save_binary_to(activationmatrix(c), joinpath(dir, ACTIVATION_MATRIX_COMP))
    save_to(failrates(c), joinpath(dir, FAIL_RATES))
    foreach(p -> CSV.write(joinpath(dir, "$(p[1])$CAT_CSV_ENDING"), p[2]), pairs(categories(c)))
end

save(i::TspInstance, expdir::String) = save(case(i), joinpath(expdir, id(i)))

const RESULT_HANDLES = [ algid systemid instanceid numtests runid iterations evaluations runtimesstr ]
const RESULT_TYPES = Type[ String String Int64 Real Int64 Int64 Int64 String ]

function initialized_frame(r::TspExperimentResult)
    entry = [ h(r) for h in RESULT_HANDLES ]
    return push!(empty_results(r), entry)
end

save(r::TspExperimentResult, expid::Integer=1) = save(r, "experiment_$expid")

function save(r::TspExperimentResult, expdir::AbstractString)
    expdir = startswith(expdir, "experiment_") ? expdir : "experiment_$expdir"
    resdir = joinpath(expdir, "results")
    if !isdir(resdir)
        mkdir(resdir)
    end

    save_binary_vectors_to(pf(r), "$resdir/$(id(r))$(PF_ENDING_COMP)")

    res_file = joinpath(resdir, RESULTS_FILE)
    df = initialized_frame(r)

    if ! (r |> history_solutions |> isempty)
        save_binary_vectors_to(history_solutions(r), "$resdir/$(id(r))$(HIST_ENDING_COMP)")
    elseif startswith(algid(r), "borg")
        save_binary_vectors_to(cached_history_solutions, "$resdir/$(id(r))$(HIST_ENDING_COMP)")
        df.runtimes[1] = cached_solving_times |> (x -> round.(x, digits=2, RoundUp)) |> string
    end
    reset_heuristic_stats()

    ((f) -> isfile(f) ? CSV.write(f, df, append=true) : CSV.write(f, df))(res_file)
end

function load_tsp_stats(expid::AbstractString,
                        objective_categories::Vector{Symbol}=Symbol[],
                        constraint_categories::Vector{Symbol}=Symbol[],
                        exp_nr::Union{String, Integer}=1,
                        c::CachedCase=case_for(expid, exp_nr),
                        results::DataFrame=load_tsp_results(exp_nr))

    sysid, iid, sc, algid, runid = split(expid, "_")

    sc = parse(Float64, sc)
    iid, runid = parse.(Int64, [ iid, runid ])

    instance = TspInstance(c, sysid, iid, sc, objective_categories, constraint_categories)
    exp = TspExperiment(instance, algid, runid)
    _runtimes, _iterations, _evaluations = filter(r -> r.runid == runid, results)[1,[:runtimes, :iterations, :evaluations]]

    pf = load_binary_vectors_from("experiment_$expid/results/$(expid)$(PF_ENDING_COMP)")

    histfile = "experiment_$expid/results/$(expid)$(HIST_ENDING_COMP)"

    local res
    if isfile(histfile)
        hs = load_binary_vectors_from(histfile)
        res = TspExperimentResult(exp, eval(Meta.parse(_runtimes)), _iterations, _evaluations, pf, hs)
    else
        res = TspExperimentResult(exp, eval(Meta.parse(_runtimes)), _iterations, _evaluations, pf)
    end

    return res
end

function current_experiment_id()
    counter = 1
    dir = "experiment_$(counter)"

    while isdir(dir)
        counter += 1
        dir = "experiment_$(counter)"
    end

    return counter-1
end

current_experiment_dir() = "experiment_$(current_experiment_id())"
next_experiment_dir() = mkdir("experiment_$(current_experiment_id()+1)")
function experiment_dir(name::String)
    full_name = "experiment_$name"
    if isdir(full_name)
        ArgumentError("\"$full_name\" directory exists already") |> throw
    else
        mkdir(full_name)
    end
    return full_name
end

using Plots

using Plots: Plot

function allpairs(n::Integer)::Vector{Vector{Int64}}

    combs = Vector{Vector{Int64}}()

    for i in 1:(n-1)
        for j in (i+1):n
            push!(combs, Integer[i, j])
        end
    end

    return combs
end

# OBS results must have same objectives; can be used for single one
function plot_pareto_fronts(rs::TspExperimentResult...)
    pf_stats = stats.(rs)
    _algids = map(r -> algid(r), rs)
    return plot_pareto_fronts(_algids, pf_stats)
end

function plot_pareto_fronts(_algids::Vector{String}, pf_stats::DataFrame...)
    dims = ncol(pf_stats[1])
    pairs = allpairs(dims)

    plots = Vector{Plot}(undef, length(pairs))

    "plotting..." |> print
    _shapes = [:rect, :star4, :rect, :star4, :diamond, :cross, :xcross]
    _alphas = [1,1,.6,.6,.6,.2,.2]

    local labels
    for (idx, pair) in enumerate(pairs)
        x = pf_stats[1][:,pair[1]]
        y = pf_stats[1][:,pair[2]]
        labels = names(pf_stats[1])[pair]
        _label = idx == 1 ? _algids[1] : :none
        plots[idx] = plot(x, y,
                          markershape = _shapes[1],
                          markersize = 5,
                          markerstrokewidth = 0,
                          markeralpha = _alphas[1],
                          seriestype = :scatter,
                          color = colorindex_for(_algids[1]),
                          xlabel = labels[1],
                          ylabel = labels[2],
                          xtickfont=font(8),
                          ytickfont=font(8),
                          xguidefontsize=12,
                          yguidefontsize=12,
                          legendfontsize=9,
                          legend = :topleft,
                          label = _label)

         for r_idx in 2:length(pf_stats)
             x = pf_stats[r_idx][:,pair[1]]
             y = pf_stats[r_idx][:,pair[2]]
             _label = idx == 1 ? _algids[r_idx] : :none
            plot!(x, y,
                markershape = _shapes[r_idx],
                markersize = 5,
                markerstrokewidth = 0,
                markeralpha = _alphas[r_idx],
                seriestype = :scatter,
                color = colorindex_for(_algids[r_idx]),
                label = _label)
         end
    end
    "done!" |> println

    return plots
end

function load_tsp_results(exp_nr::Union{String, Integer}=1)::DataFrame
    res = DataFrame!(CSV.File("experiment_$(exp_nr)/results/results.csv"))
    return sort(res, :runid)
end

function caseid(expid::AbstractString)
    sysid, iid, sc, algid, runid = split(expid, "_")
    return "$(sysid)_$(iid)_$(sc)"
end

function load_case(caseid::AbstractString, exp_nr::Union{String, Integer}=1)
    load_case("experiment_$exp_nr/$caseid", true)
end

function case_ids_for(expr_nr::Union{String, Integer}=1)
    pdir = "experiment_$expr_nr"
    ids = filter(d -> isdir(joinpath(pdir, d)) && basename(d) != "results", readdir(pdir))
    return sort(ids, by= id -> parse(Int64, split(id, '_')[3]))
end

function load_instances(objective_categories::Vector{Symbol}=Symbol[],
                        constraint_categories::Vector{Symbol}=Symbol[],
                        expr_nr::Union{String, Integer}=1)
    case_ids = case_ids_for(expr_nr)
    _cases = map(c -> load_case(c, expr_nr), case_ids)

    return map(i -> TspInstance(_cases[i], case_ids[i], objective_categories, constraint_categories), eachindex(_cases))
end

function load_instance(iid::AbstractString,
                            objective_categories::Vector{Symbol}=Symbol[],
                            constraint_categories::Vector{Symbol}=Symbol[],
                            expr_nr::Union{Integer, String}=1)
    TspInstance(load_case(iid, expr_nr), iid, objective_categories, constraint_categories)
end

function load_experimental_results(runid::Integer, exp_nr::Union{String, Integer}=1, results=load_tsp_results(exp_nr))
    println("runid: $runid")
    pdir = "experiment_$exp_nr/results/"
    pf_file = filter(f -> endswith(f, "_$runid$PF_ENDING_COMP"), readdir(pdir))[1]
    hs_file = filter(f -> endswith(f, "_$runid$HIST_ENDING_COMP"), readdir(pdir))[1]
    case_id = join(split(pf_file, "_")[1:3], "_")
    c = load_case(case_id, exp_nr)
    i = TspInstance(c, case_id)

    algid = join(split(pf_file, "_")[4:end-2], "_")
    exp = TspExperiment(i, algid, runid)

    run_row = filter(r -> r.runid == runid, eachrow(results))[1]
    pf = load_binary_vectors_from(joinpath(pdir, pf_file))
    hs = load_binary_vectors_from(joinpath(pdir, hs_file))

    return TspExperimentResult(exp,
                               eval(Meta.parse.(run_row[:runtimes])),
                               run_row[:iterations],
                               run_row[:evaluations],
                               pf, hs)
end

function load_hist(expid::Union{String, Integer}, runid::Integer)::Vector{SparseVector{UInt8}}
    res_dir = joinpath("experiment_$(expid)","results")
    hists = filter(f -> endswith(f, "_$runid$HIST_ENDING_COMP") , readdir(res_dir))
    @assert length(hists) == 1
    hist_path = joinpath(res_dir, hists[1])
    return load_binary_vectors_from(hist_path)
end

function load_pf(expid::Union{String, Integer}, runid::Integer)::Vector{SparseVector{UInt8}}
    res_dir = joinpath("experiment_$(expid)","results")
    pfs = filter(f -> endswith(f, "_$runid$PF_ENDING_COMP") , readdir(res_dir))
    @assert length(pfs) == 1
    pf_path = joinpath(res_dir, pfs[1])
    return load_binary_vectors_from(pf_path)
end

function caseid(expid::Union{String, Integer}, runid::Integer)::AbstractString
    res_dir = joinpath("experiment_$(expid)","results")
    hists = filter(f -> endswith(f, "_$runid$HIST_ENDING_COMP") , readdir(res_dir))
    @assert length(hists) == 1
    parts = split.(hists[1], "_")[1:3]
    return join(parts, "_")
end

function runids(expid::Union{String, Integer}, custom_filter::Function=(f)->true)::Vector{Int64}
    res_dir = joinpath("experiment_$(expid)","results")
    hists = filter(f -> endswith(f, HIST_ENDING_COMP) , readdir(res_dir))
    hists = filter(custom_filter , hists)
    parse.(Int64, map(f -> f[end-1], split.(hists, "_")))
end

function load_hist_stats(expid::Union{String, Integer},
                        runid::Integer)
    _stats_file = joinpath("experiment_$expid", "results", "stats", "$runid$HIST_ENDING_CSV")
    if isfile(_stats_file)
        return CSV.read(_stats_file, DataFrame)
    end

    _caseid = caseid(expid, runid)
    case = load_case(_caseid, expid)
    hist = load_hist(expid, runid)
    return stats(case, hist)
end

function load_pf_stats(expid::Union{String, Integer},
                        runid::Integer)
    _caseid = caseid(expid, runid)
    case = load_case(_caseid, expid)
    pf = load_pf(expid, runid)
    return stats(case, pf)
end

function load_hists(expid::Union{String, Integer})::Dict{Int64, Vector{SparseVector{UInt8}}}
    _runids = runids(expid)
    hist_bins = map(runid -> load_hist(expid, runid), _runids)
    Dict(_runids .=> hist_bins)
end

function load_experimental_summaries(i::TspInstance, exp_nr::Union{String, Integer}=1, results=load_tsp_results(expr_nr))
    println("instance: $(instanceid(i))")
    pdir = "experiment_$exp_nr/results/"
    all_results = readdir(pdir)

    pf_files = filter(f -> isfile(pdir * f) && startswith(f, id(i)) && endswith(f, PF_ENDING_COMP), all_results)

    "load pfs..." |> print
    pfs = map(load_binary_vectors_from, pdir .* pf_files)
    "done!" |> println

    algids = map(f ->
                    begin
                        entries = split(f, "_")
                        algstart = 4 + length(entries) - 6
                        return join(entries[4:algstart], "_")
                    end
        , pf_files)

    runids = map(f ->
                    begin
                        entries = split(f, "_")
                        runidstart = 5 + length(entries) - 6
                        return parse(Int64,entries[runidstart])
                    end
    , pf_files)

    df = DataFrame(algid = algids, runid = runids, pfs = pfs)
    df = sort(df, :runid)

    "extract stats..." |> print
    filtered_results = DataFrame(filter(r -> r.runid ∈ df[:, :runid], eachrow(results)))
    if filtered_results |> isempty
        return nothing
    end

    filtered_results = sort(filtered_results, :runid) # sort for correct assignment
    _runtimes = eval.(Meta.parse.(filtered_results[:, :runtimes]))
    _iterations = filtered_results[:,:iterations]
    _evaluations = filtered_results[:,:evaluations]
    "done!" |> println

    "create experiments..." |> print
    df = DataFrame(filter(e -> e.runid ∈ filtered_results.runid, eachrow(df)))
    df = sort(df, :runid)
    experiments = [ TspExperiment(i, df[j, :algid], df[j, :runid]) for j in 1:size(df, 1)]
    "done!" |> println

    "create results..." |> print
    exp_results = [ TspExperimentResult(experiments[j], _runtimes[j], _iterations[j], _evaluations[j], df[j, :pfs]) for j in 1:nrow(filtered_results)]
    "done!" |> println

    "extract summaries..." |> println

    toti = length(exp_results)
    _ts = map(r -> begin
                        println("$(r[1]) / $toti")
                        return TspExperimentResultSummary(r[2])
                    end
                    , enumerate(exp_results))

    "done!" |> println
    return _ts
end

function load_experimental_summary(exp_nr::Union{String, Integer}=1, results = load_tsp_results(exp_nr))
    case_ids = case_ids_for(exp_nr)

    local summaries
    for _case_id in case_ids
        _case =  load_case(_case_id, exp_nr)
        _instance = TspInstance(_case, _case_id)
        _summary = load_experimental_summaries(_instance, exp_nr, results)

        if isnothing(_summary)
            continue
        end

        if @isdefined summaries
            summaries = vcat(summaries, _summary)
        else
            summaries = _summary
        end
    end

    return sort!(summaries, by=runid)
end

function setup_new_experiment_with(experiment_dir::String, _systems::AbstractVector{System}, _scales::AbstractVector{<:Real}, variants::Integer)
    for system in _systems
        for _scale in _scales
            for variant in 1:variants
                instance = TspInstance(system, _scale)
                save(instance, experiment_dir)
            end
        end
    end
end

function setup_new_experiment_from(exp_nr::Union{Integer, String}, instance_ids::String...=case_ids_for(exp_nr)...)
    instance_handles = joinpath.("experiment_$exp_nr", instance_ids)
    ndir = next_experiment_dir()
    foreach(i -> cp(instance_handles[i], joinpath(ndir, instance_ids[i])), eachindex(instance_handles))
end

function save_category_for(exp_dir::String, caseid::String, criterion::DataFrame, criterion_name::String)
    file_name = joinpath(exp_dir, caseid, "$(criterion_name)$CAT_CSV_ENDING")
    println(file_name)
     CSV.write(file_name, criterion)
end

function load_result_stats(expid::Union{String, Integer})
    return DataFrame!(CSV.File("experiment_$(expid)/results/stats/stats.csv"))
end

function plot_2ds(expnr::Union{String, Integer}, filename::String, titlename::String, _filter, runids::Integer...)
    title = plot(title = titlename,
            grid = false,
            showaxis = false,
            ticks = false)

    _p = plot_2ds(expnr, runids...)[_filter]
    plots = _p

    plots = vcat(_p, plot(framestyle = :none))
    l = length(plots)
    _nrows = div(l, 3) + (l % 3 > 0 ? 1 : 0)

    if l % 3 > 0
        plots = vcat(plots, fill(plot(framestyle = :none), 3 - (l % 3)))
    end

    _p = plot(plots..., layout = (_nrows, 3), size=(1500,1500))
    #_p = plot(plots..., layout = (2, 2))
    _p = plot(title, _p, layout = @layout([A{0.01h}; B]))

    statsdir = joinpath("experiment_$expnr","results", "stats")
    mkpath(statsdir)

    savefig(joinpath(statsdir, filename))
    return _p
end

function algids(expnr, _runids, _results = CSV.read("experiment_$expnr/results/results.csv", DataFrame))
    _relevant = filter(r -> r[:runid] ∈ _runids, _results)
    return collect(_relevant[!, :algid]), _relevant[!, :runid]
end

function plot_2ds(expnr::Union{String, Integer}, runids::Integer...)
    _results = CSV.read("experiment_$expnr/results/results.csv", DataFrame)
    _algids, runids = algids(expnr, runids, _results)

    runids = sort(runids, by=x->_ALG_NAMES[_algids[findfirst(isequal(x), runids)]])
    _algids = sort(_algids, by=x->_ALG_NAMES[x])

    _stats = map(r -> load_pf_stats(expnr, r), runids)
    return plot_pareto_fronts(_algids, _stats...)
end

function extract_stats(expid::Union{String, Integer}, _caseid::String)
    _case = load_case(_caseid, expid)

    extract_stats_hist(expid, _case)
    extract_stats_pf(expid, _case)
end

function extract_stats(expid::Union{String, Integer}, _case::AbstractCase, stats_load_func::Function, file_ext::String)
    "extracting stats..." |> println
    _c = 0
    _runids = runids(expid)
    for runid in _runids
        _c += 1
        "$_c/$(length(_runids))" |> println
        _file = joinpath("experiment_$expid", "results", "stats", "$runid$file_ext")
        if isfile(_file)
            continue
        end
        _stats = stats_load_func(expid, runid)
        CSV.write(_file, _stats)
    end
    "done!" |> println
end

function extract_stats_hist(expid::Union{String, Integer}, _case::AbstractCase)
    extract_stats(expid, _case, load_hist_stats, HIST_ENDING_CSV)
end

function extract_stats_pf(expid::Union{String, Integer}, _case::AbstractCase)
    extract_stats(expid, _case, load_pf_stats, PF_ENDING_CSV)
end


function coverage_criteria_labels(expid, caseid)
    dir = joinpath("experiment_$expid", caseid)
    _criteria = filter(f -> endswith(f, CAT_CSV_ENDING), readdir(dir))
    _criteria = map(c -> split(c, "_")[1], _criteria)
    return push!(_criteria, "coverage")
end
