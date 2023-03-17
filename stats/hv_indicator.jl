using DataFrames, LinearAlgebra, Serialization

include("../normalized_stats.jl")
include("app_hvi.jl")

HVI_GRANULARITY = 60

hvi(m, ref) = approximate_hypervolume_ms(m, zeros(size(m, 1)), ref)

# since these are all historical entries into pareto front, make sure to reduce to pareto entries at each stage.
# objectives must be max aligned
# each column in df is an objective with its readings in order
function hv_indicators(df::DataFrame, min_refs=nothing, max_refs=nothing)::Vector{<:Real}
    "start!" |> println

    prep = Matrix{Float64}(df)
    prep = Matrix(prep')

    min_ref = isnothing(min_refs) ? minimum(prep, dims=2)[:,1] : copy(min_refs)
    max_ref = isnothing(max_refs) ? maximum(prep, dims=2)[:,1] : copy(max_refs)

    dims = size(prep, 1)
    # normalize datapoints
    for dim in 1:dims
        # ensure no divide by zero (if min and max equal)
        if max_ref[dim] == min_ref[dim]
            prep[dim,:] .= 1.0
        else
            prep[dim,:] = (prep[dim,:] .- min_ref[dim]) ./ ( max_ref[dim] - min_ref[dim])
        end
    end

    # set default to 0 (1), for better overview, and since everything below 0 (above 1) gives the same rating, 0
    prep[prep.<=0] .= 0
    prep[prep.>=1] .= 1

    prep = 1 .- prep

    ref = max_refs |> length |> ones

    num_evals = size(prep,2)

    discretizer = num_evals > HVI_GRANULARITY ? div(num_evals, HVI_GRANULARITY) : 1
    "discretizer $discretizer" |> println

    contribs = Vector{Float64}(undef, num_evals)
    current = 0
    "evals $num_evals" |> println
    _mem_part_prep = Array{Float64}(undef, size(prep, 1), 0)
    former_index = 1
    for i in 1:num_evals
        if i % discretizer != 0
            contribs[i] = current
            continue
        end

        partprep = prep[:, former_index:i]
        "index: $i" |> println

        if false
            _mem_part_prep = reduce_to_non_dominating(partprep)
        else
            former_index = i+1
            _mem_part_prep = reduce_to_non_dominating(_mem_part_prep, partprep)
        end

        contrib = hvi(_mem_part_prep, ref)

        # since the hvi is approximative, but hvi is monotonic, correct for smaller values
        if current < contrib
            current = contrib
        end
        contribs[i] = current
    end

    if discretizer != 1
        _mem_part_prep = reduce_to_non_dominating(prep)
        contrib = hvi(_mem_part_prep, ref) # make sure to register final result
        if current < contrib
            current = contrib
        end
        contribs[end] =  current
    end

    return contribs
end

function hypervolumes(expid::Union{String, Integer},
                _caseid::String,
                _objectives,
                _hv_objectives_tomax=fill(identity, length(objectives)),
                worst_refs=nothing,
                best_refs=nothing
                )::Dict{Integer, Vector{Float64}}
    _case = load_case(_caseid, expid)
    case_id_filter = (f) -> startswith(f, _caseid)

    set_refs = (refs, obj_max, objs) ->
        begin
            refs = isnothing(refs) ? nothing : copy(refs)
            if !isnothing(refs)
                refs = [ obj_max[i](refs[i]) for i in eachindex(objs)]
            end
            return refs
        end

    max_refs = set_refs(best_refs, _hv_objectives_tomax, _objectives)
    min_refs = set_refs(worst_refs, _hv_objectives_tomax, _objectives)

    hvis = Dict{Integer, Vector{Float64}}()

    count = 0
    _runids = runids(expid, case_id_filter)
    for runid in _runids
        count += 1
        _hvi_file = joinpath("experiment_$expid", "results", "stats", "$(runid)_hvi.ser")
        if isfile(_hvi_file) # && runid != 1
            hvis[runid] = deserialize(_hvi_file)
            continue
        end

        "hypervolume import for $runid ... ($count/$(length(_runids))) " |> print
        "hist stats... " |> print
        _stats = load_hist_stats(expid, runid)
        _performance = _stats[:, _objectives]
        for obj in 1:length(_objectives)
            _performance[:, obj] = _hv_objectives_tomax[obj].(_performance[:, obj])
        end

        "indicators... " |> print

        hvis[runid] = hv_indicators(_performance, min_refs, max_refs)
        serialize(_hvi_file, hvis[runid])
        "done!" |> println
    end

    return hvis
end

function accumulated_runtimes(_runtimes::Vector{Float64})
    _acc_runtimes = Vector{Float64}(undef, length(_runtimes))
    _tsf = 0
    for run in eachindex(_runtimes)
        val = _runtimes[run]
        _tsf += val
        _acc_runtimes[run] += _tsf
    end
    return _acc_runtimes
end

# in the heuristic case, the numbers are already accumulated
function hypervolume_per_run(_hv::Vector{Float64}, _results::DataFrame, runid::Int64, _is_heuristic::Bool=false)::Tuple
    _runtimes = filter(r -> r.runid == runid, _results).runtimes[1]
    _runtimes = (eval ∘ Meta.parse)(_runtimes)
    _acc_runtimes = _is_heuristic ? _runtimes : accumulated_runtimes(_runtimes)

    # normalize readings to start by 0 in graph
    num_zeros = 2 + (length(_acc_runtimes) - length(_hv))
    _hv = vcat(zeros(num_zeros), _hv)
    _acc_runtimes = vcat([0, _acc_runtimes[1]], _acc_runtimes)

    return _hv, _acc_runtimes
end

function discretize_hvis(_acc_runtimes::Vector{Float64}, _hvis::Vector{Float64}, steps::Integer=20, boundary::Real = maximum(_runtimes))
    @assert length(_acc_runtimes) == length(_hvis)

    _times = Vector{Float64}(undef, steps)
    _values = Vector{Float64}(undef, steps)

    _stepsize = boundary / steps

    _current = _stepsize
    for step in 1:steps
        _times[step] = _current

        # set correct value
        _index = 1
        while length(_hvis) > _index && _acc_runtimes[_index] < _current
            _index += 1
        end
        _values[step] = _hvis[_index]

        _current += _stepsize
    end

    if mean(_values) == 0
        _values |> println
    end

    #@assert mean(_values) > 0
    return _times, _values
end

using Plots

struct Hypervolume
    progress
    times
    original_progress
    original_acc_times
end

function hypervolumes(expid::Union{String, Integer},
                        _caseid::String,
                        _hv_objectives,
                        _total_time_s::Int64,
                        _hv_objectives_tomax=fill(identity, length(_hv_objectives)),
                        _worst_refs=nothing,
                        _best_refs=nothing,
                        _discr_steps::Int64=60,
                        _results::DataFrame=load_tsp_results(expid))::Dict{String, Dict{Int64, Hypervolume}}

    _hvs = hypervolumes(expid, _caseid, _hv_objectives, _hv_objectives_tomax, _worst_refs, _best_refs)

    _hypervolumes = Dict()
    for _algid in unique(_results.algid)
        _alg_runids = filter(r -> r.algid == _algid && r.runid ∈ keys(_hvs), _results).runid

        _hvs_per_run = Dict{Int64, Hypervolume}()
        for runid in _alg_runids
            runid |> println
            _is_heuristic = _algid ∈ HEURISTIC_ALG_NAMES
            _hv, _acc_runtimes = hypervolume_per_run(_hvs[runid], _results, runid, _is_heuristic)
            _progress, _times = discretize_hvis(_acc_runtimes, _hv, _discr_steps, _total_time_s)
            _hvs_per_run[runid] = Hypervolume(_times, _progress, _hv, _acc_runtimes)
        end

        _hypervolumes[_algid] = _hvs_per_run
    end

    return _hypervolumes
end

max_hvi_val(_hvis, _algids = collect(keys(_hvis))) = map(a -> map(v -> v.progress[end], values(_hvis[a])) |> maximum, _algids) |> maximum

function plot_hypervolumes(_hvis::Dict{String, Dict{Int64, Hypervolume}})
    foreach(k -> isempty(_hvis[k]) ? delete!(_hvis, k) : nothing, keys(_hvis)) # clean from empty garbage
    _algids = _hvis |> keys |> collect
    _algids = sort(_algids, by=x->_ALG_NAMES[x]) # ensure best possible order here
    _first_alg = _hvis[ _algids |> first ]
    _some_hv = _first_alg |> values |> collect |> first

    _combined_progress = foldl(hcat, map(a -> a.progress, values(_first_alg)))

    _mean_progress = mean(_combined_progress, dims=2)
    _std_progress = std(_combined_progress, dims=2)

    _max_val = max_hvi_val(_hvis, _algids)

    _p = plot(_some_hv.times, _mean_progress, ribbon = _std_progress,
                                        ylims = [0,_max_val*1.1],
                                        labels = _algids |> first,
                                        legend = :topleft,
                                        _cc = colorindex_for(_algids |> first),
                                        ylabel="hypervolume indicator",
                                        xlabel="time (s)",
                                        fillalpha = 0.2,
                                        xlims=[0,_some_hv.times[end]])

    lss = [:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
    for _algpos in 2:length(_algids)
        _algid = _algids[_algpos]
        _cc = colorindex_for(_algid)
        _combined_progress = foldl(hcat, map(a -> a.progress, values(_hvis[_algid])))

        _mean_progress = mean(_combined_progress, dims=2)
        _std_progress = std(_combined_progress, dims=2)
        _p = plot!(_some_hv.times,
                    _mean_progress,
                    ribbon = _std_progress,
                    color = _cc,
                    linestyle = lss[_algpos % (length(lss)-1)+1],
                    fillalpha = 0.2,
                    labels = _algid)
    end

    return _p
end

default_caseid(expid) = caseid(expid, 1) # assuming there is always a run 1

function coverage_criteria_indices(expid, _objective_labels, caseid=default_caseid(expid))
    _criteria_labels = coverage_criteria_labels(expid, caseid)
    return findall(o -> o ∈ _criteria_labels, _objective_labels)
end

function extreme_values(expid, _labels, obj_transf, expfilter=(x)->true)

    _runids = runids(expid, expfilter)

    local _worst
    local _best
    "retrieving extreme values..." |> println
    count = 1
    for runid in _runids
        "runid $(runid) : $count / $(length(_runids))" |> println
        hist_stats = load_hist_stats(expid, runid)

        current_worst = map(x -> minimum(obj_transf[x[1]](hist_stats[!, x[2]])), enumerate(_labels))
        current_best = map(x -> maximum(obj_transf[x[1]](hist_stats[!, x[2]])), enumerate(_labels))

        if @isdefined _worst
            _worst = map(x -> min(x[2], current_worst[x[1]]), enumerate(_worst))
            _best = map(x -> max(x[2], current_best[x[1]]), enumerate(_best))
        else
            _worst = current_worst
            _best = current_best
        end
        count += 1
    end

    "done!" |> println

    _worst = map(i -> obj_transf[i](_worst[i]), eachindex(_worst))
    _best = map(i -> obj_transf[i](_best[i]), eachindex(_best))

    _indices = coverage_criteria_indices(expid, _labels)
    for i in _indices
        _worst[i] = 0.0
        _best[i] = 1.0
    end

    return _worst, _best
end

# must be in same order as on plot!
MC_TOA_ALGS = [ "gurobi", "Cbc", "Clp", "Cbc_relaxed", "gurobi_relaxed"]

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

function plotphasemarkeron(_p, _run_hvis, _pos, _phase, _colorindex, _runtime, δ)
    shortest =  map(h -> length(h.original_acc_times), _run_hvis) |> minimum
    if _pos > shortest
        return
    end

    x = map(h -> h.original_acc_times[_pos], _run_hvis) |> mean
    y = map(h -> h.original_progress[_pos], _run_hvis) |> mean

    if x ≤ _runtime
        plot!(_p, rectangle(8,2δ,x-4,y-δ),
            linecolor = _colorindex,
            color = _colorindex,
            label="")
        annotate!( [ (x, y, text("$_phase",8,:white,:center)) ] )
    end
end

function add_phase_markers(_hvis, _p, _numobj, _runtime)

    _max_val = max_hvi_val(_hvis)
    δ = .025*_max_val

    # Dict{String, Dict{Int64, Hypervolume}}
    for _algid in keys(_hvis)
        if (_algid ∈ MC_TOA_ALGS)
            _run_hvis = values(_hvis[_algid])
            _cc = colorindex_for(_algid)
            # third index, the first two are setting the initial trajectory
            plotphasemarkeron(_p, _run_hvis, 3, 2, _cc, _runtime, δ)
            plotphasemarkeron(_p, _run_hvis, 3+_numobj, 3, _cc, _runtime, δ)
        end
    end
end

function theoretical_extremes(expid, _caseid, _objectives_labels)
    _case = load_case(_caseid, expid)

    _best = Vector{Float64}(undef, length(_objectives_labels))
    _worst = copy(_best)

    local _boundary
    for (i, l) in enumerate(_objectives_labels)
        if l == "size"
            _boundary = [numtests(_case), 1]
        elseif l == "failrate"
            _boundary = [0,maxfailrate(_case)]
        elseif l == "exercises"
            _boundary = [1,maxexercises(_case)]
        elseif l ∈ ["features", "morphologies", "coverage"]
            _boundary = [0,1]
        end

        _worst[i], _best[i] = _boundary
    end

    return _worst, _best
end

function plot_hypervolumes(expid, _objectives_labels, _objectives_tomax, _runtime, _caseid, _results=load_tsp_results(expid), expfilter = (f) -> startswith(f, _caseid))
    _worst_values, _best_values = extreme_values(expid, _objectives_labels, _objectives_tomax, expfilter)

    _hvis = hypervolumes(expid, _caseid, _objectives_labels, _runtime, _objectives_tomax, _worst_values, _best_values, 60, _results)
    _p = plot_hypervolumes(_hvis)

    _numobj = length(_objectives_labels)
    add_phase_markers(_hvis, _p, _numobj, _runtime)
    return _p
end
