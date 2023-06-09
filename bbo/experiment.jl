include("problem.jl")

using DataFrames

instanceidcounter = 0
next_instanceid() = global instanceidcounter += 1

struct TspInstance
    problem::TspProblem
    systemid::String
    instanceid::Int64
    scale::Real

    function TspInstance(original::System, _scale::Real,
                            objective_categories::Vector{Symbol}=Symbol[],
                            constraint_categories::Vector{Symbol}=Symbol[])
        actual = scale(original, _scale)
        case = generate_case_model(actual)
        instanceid = next_instanceid()
        return TspInstance(case, systemid(original), instanceid, _scale, objective_categories, constraint_categories)
    end

    function TspInstance(case::AbstractCase,
                            caseid::AbstractString,
                            objective_categories::Vector{Symbol}=Symbol[],
                            constraint_categories::Vector{Symbol}=Symbol[])
        systemid, iid, _scale = split(caseid, '_')

        scale_num = occursin(".", _scale) ? parse(Float64, _scale) : parse(Int64, _scale)
        return TspInstance(case, systemid, parse(Int64, iid), scale_num, objective_categories, constraint_categories)
    end

    function TspInstance(case::AbstractCase, systemid::AbstractString, iid::Int64, scale::Real,
                            objective_categories::Vector{Symbol}=Symbol[],
                            constraint_categories::Vector{Symbol}=Symbol[])
        problem = TspProblem(case, objective_categories, constraint_categories)
        TspInstance(problem, systemid, iid, scale)
    end

    function TspInstance(problem::TspProblem, systemid::AbstractString, iid::Int64=1, scale::Real=1.0)
        new(problem, systemid, iid, scale)
    end
end

problem(i::TspInstance) = i.problem
systemid(i::TspInstance) = i.systemid
instanceid(i::TspInstance) = i.instanceid
scale(i::TspInstance) = i.scale
id(i::TspInstance) = "$(systemid(i))_$(instanceid(i))_$(scale(i))"

@forward TspInstance.problem case, numtests

runidcounter = 0
next_runid() = global runidcounter += 1

struct TspExperiment
    instance::TspInstance
    algid::String
    runid::Int64

    function TspExperiment(instance, algid, runid = next_runid())
        new(instance, algid, runid)
    end
end

@forward TspExperiment.instance case, systemid, instanceid, scale, numtests

algid(e::TspExperiment) = e.algid
runid(e::TspExperiment) = e.runid
id(e::TspExperiment) = "$(id(e.instance))_$(algid(e))_$(runid(e))"

struct TspExperimentResult
    experiment::TspExperiment
    runtimes
    iterations::Integer
    evaluations::Integer
    pareto_front::AbstractVector{<:AbstractVector{<:Integer}}
    history_solutions::AbstractVector{<:AbstractVector{<:Integer}}

    function TspExperimentResult(exp, r, i, e, pf, hs = Vector{Vector{UInt8}}())
        return new(exp, r, i, e, pf, hs)
    end

    function TspExperimentResult(e::TspExperiment, r::BlackBoxOptim.OptimizationResults)
        nds = map(c -> c.inner.params, pareto_frontier(r)) # non dominated set
        _solving_times = solving_times(r)
        _iterations = BlackBoxOptim.iterations(r)
        _evaluations = BlackBoxOptim.f_calls(r)
        _history_solutions = history_solutions(r)

        return new(e, _solving_times,_iterations, _evaluations, nds, _history_solutions)
    end
end

@forward TspExperimentResult.experiment case, id, systemid, instanceid, scale, algid, runid, numtests

runtimes(r::TspExperimentResult) = r.runtimes
runtimesstr(r::TspExperimentResult) = string(r.runtimes)
pf(r::TspExperimentResult) = r.pareto_front
history_solutions(r::TspExperimentResult) = r.history_solutions
iterations(r::TspExperimentResult) = r.iterations
evaluations(r::TspExperimentResult) = r.evaluations

empty_stats(::CachedCase) = DataFrame(size = Int64[],
                              exercises = Float64[],
                              failrate = Float64[],
                              coverage = Float64[])

function stats(case::CachedCase, solutions::AbstractVector{<:AbstractVector{<:Integer}})::DataFrame
    _stats = empty_stats(case)

    for ind in solutions
        _size = tsp_size(ind)
        _exercises = tsp_exercises(ind, funcspertest(case))
        _failrate = tsp_failrate(ind, failrates(case))
        _coverage = tsp_coverage(ind, activationvectors(case), numfuncs(case))

        _fitness = [ _size, _exercises, _failrate, _coverage ]
        push!(_stats, _fitness)
    end

    for cat in case |> categories |> keys
        cat_df = categories(case)[cat]
        _stats[!, cat] = map(ind -> tsp_criterion_coverage(ind, cat_df), solutions)
    end

    return _stats
end

stats(r::TspExperimentResult) = stats(case(r), pf(r))
stats_hist(r::TspExperimentResult) = stats(case(r), history_solutions(r))

struct TspExperimentResultSummary
    systemid::String
    instanceid::Int64
    scale::Float64
    algid::String
    runid::Int64
    runtimes
    iterations::Integer
    evaluations::Integer
    stats::DataFrame
    categories::Vector{Symbol}
    numtests
    numfuncs
    avgfailrate
    avgexercises

    function TspExperimentResultSummary(r::TspExperimentResult)
        c = r |> case
        nt = c |> numtests
        nf = c |> numfuncs
        fr = c |> avgfailrate
        er = c |> avgexercises
        cts = c |> categories |> keys |> collect
        new(systemid(r), instanceid(r), scale(r),
            algid(r), runid(r), runtimes(r), iterations(r),
            evaluations(r), stats(r), cts, nt, nf, fr, er)
    end
end

systemid(s::TspExperimentResultSummary) = s.systemid
instanceid(s::TspExperimentResultSummary) = s.instanceid
scale(s::TspExperimentResultSummary) = s.scale
algid(s::TspExperimentResultSummary) = s.algid
runid(s::TspExperimentResultSummary) = s.runid
runtimes(s::TspExperimentResultSummary) = s.runtimes
iterations(s::TspExperimentResultSummary) = s.iterations
evaluations(s::TspExperimentResultSummary) = s.evaluations
stats(s::TspExperimentResultSummary) = s.stats
numtests(s::TspExperimentResultSummary) = s.numtests
numfuncs(s::TspExperimentResultSummary) = s.numfuncs
avgfailrate(s::TspExperimentResultSummary) = s.avgfailrate
maxfailrate(s::TspExperimentResultSummary) = s.maxfailrate
avgexercises(s::TspExperimentResultSummary) = s.avgexercises
maxexercises(s::TspExperimentResultSummary) = s.maxexercises
categories(s::TspExperimentResultSummary) = s.categories

struct TspExperimentalSeries
    instanceids::Vector{String}
    optimizer_setup::Vector{<:ParamsDict}
    runs::Int64
    exectime::Int64
    objective_categories::Vector{Symbol}
    constraint_categories::Vector{Symbol}
    expnr::Union{Integer, String}

    function TspExperimentalSeries(o::Vector{<:Dict},
                                r::Int64,
                                e::Int64,
                                oc::Vector{Symbol}=Symbol[],
                                cc::Vector{Symbol}=Symbol[],
                                exp_nr::Union{Integer, String} = current_experiment_id())
        i = case_ids_for(exp_nr)
        TspExperimentalSeries(i, o, r, e, oc, cc, exp_nr)
    end

    function TspExperimentalSeries(i::Vector{String},
                                o::Vector{<:Dict},
                                r::Int64,
                                e::Int64,
                                oc::Vector{Symbol}=Symbol[],
                                cc::Vector{Symbol}=Symbol[],
                                exp_nr::Union{Integer, String} = current_experiment_id())
        new(i, o, r, e, oc, cc, exp_nr)
    end
end

instaneceids(s::TspExperimentalSeries) = s.instanceids
optimizer_setup(s::TspExperimentalSeries) = s.optimizer_setup
runs(s::TspExperimentalSeries) = s.runs
exectime(s::TspExperimentalSeries) = s.exectime
expnr(s::TspExperimentalSeries) = s.expnr
objective_categories(s::TspExperimentalSeries) = s.objective_categories
constraint_categories(s::TspExperimentalSeries) = s.constraint_categories

function totalexperiments(s::TspExperimentalSeries)::Integer
    length(s.instanceids) * length(s.optimizer_setup) * s.runs
end

function totalexectime(s::TspExperimentalSeries)::String
    time =totalexperiments(s) * s.exectime / 60 / 60
    return "$(round(time, sigdigits = 2)) hours"
end

function experiment_already_executed(runid::Integer, expid::Union{Integer, AbstractString})
    result_file = "experiment_$expid/results/results.csv"
    if !isfile(result_file)
        return false
    end

    wc_text = read(`wc -l $result_file`, String)
    runssofar = parse(Int64, split(wc_text, " ")[1]) - 1 # first row is header
    return runid ≤ runssofar
end

function run(s::TspExperimentalSeries)
    expcounter = 0
    global runidcounter = 0 # have experiment runs start at count 1
    for iid in instaneceids(s)
        instance = load_instance(iid, objective_categories(s), constraint_categories(s), expnr(s))
        for (idx, optimizer_setup) in enumerate(optimizer_setup(s))
            opt_name = optimizer_setup[:AlgorithmName]
            optimizer_setup[:MaxTime] = exectime(s)
            for run in 1:runs(s)
                expcounter += 1
                if experiment_already_executed(expcounter, expnr(s))
                    global runidcounter += 1
                    continue
                end

                println("===============optimizer $opt_name on instance $iid (exp $expcounter/$(totalexperiments(s)), total time = $(totalexectime(s)))========================")
                experiment = TspExperiment(instance, opt_name)
                res = bboptimize(problem(instance); optimizer_setup...)
                if expcounter == 1 # very first experiment, do it twice to not get unfair comparison for jit compiler
                    res = bboptimize(problem(instance); optimizer_setup...)
                end
                tsp_results = TspExperimentResult(experiment, res)
                save(tsp_results, expnr(s))
            end
        end
    end
end
