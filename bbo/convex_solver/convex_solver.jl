include("solver.jl")
include("../history_optimizer.jl")

using LinearAlgebra, Lazy

@enum SolverStage phase1 phase2 phase3

# contains an archive with pareto set
mutable struct TspSolver <: HistoryOptimizer
    stage::SolverStage
    problem::TspProblem
    archive
    evaluator::BlackBoxOptim.Evaluator
    options
    solving_times::Vector{Float64}
    objectives::Vector{Function}
    constraints::Vector{Function}
    history_solutions::AbstractVector{<:AbstractVector{<:Integer}} # keep track of solutions
    start_time::Float64

    function TspSolver(tsp::TspProblem, options)
        fit_scheme = EpsBoxDominanceFitnessScheme(fitness_scheme(tsp), options[:ϵ])
        archive = EpsBoxArchive(fit_scheme)
        evaluator = BlackBoxOptim.make_evaluator(tsp, archive, options)
        _objectives = solver_objectives(tsp)
        _constraints = map(cid -> solver_constraints(tsp, cid), constraints(tsp))
        return new(phase1::SolverStage, tsp, archive, evaluator, options, Float64[], _objectives, _constraints, Vector{Vector{UInt8}}())
    end
end

function time_limit(s::TspSolver)::Float64
    totaltime = options(s)[:MaxTime]
    timeleft = totaltime - (time() - s.start_time)
    if 60 > timeleft > 0
        timeleft = 60
    end
    return timeleft
end

function solver_constraints(tsp::TspProblem, cid::Symbol)
    if cid ∈ tsp |> categories |> keys
        return (m, c) -> tsp_add_category_constraints(m, c, cid)
    else
        throw(ArgumentError(cid, "no feature with name '$cid' could be found. Execution stops here. Revise argument for rerun."))
    end
end

function solver_objective(tsp::TspProblem, cid::Symbol)
    pool = TSP_NORMALIZED_EXPRESSIONS
    if haskey(pool, cid)
        return pool[cid]
    elseif cid ∈ tsp |> categories |> keys
        return (m, c) -> tsp_normalized_category_expression(m, c, cid)
    else
        throw(ArgumentError(cid, "no feature with name '$cid' could be found. Execution stops here. Revise argument for rerun."))
    end
end

function solver_objectives(tsp::TspProblem)
    return map(oid -> solver_objective(tsp, oid), objectives(tsp))
end

stage(s::TspSolver) = s.stage
BlackBoxOptim.problem(s::TspSolver) = s.problem
evaluator(s::TspSolver) = s.evaluator
options(s::TspSolver) = s.options
solving_times(s::TspSolver) = s.solving_times
objectives(s::TspSolver) = s.objectives
constraints(s::TspSolver) = s.constraints
history_solutions(s::TspSolver) = s.history_solutions
problem(s::TspSolver) = s.problem

@forward TspSolver.problem case

function fitness(s::TspSolver, indiv)
    push!(history_solutions(s), indiv)# implicitly add to archive
    return BlackBoxOptim.fitness(indiv, evaluator(s))
end

function set_phase1_fitness(s::TspSolver, indiv)
    fitness(s, indiv)
    options(s)[:OptSize] = sum(indiv)

    manual_cap = round(options(s)[:λᵤ] * numtests(problem(s)), RoundUp)
    options(s)[:Cap] = max(options(s)[:OptSize] * options(s)[:ζ], manual_cap)

    s.stage = phase2 # upgrade search stage
end

function solve_phase1(s::TspSolver)
    indiv = solve_min_size(case(s), options(s), constraints(s)) # get solution
    set_phase1_fitness(s, indiv)
end

function step!(s::TspSolver, ::Type{Val{phase1}})
    s.start_time = round_start = time()

    # find global optimal min size with 100% coverage
    try
        solve_phase1(s)
    catch e
        println(e)
        # timeout -> fake minimum fitness to the archive for smooth takedown 
        indiv = ones(Int8, numtests(problem(s)))
        set_phase1_fitness(s, indiv)

        # since time is up, no results will be reported, so nothing must be done here.
    end

    push!(solving_times(s), time() - round_start)
    return s
end

function step!(s::TspSolver, ::Type{Val{phase2}})

    cap = options(s)[:Cap]
    for obj_func in objectives(s)
        if obj_func == tsp_normalized_expression_min_size # minsize has been solved at this stage
            continue
        end

        round_start = time()
        timelimit = time_limit(s)
        if timelimit > 0
            relaxed = options(s)[:Relaxed]
            tsp = optimize!(case(s), cap, relaxed, options(s)[:Solver], obj_func, timelimit, constraints(s))

            push!(solving_times(s), time() - round_start)
            tsp |> model |> termination_status |> println
            if tsp |> model |> termination_status |> string == "OPTIMAL"
                indiv = solution_binary(tsp, relaxed, cap)
                fitness(s, indiv)
            end
        end
    end

    s.stage = phase3 # upgrade search stage
    return s
end

function randomized_solution_binary(tsp::TestSelectionProblemInstance, size::Real)
    size = trunc(Int64, size)
    vvv = value.(model(tsp)[:t])
    res = Vector{Float64}(undef, length(vvv))
    foreach(i -> res[i] = vvv[i], eachindex(res))
    highest_bunch = sortperm(res, rev=true)[1:size]
    result = zeros(UInt8, length(res))
    foreach(i -> result[i] = one(UInt8), highest_bunch)
    return result
end

function solution_binary(tsp::TestSelectionProblemInstance, relaxed::Bool, size::Real)
    if relaxed
        println("relaxed!")

        if size < 0
            vvv = value.(model(tsp)[:t])
            res = Vector{Float64}(undef, length(vvv))
            foreach(i -> res[i] = vvv[i], eachindex(res))
            size = sum(res[ res .> 0 ])
        end

        return randomized_solution_binary(tsp, size)
    else
        println("strict!")
        return solution_binary(tsp)
    end
end

function step!(s::TspSolver, ::Type{Val{phase3}})

    round_start = time()

    relaxed = options(s)[:Relaxed]
    rand_cap = rand(options(s)[:OptSize]:options(s)[:Cap])
    tsp = create_tsp_second_order(case(s), rand_cap, relaxed, options(s)[:MaxTime], options(s)[:Solver])

    foreach(f -> f(model(tsp), case(s)), constraints(s)) # activate all constraints

    N = evaluator(s) |> fitness_scheme |> numobjectives
    objective_idx = sample(1:N,rand(2:N), replace = false) # draw the involved objectives
    obj_funcs = objectives(s)[objective_idx]
    obj_exprs = map(f -> f(model(tsp), case(s)), obj_funcs)
    rand_weights = normalize!(rand(length(obj_exprs)), 1) # random objective weights to vary landscape
    expr = sum(rand_weights .* obj_exprs)

    @objective(model(tsp), Min, expr)

    timelimit = time_limit(s)
    if timelimit > 0
        optimize!(tsp, timelimit)

        tsp |> model |> termination_status |> println
        push!(solving_times(s), time() - round_start)

        if tsp |> model |> termination_status |> string == "OPTIMAL"
            indiv = solution_binary(tsp, relaxed, rand_cap)
            fitness(s, indiv)
        end
    end

    return s
end

BlackBoxOptim.step!(s::TspSolver) = step!(s, Val{stage(s)})

function convex_solver(problem::TspProblem, options::Parameters = EMPTY_PARAMS)
    opts = chain(BlackBoxOptim.EMPTY_PARAMS, options)
    return TspSolver(problem, opts)
end
