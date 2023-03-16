include("../../model/case.jl")
include("../../model/fitness.jl")

using JuMP, Cbc, Gurobi, Clp

# ------------------------PROBLEM DEFINITIONS-------------------- #

tsp_expression_min_size(m::Model, ::AbstractCase) = tsp_size(m[:t])
tsp_expression_max_exercise(m::Model, c::AbstractCase) = tsp_exercises(m[:t], funcspertest(c))
tsp_expression_max_failrate(m::Model, c::AbstractCase) = tsp_failrate(m[:t], failrates(c))

tsp_normalized_expression_min_size(m, c) = tsp_expression_min_size(m,c) / numtests(c)
tsp_normalized_expression_max_exercise(m, c) = 1 - tsp_expression_max_exercise(m,c) / maxexercises(c)
tsp_normalized_expression_max_failrate(m, c) = 1 - tsp_expression_max_failrate(m,c) / maxfailrate(c)

function tsp_set_objective!(m::Model, case::AbstractCase, expr_func::Function)
    expr = expr_func(m, case)
    @objective(m, Min, expr)
end

tsp_objective_min_size(m, case) = tsp_set_objective!(m, case, tsp_normalized_expression_min_size)
tsp_objective_max_exercise(m, case) = tsp_set_objective!(m, case, tsp_normalized_expression_max_exercise)
tsp_objective_max_failrate(m, case) = tsp_set_objective!(m, case, tsp_normalized_expression_max_failrate)

TSP_NORMALIZED_EXPRESSIONS = Dict{Symbol, Function}(
                            :minSize => tsp_normalized_expression_min_size,
                            :maxExercise => tsp_normalized_expression_max_exercise,
                            :maxFailRate => tsp_normalized_expression_max_failrate)

TSP_OBJECTIVES = Dict{Symbol, Function}(
                            :minSize => tsp_objective_min_size,
                            :maxExercise => tsp_objective_max_exercise,
                            :maxFailRate => tsp_objective_max_failrate)

function tsp_normalized_category_expression(m, case, cid::Symbol)
    category_df = categories(case)[cid]
    slack_var_ids = map(p -> slack_var_name(cid, p), names(category_df))
    slack_vars = map(v_id -> variable_by_name(m, v_id), slack_var_ids)

    # Ensure that binary slack variables activate only if property of category is covered, i.e. ≥ 1
    # Constraint: sum (property_vector .* solution) - EXTREME_BOUNDARY * category_property_slack ≤ 0
    for idx in eachindex(slack_vars)
        p = category_df[:,idx]
        slack_var = slack_vars[idx]
        @constraint(m, sum(p[i] * m[:t][i] for i in eachindex(p)) - sum(p) * slack_var ≤ 0)
    end

    # Originally Max objective: 1/|properties| * sum(slack_property) <<-- always in range [0,1]
    # here expressed as "Min objective" in normalized form to be combined with other objectives : 1 - OBJ
    @expression(m, 1 - (1/length(slack_vars)) .* sum(slack_vars))
end

function tsp_add_category_constraints(m, c, cid::Symbol)
    category_df = categories(c)[cid]

    for property in names(category_df)
        p = category_df[:,property]
        @constraint(m, sum(p[i] * m[:t][i] for i in eachindex(p)) ≥ 1)
    end
end

# --------------------------------------------------------------- #

function my_callback_function(cb_data)
    #"start callback" |> println
    #global gl_cb_data = cb_data
    #cb_data |> println
    #x_val = callback_value(cb_data, t)
    #status = MOI.submit(
    #    model, MOI.HeuristicSolution(cb_data), [t], [floor(Int, x_val)]
    #)
    #println("I submitted a heuristic solution, and the status was: ", status)
end

mutable struct TestSelectionProblemInstance
    model::Model
    case::AbstractCase

    function TestSelectionProblemInstance(case, opt=Gurobi.Optimizer, timelimit = 60)
        m = Model(opt)
        set_time_limit_tsp(m, timelimit)
        set_silent(m)
        #MOI.set(m, MOI.HeuristicCallback(), my_callback_function)
        new(m, case)
    end
end

model(t::TestSelectionProblemInstance) = t.model
case(t::TestSelectionProblemInstance) = t.case

struct TestSelectionProblemInstanceResult end

slack_var_name(cat::Symbol, prop::String)::String = "$(string(cat))_$(prop)"

function variables!(tsp::TestSelectionProblemInstance, relaxed=false::Bool)
    var_range = 1:numtests(case(tsp))

    if relaxed
        println("relaxed!")
        @variable(model(tsp), 0 ≤ t[var_range] ≤ 1)
    else
        println("strict!")
        @variable(model(tsp), t[var_range], Bin)
    end

    # add slack variables per category and property combination with unique names
    for (c, df) in (pairs ∘ categories ∘ case)(tsp)
        for p in names(df)
            if relaxed
                v = @variable(model(tsp), base_name = slack_var_name(c, p))
                set_lower_bound(v, 0.0)
                set_upper_bound(v, 1.0)
            else
                @variable(model(tsp), base_name = slack_var_name(c, p), binary = true)
            end
        end
    end
end

function objective!(tsp::TestSelectionProblemInstance, objective::Symbol=:minSize)
    TSP_OBJECTIVES[objective](model(tsp), tsp.case)
end

function constraint_guarantee_all_functions_covered!(tsp::TestSelectionProblemInstance)
    m = model(tsp)
    c = case(tsp)

    fs = map(f -> activationmatrix(c)[f,:], 1:numfuncs(c))
    foreach(f -> @constraint(m,  sum(f[i] * m[:t][i] for i in eachindex(f)) ≥ 1), fs)
end

function constraint_cover_category!(tsp::TestSelectionProblemInstance, category::Vector{<:Bool})
    m = model(tsp)

    @constraint(m, sum(m[:t] * category) ≥ 1)
end

function constraint_cap_tests_suite!(tsp::TestSelectionProblemInstance, cap::Number)
    @constraint(model(tsp), sum(model(tsp)[:t]) ≤ cap)
end

function set_time_limit_tsp(model, timelimit)
    if solver_name(model) == "COIN Branch-and-Cut (Cbc)"
        set_optimizer_attribute(model, "seconds", timelimit)
    else
        set_time_limit_sec(model, timelimit)
    end
    "time limit set to $timelimit seconds" |> println
end

optimize!(tsp::TestSelectionProblemInstance) = JuMP.optimize!(model(tsp))
function optimize!(tsp::TestSelectionProblemInstance, timelimit::Real)
    set_time_limit_tsp(model(tsp), timelimit)
    optimize!(tsp)
end

function solution(tsp::TestSelectionProblemInstance)
    collect(1:numtests(case(tsp)))[value.(model(tsp)[:t]) .== 1]
end

function solution_binary(tsp::TestSelectionProblemInstance)
    vvv = value.(model(tsp)[:t])
    res = Vector{UInt8}(undef, length(vvv))
    #Cbc can make rounding errors, which is why rounding necessary here, generally not required though.
    foreach(i -> res[i] = round(vvv[i]), eachindex(res))
    return res
end

# this means for relaxed version the outcome is not optimal, it can be better than that! (consider in results)
function solve!(tsp::TestSelectionProblemInstance, objective_type::Symbol, relaxed::Bool = false, constraints...)
    variables!(tsp, relaxed)
    objective!(tsp, objective_type)
    foreach(c -> c(tsp), constraints)
    optimize!(tsp)
end

function solve_indices!(tsp::TestSelectionProblemInstance, objective_type::Symbol, constraints...)
    solve!(tsp, objective_type, constraints...)
    return solution(tsp)
end

function solve_binary!(tsp::TestSelectionProblemInstance, cap::Real, objective_type::Symbol, relaxed::Bool = false, constraints...)
    solve!(tsp, objective_type, relaxed, constraints...)

    if ! (tsp |> model |> has_values)
        throw(DomainError("No solution could be found in the setout time window."))
    end

    tsp |> model |> termination_status |> string |> println

    return solution_binary(tsp, relaxed, cap)
end

# relaxed case: encodes variables as floats, OBS result is not optimal, its an approximation!
function solve_min_size(tsp::TestSelectionProblemInstance, relaxed::Bool, _constraints::Vector{Function} = Function[])

    _suitable_constraints = map(c -> ((_tsp) -> c(model(_tsp), case(_tsp))), _constraints)

    return solve_binary!(tsp, -1, :minSize, relaxed, # The size will never activate, so use an "invalid" entry to raise error in case it does
                    constraint_guarantee_all_functions_covered!, _suitable_constraints...)
end

function solve_min_size(case::AbstractCase, options, _constraints::Vector{Function})
    solver = options[:Solver]
    relaxed = options[:Relaxed]

    half_of_budget = div(options[:MaxTime], 2)
    tsp = TestSelectionProblemInstance(case, solver, half_of_budget)
    return solve_min_size(tsp, relaxed, _constraints)
end

function create_tsp_second_order(case::AbstractCase, cap::Real, relaxed::Bool, opt = Gurobi.Optimizer)
    tsp = TestSelectionProblemInstance(case, opt)
    cgts! = (tsp) -> constraint_cap_tests_suite!(tsp, cap)

    constraints = [ constraint_guarantee_all_functions_covered!, cgts! ]

    variables!(tsp, relaxed)
    foreach(c -> c(tsp), constraints)

    return tsp
end

function optimize!(case::AbstractCase, cap::Real, relaxed::Bool, solver, expr_func::Function, timelimit::Real, constraint_funcs::Vector{Function}=[])
    tsp = create_tsp_second_order(case, cap, relaxed, solver)

    foreach(cf -> cf(model(tsp), case), constraint_funcs) # install constraints

    expr = expr_func(model(tsp), case)  # create objective expression
    @objective(model(tsp), Min, expr)   # install objective expression
    set_time_limit_tsp(model(tsp), timelimit)
    optimize!(tsp)
    return tsp
end

function solve_second(case::AbstractCase, cap::AbstractFloat, objective::Symbol, opt = Gurobi.Optimizer)
    tsp = TestSelectionProblemInstance(case, opt)
    cgts! = (tsp) -> constraint_cap_tests_suite!(tsp, cap)

    return solve_binary!(tsp, objective, cap,
                    constraint_guarantee_all_functions_covered!,
                    cgts!)
end

function solve_max_exercise(case::AbstractCase, cap::AbstractFloat, opt = Gurobi.Optimizer)
    solve_second(case::AbstractCase, cap::AbstractFloat, :maxExercise, opt)
end

function solve_max_exercise(case::AbstractCase, options)
    solve_max_exercise(case::AbstractCase, options[:Cap])
end

function solve_max_failrates(case::AbstractCase, cap::AbstractFloat, opt = Gurobi.Optimizer)
    solve_second(case::AbstractCase, cap::AbstractFloat, :maxFailRate, opt)
end

function solve_max_failrates(case::AbstractCase, options)
    solve_max_failrates(case::AbstractCase, options[:Cap])
end
