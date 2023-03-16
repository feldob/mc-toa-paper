include("../../model/case.jl")
include("../../model/fitness.jl")

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

abstract type Criterion end

struct CompleteCrtierion
    objective_expression::Function

end

expression(c::Criterion) = c.objective_expression
normalizedexpression()

min_size_criterion = Criterion(
                            norm_obj_expr = (m,c) -> tsp_size(m[:t]) / numtests(c),
                            )

as_objective(m::Model,c::Case) = tsp_set_objective!(m, case, c.objective_expression)
