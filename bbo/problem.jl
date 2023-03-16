include("../model/case.jl")
include("../model/fitness.jl")

using Lazy

# mop, all minimizing
function tsp_mop(testsactive::AbstractVector{<:Integer}, c::CachedCase)
    _coverage = tsp_coverage(testsactive, activationvectors(c), numfuncs(c))
    if _coverage < 1.0
        return convert(Float64, numtests(c)), 0.0, 0.0
    end

    _size = convert(Float64, tsp_size(testsactive))
    _exercises = convert(Float64, tsp_exercises(testsactive, funcspertest(c)))
    _failrate = tsp_failrate(testsactive, failrates(c))

    return _size, -_exercises, -_failrate
end

using BlackBoxOptim

function sum_up(fitness)
    f1 = 1-1/fitness[1]
    f2 = 1/fitness[2]
    f3 = 1/fitness[3]

    return f1 + f2 + f3
end

struct TspObjectives
    fitness_closure::Function # 2 inputs, solution and case
    objectives::Vector{Symbol}
    aggregator::Function
end

fitness_closure(o::TspObjectives) = o.fitness_closure
objectives(o::TspObjectives) = o.objectives
aggregator(o::TspObjectives) = o.aggregator

const REGULAR_TSP_OBJECTIVES = [ :minSize, :maxExercise, :maxFailRate ]

function category_fitness_aggregation_wrapper(numregular::Int64, fitness_vector, aggregator::Function)
    _reg_agg_fitness = aggregator(fitness_vector)
    _cat_fitnesses = 1 .- map(i -> fitness_vector[i], numregular+1:length(fitness_vector))
    return _reg_agg_fitness + sum(_cat_fitnesses)
end

function category_fitness_wrapper(_categories::Vector{Symbol}, fitness_func::Function, indiv, case)
    cat_fitness = map(c -> tsp_criterion_coverage(indiv, categories(case)[c]), _categories)
    regular_mop_fitness = fitness_func(indiv, case)
    n = length(cat_fitness) + length(regular_mop_fitness)
    NTuple{n, Float64}(vcat(regular_mop_fitness..., cat_fitness...))
end

function tspObjectives(_categories::Vector{Symbol})
    if _categories |> isempty
        return TspObjectives(tsp_mop, REGULAR_TSP_OBJECTIVES, sum_up)
    end

    _objectives = vcat(REGULAR_TSP_OBJECTIVES, _categories...)
    _numregular = length(REGULAR_TSP_OBJECTIVES)
    _fitness_closure = (indiv, c) -> category_fitness_wrapper(_categories, tsp_mop, indiv, c)
    _aggregator = (fitness -> category_fitness_aggregation_wrapper(_numregular, fitness, sum_up))
    return TspObjectives(_fitness_closure, _objectives, _aggregator)
end

mutable struct TspProblem{FS<:FitnessScheme} <: OptimizationProblem{FS}
    fitness_func::Function
    fitness_scheme::FS
    search_space::SearchSpace
    case::CachedCase
    objectives::Vector{Symbol}
    constraints::Vector{Symbol}

    function TspProblem(case::CachedCase,
                            cat_objectives::Vector{Symbol}=Symbol[],
                            cat_constraints::Vector{Symbol}=Symbol[])
        @assert intersect(cat_objectives, cat_constraints) |> isempty
        o::TspObjectives = tspObjectives(cat_objectives)
        tsp_fitness(t) = fitness_closure(o)(t, case)
        DType = UInt8
        _bounds = fill((zero(DType), one(DType)), numtests(case))
        _search_space = RectSearchSpace(_bounds, DType)

        _objectives = objectives(o)
        _fitness_scheme = ParetoFitnessScheme{length(_objectives)}(is_minimizing = true, aggregator = aggregator(o))

        new{typeof(_fitness_scheme)}(tsp_fitness, _fitness_scheme, _search_space, case, _objectives, cat_constraints)
   end
end

TspProblem(c::Case, os=Symbol[], cs=Symbol[]) = TspProblem(CachedCase(c), os, cs)

name(tsp::TspProblem) = "Test Selection Problem"
case(tsp::TspProblem) = tsp.case
objectives(tsp::TspProblem) = tsp.objectives
constraints(tsp::TspProblem) = tsp.constraints

BlackBoxOptim.fitness(input::AbstractVector, p::TspProblem) = p.fitness_func(input)

@forward TspProblem.case numtests, categories
@forward TspProblem.search_space bounds

if (@isdefined SPARSE) && SPARSE
    include("sparse.jl")
end

tspMutation(tsp, options) = BlackBoxOptim.BinaryFlipMutation(.01)
tspEmbedding(tsp, options) = BlackBoxOptim.NoEmbedding()

function tspCrossover(tsp::TspProblem, options::AbstractDict)
    points = trunc(Int64, numtests(tsp) * .01)
    return CrossoverOperator[ BlackBoxOptim.MultiPointCrossover(points) ]
end
