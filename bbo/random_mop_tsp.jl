include("problem.jl")
include("../heuristics.jl")
include("history_optimizer.jl")

using BlackBoxOptim

abstract type HeuristicTspStrategy end

struct RandomTspStrategy <: HeuristicTspStrategy
    RandomTspStrategy(c::AbstractCase) = new()
end

mutable struct GreedyTspStrategy <: HeuristicTspStrategy
    case::AbstractCase
    first::Bool

    function GreedyTspStrategy(c::AbstractCase)
        new(c, true)
    end
end

case(g::GreedyTspStrategy) = g.case
thefirst(g::GreedyTspStrategy) = g.first
toheuristic(g::GreedyTspStrategy) = g.first = false

mutable struct RandomMopTsp <: HeuristicHistoryOptimizer
    problem::TspProblem
    evaluator::BlackBoxOptim.Evaluator
    strategy::HeuristicTspStrategy
    options

    function RandomMopTsp(tsp::TspProblem, options)
        fit_scheme = EpsBoxDominanceFitnessScheme(fitness_scheme(tsp), options[:ϵ])
        archive = EpsBoxArchive(fit_scheme)
        evaluator = BlackBoxOptim.make_evaluator(tsp, archive, options)
        strategy = options[:Strategy]
        return new(tsp, evaluator, strategy(case(tsp)), options)
    end
end

BlackBoxOptim.problem(s::RandomMopTsp) = s.problem
evaluator(s::RandomMopTsp) = s.evaluator
options(s::RandomMopTsp) = s.options
case(s::RandomMopTsp) = case(s.problem)
strategy(r::RandomMopTsp) = r.strategy

function greedy_full_coverage_subset_slow_opt_prep(s::GreedyTspStrategy)::Vector{<:Integer}
    randomized = !thefirst(s)

    if thefirst(s)
        toheuristic(s)
    end

    return greedy_full_coverage_subset_slow_opt(case(s), randomized)
end

function apply_strategy!(g::GreedyTspStrategy, r::RandomMopTsp)
    return greedy_full_coverage_subset_slow_opt_prep(g)
end

function apply_strategy!(::HeuristicTspStrategy, r::RandomMopTsp)
    return random_full_coverage
end

function BlackBoxOptim.step!(r::RandomMopTsp)
    result = apply_strategy!(strategy(r), r)
    indiv = to_solution(result, numtests(case(r)))
    BlackBoxOptim.fitness(indiv, evaluator(r))
    return r
end

BlackBoxOptim.archive(alg::RandomMopTsp) = alg.evaluator.archive

function random_mop_tsp(problem::TspProblem, options::Parameters = EMPTY_PARAMS)
    opts = chain(BlackBoxOptim.EMPTY_PARAMS, options)
    return RandomMopTsp(problem, opts)
end

BlackBoxOptim.add_method_to_bbo(:random_mop_tsp, random_mop_tsp)

random_mop_tsp_setup = ParamsDict(
                     :ϵ => .1,
                     :AlgorithmName => "random",
                     :CallbackFunction => trace_results_heuristic,
                     :CallbackInterval => 5.0,
                     :Strategy => RandomTspStrategy,
                     #:NThreads => Threads.nthreads() - 1,
                     :Method => :random_mop_tsp)

random_greedy_mop_tsp_setup = deepcopy(random_mop_tsp_setup)
random_greedy_mop_tsp_setup[:Strategy] = GreedyTspStrategy
random_greedy_mop_tsp_setup[:AlgorithmName] = "random_greedy"
