include("../problem.jl")
include("../../heuristics.jl")

# Override of frontierindividual which enforces sparse entries in the archive only
BlackBoxOptim.FrontierIndividual(fitness::F,
               params, tag, num_fevals, n_restarts, timestamp=time()) where F =
    BlackBoxOptim.FrontierIndividual{F}(fitness, sparse(params), tag, num_fevals, n_restarts, timestamp)

abstract type PopulationStrategy end
struct RandomPopulationStrategy <: PopulationStrategy end
struct RandomCoveragePopulationStrategy <: PopulationStrategy end
struct GreedyCoveragePopulationStrategy <: PopulationStrategy end
struct GreedyOptCoveragePopulationStrategy <: PopulationStrategy end

function population(c::AbstractCase, popsize::Integer, _issparse::Bool, subset_creator::Function)
    candidates = [ to_solution(subset_creator(c), numtests(c)) for _ in 1:popsize]
    pop_matrix = reduce(hcat, candidates)
    return _issparse ? pop_matrix : Matrix(pop_matrix)
 end

function population(::GreedyCoveragePopulationStrategy, c::CachedCase, _popsize::Integer, options)
    return population(c, _popsize, options[:Sparse], greedy_full_coverage_subset_fast)
end

# creates randomized initial population of greedy candidates to have a good variability in the inital population
function population(::GreedyOptCoveragePopulationStrategy, c::CachedCase, _popsize::Integer, options)
    pop_matrix = population(c, _popsize, options[:Sparse], (_c) -> greedy_full_coverage_subset_slow_opt(_c, true))
    #pop_matrix = population(RandomCoveragePopulationStrategy(), c, _popsize, options) # random with coverage init
    pop_matrix[:,1] = to_solution(greedy_full_coverage_subset_slow_opt(c), numtests(c)) # add opt greedy candidate
    return options[:Sparse] ? pop_matrix : Matrix(pop_matrix)
end

function population(::RandomCoveragePopulationStrategy, c::CachedCase, _popsize::Integer, options)
    return population(c, _popsize, options[:Sparse], random_full_coverage_subset)
end

function population(::RandomPopulationStrategy, c::CachedCase, _popsize::Integer, options)
    percent_initial_activations = options[:ProportionInitialActivations]

    activation_vector_length = _popsize * numtests(c)
    ones = randsubseq(1:activation_vector_length, percent_initial_activations)

    if options[:Sparse]
        pop_vector = sparsevec(ones, 0x01, activation_vector_length)
        return sparse(reshape(pop_vector, numtests(c), _popsize))
    else
        pop_vector = zeros(UInt8, activation_vector_length)
        pop_vector[ones] .= 1
        return reshape(pop_vector, numtests(c), _popsize)
    end
 end

function BlackBoxOptim.population(tsp::TspProblem, options::AbstractDict{Symbol,Any})
    fs = fitness_scheme(tsp)
    _nafitness = nafitness(IndexedTupleFitness{numobjectives(fs),fitness_eltype(fs)})
    strategy = options[:PopulationStrategy]
    popsize = options[:PopulationSize]
    ntransient = options[:Ntransient]
    pop_matrix = population(strategy, case(tsp), popsize + ntransient, options)
    return FitPopulation(pop_matrix, _nafitness, ntransient=ntransient)
end

function BlackBoxOptim.borg_moea(problem::TspProblem, options::AbstractDict{Symbol,Any})
    options[:Population] = get(options, :Population, BlackBoxOptim.population)
    return invoke(BlackBoxOptim.borg_moea, Tuple{OptimizationProblem, AbstractDict{Symbol,Any}}, problem, options)
end
