#OBS reset requires restart of environment
global SPARSE = false

include("population.jl")
include("../history_tracing.jl")

borg_setup = ParamsDict(
                    :Sparse => SPARSE,
                    :ϵ => .1,
                    :PopulationStrategy => RandomPopulationStrategy(),
                    :CallbackFunction => trace_results_heuristic,
                    :CallbackInterval => 5.0,
                    :PopulationSize => 10,
                    :Ntransient => 1,
                    :ProportionInitialActivations => .3,
                    :CrossoverOperator => tspCrossover,
                    :MutationOperator => tspMutation,
                    :EmbeddingOperator => tspEmbedding,
                     :AlgorithmName => "borg",
                    :Method => :borg_moea)

borg_greedy_setup = deepcopy(borg_setup)
borg_greedy_setup[:AlgorithmName] = "borg_greedy"
borg_greedy_setup[:PopulationStrategy] = GreedyOptCoveragePopulationStrategy()
