include("bbo/random_mop_tsp.jl")
include("bbo/borg/borg_setup.jl")
include("bbo/convex_solver/convex_solver.jl")
include("systems/testcase_selection_optimization_problem.jl")
include("bbo/experiment.jl")

BlackBoxOptim.add_method_to_bbo(:convex_solver, convex_solver) # making our convex solver algorithm available to BlackBoxOptim

# --------------------- alg setup ------------------------

alg_setup = ParamsDict(
                     :λᵤ => .1,                     # test set size limit expressed as in article (in percentage of original test set size)
                     :ζ => 1.3,                     # slack factor expressed as in article (multiplier of min size set)
                     :ϵ => .1,                      # Borg specific parameter for pareto archive granularity
                     :Solver => Cbc.Optimizer,      # Solver alternatives from study are Gurobi.Optimizer and Clp.Optimizer
                     :AlgorithmName => "Cbc",       # just the name of the solver, OBS set correctly here
                     :Relaxed => false,             # whether to solve as relatex version (linear program - approximate) or integer program (optimal)
                     :MaxTime => 60,                # runtime for entire experiment (depends on problem size, but likely at least above 60 s)
                     :Method => :convex_solver)     # to tell BlackBoxOptim to use MC-TOA.
                     
# --------------------- problem definition ------------------------
# fictious failure rate distribution (based on experience). For some companies lower, i.e. Normal(0.03, 0.02). Not used for the study, but can be used to create artificial input.
H_FAILURE_RATE_DISTRIBUTION = truncated(Normal(0.08, 0.05), 0.0, 1.0)
E = sparse(rand(Bool, 50, 100))         # randomly create execution matrix
_failrates = rand(H_FAILURE_RATE_DISTRIBUTION, size(E, 2))  # randomly create failure distribution

mod1 = rand(Bool, 100) # mututal exclusive module coverage modeled by binary flags
mod2 = map(!, mod1)
c_df = DataFrame(:module1 => mod1, :module2 => mod2)
o_df = DataFrame(:feature1 => rand(Bool, 100), :feature2 => rand(Bool, 100), :feature3 => rand(Bool, 100))

criteria = Dict{Symbol, DataFrame}(:objective_categorization => o_df, :constraint_categorization => c_df)

_problem = TspProblem(Case(sparse(E),
                            criteria,          # all criteria
                            _failrates),       # failure distribution
                            Symbol[:objective_categorization],          # choice of categories as objectives
                            Symbol[:constraint_categorization])          # choices of categories as constraints
instance = TspInstance(_problem, "example_problem") # give instance a name

# --------------------- run experiment ------------------------

opt_name = alg_setup[:AlgorithmName]
exec_time = alg_setup[:MaxTime]
println("===============optimizer $opt_name on problem $(id(instance)), total time = $exec_time)========================")
experiment = TspExperiment(instance, opt_name)  # create experiment

res = bboptimize(problem(instance); alg_setup...)   # run optimization
tsp_results = TspExperimentResult(experiment, res)  # retrieve results

#printing the indexes for one of the resulting regression test sets from the pareto front:
nondominatedtestset = pf(tsp_results)[1].nzind
