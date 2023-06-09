include("convex_solver.jl")

BlackBoxOptim.add_method_to_bbo(:convex_solver, convex_solver)

gurobi_setup = ParamsDict(
                     #:CallbackInterval => 2.0,
                     :λᵤ => .1,
                     :ζ => 1.3,
                     :ϵ => .1,
                     :Solver => Gurobi.Optimizer,
                     :Relaxed => false,
                     :AlgorithmName => "gurobi",
                     :Method => :convex_solver)

gurobi_setup_relaxed = copy(gurobi_setup)
gurobi_setup_relaxed[:Relaxed] = true
gurobi_setup_relaxed[:AlgorithmName] = "gurobi_relaxed"

cbc_setup = copy(gurobi_setup)
cbc_setup[:Solver] = Cbc.Optimizer
cbc_setup[:AlgorithmName] = "Cbc"

cbc_setup_relaxed = copy(gurobi_setup_relaxed)
cbc_setup_relaxed[:Solver] = Cbc.Optimizer
cbc_setup_relaxed[:AlgorithmName] = "Cbc_relaxed"

clp_setup = copy(gurobi_setup_relaxed)
clp_setup[:Solver] = Clp.Optimizer
clp_setup[:Relaxed] = true
clp_setup[:AlgorithmName] = "Clp"
