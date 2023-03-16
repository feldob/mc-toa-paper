include("history_tracing.jl")

abstract type HistoryOptimizer <: SteppingOptimizer end
abstract type HeuristicHistoryOptimizer <: HistoryOptimizer end

struct HistoryOptimizerOutput <: BlackBoxOptim.MethodOutput
    solving_times::AbstractVector{Float64}
    history_solutions::AbstractVector{<:AbstractVector{<:Integer}}

    HistoryOptimizerOutput(s::HistoryOptimizer)= new(solving_times(s), history_solutions(s))
    HistoryOptimizerOutput(s::HeuristicHistoryOptimizer)= new(cached_solving_times, cached_history_solutions)
end

# The generic variant that simply wraps the elapsed time into a vector
solving_times(r::BlackBoxOptim.OptimizationResults, ::BlackBoxOptim.MethodOutput) = [ BlackBoxOptim.elapsed_time(r) ]
solving_times(r::BlackBoxOptim.OptimizationResults) = solving_times(r, r.method_output)

history_solutions(r::BlackBoxOptim.OptimizationResults, ::BlackBoxOptim.MethodOutput) = Vector{Vector{UInt8}}()
history_solutions(r::BlackBoxOptim.OptimizationResults) = history_solutions(r, r.method_output)

BlackBoxOptim.MethodOutput(s::HistoryOptimizer) = HistoryOptimizerOutput(s)
solving_times(r::BlackBoxOptim.OptimizationResults, o::HistoryOptimizerOutput) = o.solving_times
history_solutions(r::BlackBoxOptim.OptimizationResults, o::HistoryOptimizerOutput) = o.history_solutions
