# to keep track of Borg history, must be cleaned up by "collector"
function reset_heuristic_stats()
    global cached_solving_times = Vector{Float64}()
    global cached_history_solutions = Vector{AbstractVector{<:Integer}}()
    global cached_archive_snapshot = Set{AbstractVector{BlackBoxOptim.FrontierIndividual}}()
end

function trace_results_heuristic(optcontroller::BlackBoxOptim.OptRunController)
    generation = BlackBoxOptim.num_steps(optcontroller)
    elapsed = BlackBoxOptim.elapsed_time(optcontroller)
    cached_pf = BlackBoxOptim.pareto_frontier(BlackBoxOptim.archive(optcontroller.optimizer))
    new_candidates = Set(cached_pf)

    #calc diff to last archive, and add with current timestamp the solutions
    diff_candidates = setdiff!(new_candidates, cached_archive_snapshot)

    # store the diff as history entries
    for c in diff_candidates
        push!(cached_solving_times, elapsed)
        push!(cached_history_solutions, c.params)
    end

    # reset borg archive
    global cached_archive_snapshot = new_candidates
end

# IMPORTANT THAT THIS IS DONE INITIALLY!
reset_heuristic_stats()
