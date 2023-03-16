include("experiments_results.jl")

function stats_exp1_sys1(_runtime)
    _objectives_tomax = [(x) -> -x, identity, identity, identity, identity]
    _objectives_labels = ["size", "exercises", "failrate", "features", "coverage"]

    _algids = ["gurobi", "gurobi_relaxed", "Cbc", "Cbc_relaxed", "Clp", "borg_greedy", "random_greedy"]

    _expnr_c = "original_system_1_constraint"
    #create_stats(_expnr_c,  _objectives_labels, _objectives_tomax, _runtime, _algids)
    plot_2ds(_expnr_c, "sys1_constr_2d.pdf", "System 1, constraints", [1,2,4,5], 2, 17, 32, 47, 62)

    _expnr_o = "original_system_1_objective"
    #create_stats(_expnr_o, _objectives_labels, _objectives_tomax, _runtime, _algids)
    plot_2ds(_expnr_o, "sys1_obj_2d.pdf", "System 1, objective", [1,2,4,5,7,9], 47, 76, 92, 62, 3, 17, 39)

    #plot_runtimes(_expnr_c, _expnr_o, "sys1_runtimes.pdf")
end

function stats_exp1_sys2(_runtime)
    _objectives_tomax = [(x) -> -x, identity, identity, identity, identity]
    _objectives_labels = ["size", "exercises", "failrate", "morphologies", "coverage"]
    _algids = ["gurobi", "gurobi_relaxed", "Cbc", "Cbc_relaxed", "Clp", "borg_greedy", "random_greedy"]

    _expnr_c = "original_system_2_constraint"
    #create_stats(_expnr_c, _objectives_labels, _objectives_tomax, _runtime, _algids)
    plot_2ds(_expnr_c, "sys2_constr_2d.pdf", "System 2, constraints", [1,2,7], 2, 17, 32, 47, 62)

    _expnr_o = "original_system_2_objective"
    #create_stats(_expnr_o, _objectives_labels, _objectives_tomax, _runtime, _algids)
    plot_2ds(_expnr_o, "sys2_obj_2d.pdf", "System 2, objective", [1,2,7], 47, 76, 92, 62, 3, 17, 39)

    #plot_runtimes(_expnr_c, _expnr_o, "sys2_runtimes.pdf")
end

function stats_exp_series_1()
    _runtime = 300
    stats_exp1_sys1(_runtime)
    stats_exp1_sys2(_runtime)
end

#stats_exp_series_1()