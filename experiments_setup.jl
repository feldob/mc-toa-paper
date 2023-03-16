include("bbo/random_mop_tsp.jl")
include("bbo/borg/borg_setup.jl")
include("bbo/convex_solver/convex_solver_setup.jl")
include("io.jl")

include("systems/testcase_selection_optimization_problem.jl")

NO_CRITERIA = Symbol[]

_runtime = 300
_runs = 15
_runs_exp2 = 5

_algs_setup_objective = [ gurobi_setup, gurobi_setup_relaxed, random_greedy_mop_tsp_setup, borg_greedy_setup, clp_setup, cbc_setup, cbc_setup_relaxed ]
_algs_setup_constraint = [ gurobi_setup, gurobi_setup_relaxed, cbc_setup,  cbc_setup_relaxed, clp_setup ]

function setup_experiment_for_huawei_system_2(exp_id::String)
    case_file = "systems/system2_testcase_2_function_call_chain.json"
    exp_dir = experiment_dir(exp_id)

    tsp_ex = TestcaseSelectionOptimizationProblem(case_file, "systems/system2_testcase_info.json")

    df_m = DataFrame()
    for fid in 1:size(tsp_ex.hwmorphologyByTestcase, 1)
        df_m[!, "M$fid"] = tsp_ex.hwmorphologyByTestcase[fid, :]
    end

    df_ccd3 = DataFrame()
    for fid in 1:size(tsp_ex.testfeatureByTestcase, 1)
        df_ccd3[!, "F$fid"] = tsp_ex.testfeatureByTestcase[fid, :]
    end

    tsp_ex = TestcaseSelectionOptimizationProblem(case_file, "systems/system2_testcase_info_featuredepth4.json")

    df_ccd4 = DataFrame()
    for fid in 1:size(tsp_ex.testfeatureByTestcase, 1)
        df_ccd4[!, "F$fid"] = tsp_ex.testfeatureByTestcase[fid, :]
    end

    cats = Dict{Symbol, DataFrame}( :morphologies => df_m,
                                    :callchaindepth3 => df_ccd3,
                                    :callchaindepth4 => df_ccd4)

    system_case = Case(tsp_ex.functionByTestcase, cats, tsp_ex.failRatePerTestcase)
    system_problem = TspProblem(system_case, Symbol[], [ :morphologies, :callchaindepth3 ])
    system_instance = TspInstance(system_problem, "System 2", 1)
    save(system_instance, exp_dir)
end

function setup_experiment_for_huawei_system_1(exp_id::String)
    exp_dir = experiment_dir(exp_id)

    tsp_ex = TestcaseSelectionOptimizationProblem("systems/system1_testcase_2_function_call_chain.json", "systems/system1_testcase_info.json")

    df = DataFrame()
    for fid in 1:size(tsp_ex.testfeatureByTestcase, 1)
        df[!, "F$fid"] = tsp_ex.testfeatureByTestcase[fid, :]
    end

    cats = Dict{Symbol, DataFrame}( :features => df)
    system_case = Case(tsp_ex.functionByTestcase, cats, tsp_ex.failRatePerTestcase)
    system_problem = TspProblem(system_case, Symbol[], [ :features ])
    system_instance = TspInstance(system_problem, "System 1", 1)
    save(system_instance, exp_dir)
end

function run_experiment_series_1_system1_as_objective(exp_id::String)
    _cat_objectives = Symbol[ :features ]
    series = TspExperimentalSeries(_algs_setup_objective, _runs, _runtime, _cat_objectives, NO_CRITERIA, exp_id)
    run(series)
end

function run_experiment_series_1_system1_as_constraint(exp_id::String)
    _cat_constraints = Symbol[ :features ]
    series = TspExperimentalSeries(_algs_setup_constraint, _runs, _runtime, NO_CRITERIA, _cat_constraints, exp_id)
    run(series)
end

function run_experiment_series_1_system2_as_objective(exp_id::String)
    _cat_objectives = Symbol[ :callchaindepth3, :morphologies ]
    series = TspExperimentalSeries(_algs_setup_objective, _runs, _runtime, _cat_objectives, NO_CRITERIA, exp_id)
    run(series)
end

function run_experiment_series_1_system2_as_constraint(exp_id::String)
    _cat_constraints = Symbol[ :callchaindepth3, :morphologies ]
    series = TspExperimentalSeries(_algs_setup_constraint, _runs, _runtime, NO_CRITERIA, _cat_constraints, exp_id)
    run(series)
end

function run_experiment_series_1_system_1()
    exp_id = "original_system_1_objective"
    setup_experiment_for_huawei_system_1(exp_id)
    run_experiment_series_1_system1_as_objective(exp_id)

    exp_id = "original_system_1_constraint"
    setup_experiment_for_huawei_system_1(exp_id)
    run_experiment_series_1_system1_as_constraint(exp_id)
end

function run_experiment_series_1_system_2()
    exp_id = "original_system_2_objective"
    setup_experiment_for_huawei_system_2(exp_id)
    run_experiment_series_1_system2_as_objective(exp_id)

    exp_id = "original_system_2_constraint"
    setup_experiment_for_huawei_system_2(exp_id)
    run_experiment_series_1_system2_as_constraint(exp_id)
end

function run_experiment_series_1()
    run_experiment_series_1_system_1()
    run_experiment_series_1_system_2()
end