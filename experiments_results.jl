include("stats/stats_table.jl")
include("stats/solver_plots.jl")
include("stats/hv_indicator.jl")

STATS_FILE = "stats.csv"
TABLE_FILE = "_table.txt"
RUNTIME_FILE = "_runtimes.pdf"
EVALS_FILE = "_evals.pdf"
HVI_FILE = "_hvi.pdf"

HEURISTIC_ALG_NAMES = ["borg", "borg_greedy", "random", "random_greedy"]

# to keep order in plots
_ALG_NAMES = Dict("gurobi" => 1,
                "gurobi_relaxed" => 2,
                "Cbc" => 3,
                "Cbc_relaxed" => 4,
                "Clp" => 5,
                "borg" => 6,
                "borg_greedy" => 7,
                "random" => 8,
                "random_greedy" => 9)

_COL_NAMES = Dict(1 => :black,# :blue,
                2 => :black,#:lightblue,
                3 => :white,#:green,
                4 => :white,#:lightgreen,
                5 => :white,#:red,
                6 => :blue,#:purple,
                7 => :lightblue,#:pink,
                8 => :orange,
                9 => :red)

colorindex_for(_algid) = _COL_NAMES[_ALG_NAMES[_algid]]

function save_stats_file(statsfolder, _exp_stats::DataFrame)
    CSV.write(joinpath(statsfolder, STATS_FILE), _exp_stats)
end

function create_stats(expid::Union{String, Integer},
                        _objectives_labels,
                        _objectives_tomax,
                        _runtime,
                        _algids,
                        _exp_stats::DataFrame=full_stats(expid))
    statsfolder = joinpath("experiment_$expid", "results", "stats")

    "creating stats for project '$expid'."|> println

    if !isdir(statsfolder)
        mkdir(statsfolder)
    end

    #save_stats_file(statsfolder, _exp_stats) # reduntant from the full_stats call, right?
    _exp_stats_norm = add_normalized_stats(_exp_stats)

    "load results. (might take a while...)"|> print
    _results = load_tsp_results(expid)

    " done!"|> println

    algids = _algids ∩ unique(_results.algid)

    _caseids = map(r -> "$(r.systemid)_$(r.instanceid)_", eachrow(_results))
    for _caseid in unique(_caseids)
        systemid = string(split(_caseid, "_")[1])
        _caseid = filter(f -> startswith(f, _caseid), readdir("experiment_$(expid)"))[1]# gotta complete it
        "---- $systemid ----" |> print
        case_startof = joinpath(statsfolder, systemid)
        _table_file = "$case_startof$TABLE_FILE"
        "extract stats..."
        extract_stats(expid, _caseid)
        "create table..." |> print
        latex_table(_exp_stats_norm, systemid, _table_file)
        relev_algids = setdiff(algids, HEURISTIC_ALG_NAMES)
        relev_algids |> println
        _runtime_file = "$case_startof$RUNTIME_FILE"
        "create runtime plots..." |> print
        scale_runtime_plots(_results, systemid, _runtime_file, relev_algids)
        " done!"|> println
        _evals_file = "$case_startof$EVALS_FILE"
        "create eval plots..." |> print
        scale_eval_plots(_results, systemid, _evals_file, relev_algids)
        " done!"|> println
        "create hypervolume plots..." |> print

        expfilter = (f) -> (startswith(f, _caseid) && !reduce(|, contains.(f, HEURISTIC_ALG_NAMES)))
        _p = plot_hypervolumes(expid, _objectives_labels, _objectives_tomax, _runtime, _caseid, _results, expfilter)
        savefig(_p, "$case_startof$HVI_FILE")
        " done!"|> println
        title!(_p, case_startof)
    end

    "Stats creation finished successfully."|> println
end