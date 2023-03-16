include("../normalized_stats.jl")

using CSV, DataFrames, Printf

# For series 1, system 2, we have to add the "morphologies" objective too. add the optional ones as strings to pass here!? naming?
function latex_table(_exp_data::DataFrame, systemid::AbstractString)::IOBuffer
    iobuffer = IOBuffer()
    for method in unique(_exp_data.algid)
        _exp_data_method = filter(row -> row.algid == method, _exp_data)

        for (i, scale) in enumerate(unique(_exp_data.scale))
            method_data = filter(row -> row.scale == scale, _exp_data_method)
            if i == 1
                print(iobuffer, replace("$method", "_" => "\\_"))
            end

            ms_mean = round(Int64,mean(method_data.minsize), RoundUp)
            ms_std = round(Int64,std(method_data.minsize), RoundUp)

            print(iobuffer, "&\$ $ms_mean \\pm $ms_std \$")

            ms_mean_gap = round(Int64, mean(method_data.abs_size_reduction), RoundUp)
            ms_std_gap = round(Int64, std(method_data.abs_size_reduction), RoundUp)

            print(iobuffer, "&\$ $ms_mean_gap \\pm $ms_std_gap \$")

            #fr_mean = @sprintf("%.2e",mean(method_data.max_fr_per_size))
            #fr_std = @sprintf("%.2e",std(method_data.max_fr_per_size))
            fr_mean = round(mean(method_data.max_fr_per_size), digits = 1)
            fr_std = std(method_data.max_fr_per_size)

            print(iobuffer, "&\$ $fr_mean\$")

            fr_mean_gap = round(Int64, mean(method_data.dist_to_best_fr_per_size) * 100, RoundUp)
            fr_std_gap = round(Int64, std(method_data.dist_to_best_fr_per_size) * 100, RoundUp)

            print(iobuffer, "&\$ $fr_mean_gap \\pm $fr_std_gap \$")

            ex_mean = round(Int64,mean(method_data.max_ex_per_size), RoundUp)
            ex_std = round(Int64,std(method_data.max_ex_per_size), RoundUp)

            print(iobuffer, "&\$ $ex_mean \\pm $ex_std \$")

            ex_mean_gap = round(Int64, mean(method_data.dist_to_best_ex_per_size) * 100, RoundUp)
            ex_std_gap = round(Int64, std(method_data.dist_to_best_ex_per_size) * 100, RoundUp)

            print(iobuffer, "&\$ $ex_mean_gap \\pm $ex_std_gap \$")

            evals_mean_gap = round(Int64, mean(method_data.evaluations), RoundUp)
            evals_std_gap = round(Int64, std(method_data.evaluations), RoundUp)

            if method ∈ HEURISTIC_ALG_NAMES
                print(iobuffer, "&---")
            else
                print(iobuffer, "&\$ $evals_mean_gap \\pm $evals_std_gap \$")
            end

            println(iobuffer, "\\\\")
        end
        println(iobuffer, "\\hline")
    end
    return iobuffer
end

function latex_table_study2(_exp_data::DataFrame, _table_file::String)
    _table_buffer = latex_table_study2(_exp_data)
    open(_table_file,"w") do f
        write(f, String(take!(_table_buffer)))
    end
end

function latex_table(_exp_data::DataFrame, systemid::AbstractString, to_file::String)
    _table_buffer = latex_table(_exp_data, systemid)
    open(to_file,"w") do f
        write(f, String(take!(_table_buffer)))
    end
end
