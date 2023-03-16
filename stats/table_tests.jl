include("stats_table.jl")

_exp_stats = full_stats("original_system_2_constraint")
_exp_stats_norm = add_normalized_stats(_exp_stats)

_table_buffer = latex_table(_exp_stats_norm, "System 2_1_1.0")

open("example.txt","w") do f
    write(f, String(take!(_table_buffer)))
end
