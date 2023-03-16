using StatsPlots, Plots, CSV, DataFrames

_st = DataFrame!(CSV.File("experiment_original_system_1_objective/results/stats/stats.csv"))

_sts = DataFrames.groupby(_st, :algid)

_p = violin()
prop = "max_ex_per_size"
#prop = "max_fr_per_size"
#prop = "minsize"
for a in keys(_sts)
    algid = a[1]
    boxplot!([algid], _sts[a][!, prop], label = algid, leg = false)
end
xlabel!(prop)

_p
