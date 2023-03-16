using Plots, StatsPlots

include("../model/system.jl")

d = import_systems_json("data/stats_per_func_and_testcase_sys1_and_sys2_20200624.json")

sys1 = System(d, "System1")
sys2 = System(d, "System2")

emp1 = empiricalTestCaseCount(sys1)
emp2 = empiricalTestCaseCount(sys2)

hist1 = histogram(emp1, title = "system 1 (all)")

ylabel!("covered functions")
xlabel!("test case quantities")

hist2 = histogram(emp2, title = "system 2 (all)")

ylabel!("covered functions")
xlabel!("test case quantities")

hist3 = histogram(emp1[emp1 .< 40], title = "system 1 (up to 40)")

ylabel!("covered functions")
xlabel!("test case quantities")

hist4 = histogram(emp2[emp2 .< 10], title = "system 2 (up to 10)")

ylabel!("covered functions")
xlabel!("test case quantities")

plot(hist1, hist2, hist3, hist4, layout=(2,2), legend = false)
