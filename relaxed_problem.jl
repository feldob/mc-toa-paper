include("io.jl")

using JuMP, Clp

lp = false

exp = load_instance("System 1_1_1.0", [:features], Symbol[], 2)

if lp
    m = Model(Clp.Optimizer)
    vars = @variable(m, 0 ≤ x[1:numtests(problem(exp))] ≤ 1)
else
    m = Model(Gurobi.Optimizer)
    vars = @variable(m, x[1:numtests(problem(exp))], Bin)
end



# constraints
_case = case(exp)

for (c, df) in (pairs ∘ categories)(_case)
    for p in names(df)
        svn = slack_var_name(c, p)
        if lp
            v = @variable(m, base_name = svn)
            set_lower_bound(v, 0.0)
            set_upper_bound(v, 1.0)
        else
            @variable(m, base_name = svn, binary = true)
        end
    end
end

cid = :features
category_df = categories(_case)[cid]
slack_var_ids = map(p -> slack_var_name(cid, p), names(category_df))
slack_vars = map(v_id -> variable_by_name(m, v_id), slack_var_ids)

for idx in eachindex(slack_vars)
    p = category_df[:,idx]
    slack_var = slack_vars[idx]
    @constraint(m, sum(p[i] * x[i] for i in eachindex(p)) - sum(p) * slack_var ≤ 0)
end

#function must be activated
fs = map(f -> activationmatrix(_case)[f,:], 1:numfuncs(_case))
foreach(f -> @constraint(m,  sum(f[i] * x[i] for i in eachindex(f)) ≥ 1), fs)

# objective
@objective(m, Min, sum(x))

JuMP.optimize!(m)

println(termination_status(m))
vvv = value.(x)
res = Vector{Float64}(undef, length(vvv))
foreach(i -> res[i] = vvv[i], eachindex(res))
res
