using Distributions, SparseArrays, Random, Lazy, DataFrames

include("system.jl")

# based on experience. For some companies lower, i.e. Normal(0.03, 0.02). Not used for the study.
H_FAILURE_RATE_DISTRIBUTION = truncated(Normal(0.08, 0.05), 0.0, 1.0)

abstract type AbstractCase end

struct Case <: AbstractCase
    activation::SparseMatrixCSC{UInt8}
    failrates::Vector{Float64}
    funcspertest::AbstractVector{UInt64}
    avgexercises::Float64
    maxexercises::Float64
    avgfailrate::Float64
    maxfailrate::Float64
    categories::Dict{Symbol, DataFrame}

    function Case(activationmatrix, categories=Dict{Symbol, DataFrame}(), failurerates = rand(H_FAILURE_RATE_DISTRIBUTION, size(activationmatrix, 2)))
        funcspertest = sum(UInt64, activationmatrix, dims=1)[1,:]
        avgexercises = mean(activationmatrix)
        avgfailrate = mean(failurerates)
        maxexercises = maximum(activationmatrix)
        maxfailrate = maximum(failurerates)
        new(activationmatrix, failurerates, funcspertest, avgexercises, maxexercises, avgfailrate, maxfailrate, categories)
    end
end

Case(am::Matrix{UInt8}, c, fr) = Case(sparse(am), c, fr)
activationmatrix(c::Case) = c.activation
failrates(c::Case) = c.failrates
numtests(c::Case) = size(c.activation,2)
numfuncs(c::Case) = size(c.activation,1)
testsperfunc(c::Case) = sum(UInt64, c.activation, dims=2)[:,1]
funcspertest(c::Case) = c.funcspertest
avgexercises(c::Case) = c.avgexercises
maxexercises(c::Case) = c.maxexercises
avgfailrate(c::Case) = c.avgfailrate
maxfailrate(c::Case) = c.maxfailrate
categories(c::Case) = c.categories

struct CachedCase <: AbstractCase
    case::Case
    activation_vectors::AbstractVector{SparseVector{UInt8}}

    function CachedCase(case::Case)
        am = activationmatrix(case)
        activation_vectors = [ am[:, t]  for t in 1:numtests(case)]
        new(case, activation_vectors)
    end
end

@forward CachedCase.case activationmatrix, failrates, numtests, numfuncs, testsperfunc,
        funcspertest, avgexercises, maxexercises, avgfailrate, maxfailrate, categories

activationvectors(c::CachedCase) = c.activation_vectors

# due to the long tail and min 1, chosen Pareto distribution to be suitable.
# OBS independent of size of system, i.e. number of functions. Further, model fit might be non-optimal.
model(s::System) = fit(Pareto, empiricalTestCaseCount(s))

# Generalizes for multiple systems simplified by Pareto with α set to mean α.
# OBS model fit might be non-optimal.
generalized_model(sys::System...) = Pareto(mean(shape.(model.(sys))), 1)

# Pareto -> Inf | guarantee possible number of tests to cover function
# convert to discrete from continuous distribution
function generate_case(model, numfunctions::Integer, numtests::Integer, max_coveringtests::Integer=numtests)::Case
    covering_numbers = map(_ -> trunc(Int64, min(rand(model), max_coveringtests)), 1:numfunctions)
    test_coverings = map(cn -> randperm(numtests)[1:cn], covering_numbers)

    J = reduce(vcat, test_coverings)
    I = reduce(vcat, map((idx_tc) -> idx_tc[2][:] .= idx_tc[1], enumerate(copy(test_coverings))))
    activationmatrix = sparse(I, J, UInt8(1), numfunctions, numtests)

    return Case(activationmatrix)
end

sparsity(m::AbstractMatrix{T}) where {T <: Integer} = 1 - sum(m) / length(m)

function generate_case(model, s::System)
    nf = numfuncs(s)
    nt = numtests(s)
    e = empiricalTestCaseCount(s)

    generate_case(model, nf, nt, maximum(e))
end

function generate_case_raw(s::System)
    generate_case(empiricalTestCaseCount(s), s)
end

function generate_case_model(s::System)
    case = generate_case(model(s), s)
end

scale(s::System, nt::Integer) = scale(s, nt / numtests(s))
function scale(s::System, ratio::AbstractFloat)::System
    etcc = round.(Int64, copy(empiricalTestCaseCount(s)) .* ratio, RoundUp)
    #etcc = empiricalTestCaseCount(s) # no scaling, @Robert?

    nf = round.(Int64, numfuncs(s) * ratio, RoundUp)
    nt = round.(Int64, numtests(s) * ratio, RoundUp)

    bootstrap_case = generate_case(etcc, nf, nt, maximum(etcc))

    tcc = testsperfunc(bootstrap_case)
    fctc = funcspertest(bootstrap_case)

    return System(systemid(s), tcc, fctc)
end

function create_category(_numtests::Int64, category_names::Vector{String}, distr=fill!(Vector{Float64}(undef, length(category_names)), 1/length(category_names)))::DataFrame
    _perm = randperm(_numtests)

    previous_cut_point = 0
    df = DataFrame()
    for i in eachindex(distr)
        next_cut_point::Int64 = previous_cut_point + round(_numtests * distr[i], RoundUp)
        cat_idx = _perm[previous_cut_point+1:min(next_cut_point, _numtests)]
        df[:, category_names[i]] = sparsevec(cat_idx, 1, _numtests)
        previous_cut_point = next_cut_point
    end

    return df
end
