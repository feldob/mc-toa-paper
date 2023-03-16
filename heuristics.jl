include("model/case.jl")
include("model/fitness.jl")

# each column is a solution
struct TspParetoFront
    front::BitMatrix
end

front(f::TspParetoFront) = f.front

struct TspSolution
    solution_vector::BitVector
end

vector(s::TspSolution) = s.solution_vector

function random_full_coverage_subset(c::AbstractCase)::Vector{<:Integer}
    shortests_full_coverage_subset_in_order(randperm(numtests(c)), c)
end

function shortests_full_coverage_subset_in_order(seq::Vector{<:Integer}, c::AbstractCase)::Vector{<:Integer}
    avs = activationvectors(c)
    coverage = zeros(Bool, numfuncs(c))

    idx = 0
    while sum(coverage) < length(coverage)
        idx += 1
        covered = avs[seq[idx]].nzind
        coverage[covered] .= 1
    end

    return seq[1:idx]
end

function greedy_full_coverage_subset_fast(c::AbstractCase)::Vector{<:Integer}
    coverings = funcspertest(c)
    indx_tupes = [ (i, coverings[i]) for i in 1:length(coverings) ]
    sorted_tuples = sort(indx_tupes, rev = true, by = x -> x[2])
    order = [ s[1] for s in sorted_tuples ]

    shortests_full_coverage_subset_in_order(order, c)
end

function heuristic_choice(cov::Vector{<:Integer})
    ordered_cov = sortperm(cov, rev=true)
    lasttoconsider = count(f-> f != 0, cov)
    range = 1:round(Int64, lasttoconsider/2, RoundUp)
    return ordered_cov[rand(range)]
end

function greedy_full_coverage_subset_slow_opt(c::AbstractCase, randomized::Bool = false)::Vector{<:Integer}

    subset = Set{Int64}()

    funccoverage = Vector{Vector{Int64}}(undef, numtests(c))

    testrelevance = Vector{Vector{Int64}}(undef, numfuncs(c))
    foreach(idx -> testrelevance[idx] = Vector{Int64}(), eachindex(testrelevance))

    @inbounds for (idx, v) in enumerate(activationvectors(c))
        t_functions = v.nzind
        funccoverage[idx] = Vector{Int64}(t_functions)
        for t in t_functions
            push!(testrelevance[t], idx)
        end
    end

    remaining_funcs = numfuncs(c)
    @inbounds while remaining_funcs > 0
        current_coverage = [ length(funccoverage[t]) for t in 1:numtests(c) ]
        next_t = randomized ? heuristic_choice(current_coverage) : argmax(current_coverage)
        next_funcs_covered = copy(funccoverage[next_t])
        remaining_funcs -= length(next_funcs_covered)
        to_update = unique!(reduce(vcat, [ testrelevance[f] for f in next_funcs_covered ] ))
        foreach(t -> setdiff!(funccoverage[t], next_funcs_covered), to_update)
        push!(subset, next_t)
    end

    return collect(subset)
end

function greedy_full_coverage_subset_slow(c::AbstractCase)::Vector{<:Integer}
    subset = Set{Int64}()

    funccoverage = Vector{Set{Int64}}(undef, numtests(c))
    testrelevance = Vector{Vector{Int64}}(undef, numfuncs(c))
    foreach(idx -> testrelevance[idx] = Vector{Int64}(), eachindex(testrelevance))

    @inbounds for (idx, v) in enumerate(activationvectors(c))
        t_functions = v.nzind
        funccoverage[idx] = Set{Int64}(t_functions)
        for t in t_functions
            push!(testrelevance[t], idx)
        end
    end

    remaining_funcs = numfuncs(c)
    @inbounds while remaining_funcs > 0
        current_coverage = [ length(funccoverage[t]) for t in 1:numtests(c) ]
        next_t = argmax(current_coverage)
        next_funcs_covered = copy(funccoverage[next_t])
        remaining_funcs -= length(next_funcs_covered)
        to_update = reduce(union!, [ testrelevance[f] for f in next_funcs_covered ] )
        foreach(t -> setdiff!(funccoverage[t], next_funcs_covered), to_update)
        push!(subset, next_t)
    end

    return collect(subset)
end

function to_solution(testsubset::AbstractVector{<:Integer}, totalnumtests::Integer)::SparseVector
    return sparsevec(testsubset, UInt8(1), totalnumtests)
end
