println("load sparse libs!")

function BlackBoxOptim.acquire_candi_fresh(pop::FitPopulation{F}) where {F}
    ones = randsubseq(1:numdims(pop), .3)
    fresh = sparsevec(ones, 0x01, numdims(pop))
    return BlackBoxOptim.Candidate{F}(fresh, -1, pop.nafitness)
end

# optimized mutation version, sparse (with 2 helper functions)
function sparsevec_float_binary(original::SparseVector, added_idx::Vector{<:Integer}, removed_idx::Vector{<:Integer})
    new_idxs = copy(original.nzind)
    setdiff!(new_idxs, removed_idx)
    new_idxs = vcat(new_idxs, added_idx)

    return sparsevec(new_idxs, 0x01, length(original))
end

function indices_swap(values::AbstractVector{<:Real})
    _ones = findall(x -> x > 0, values)
    _zeros = findall(x -> x < 1, values)

    return _ones, _zeros
end