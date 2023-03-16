# testactive vector kept generic to allow for JuMP to pass VariableRef instance

#OBS this sum only works for this problem because the number of entries in array same as sum
function tsp_size(testsactive::SparseVector)::Integer
    length(testsactive.nzval)
end

function tsp_size(testsactive::AbstractVector)
    sum(testsactive)
end

function tsp_size(testsactive::AbstractVector{<:Number})::Integer
    sum(testsactive)
end

# function tsp_exercises(testsactive::SparseVector, funcspertest::AbstractVector{<:Integer})::AbstractFloat
#     return mean(view(funcspertest, testsactive.nzind))
# end

function tsp_exercises(testsactive::AbstractVector, funcspertest::AbstractVector{<:Integer})
    return sum(testsactive[i] * funcspertest[i] for i in eachindex(testsactive))
end

# function tsp_failrate(testsactive::SparseVector, failrates::AbstractVector{<:AbstractFloat})::AbstractFloat
#     return mean(view(failrates, testsactive.nzind))
# end

function tsp_failrate(testsactive::AbstractVector, failrates::AbstractVector{<:AbstractFloat})
    return sum(testsactive[i] * failrates[i] for i in eachindex(testsactive))
end

function tsp_coverage(test_indices::AbstractVector{<:Integer}, activationVectors::AbstractVector{SparseVector{UInt8}}, nfunctions::Integer)::AbstractFloat
    if length(test_indices) == 0
        return 0.0
    end

    relevant_test_vectors = view(activationVectors, test_indices)
    pos = zeros(Bool, nfunctions)
    @inbounds foreach( v -> pos[v.nzind] .= 1, relevant_test_vectors)
    unique_functions = sum(pos)

    return unique_functions / nfunctions
end

function tsp_coverage(test_indices::AbstractVector{<:Integer}, activationMatrix::AbstractMatrix{UInt8}, nfunctions::Integer)::AbstractFloat
    _m = view(activationMatrix, :, test_indices)
    _s = sum(_m, dims = 2)[:,1]
    return (nfunctions - length(_s[_s .< 1])) / nfunctions
end

function tsp_coverage(testsactive::SparseVector, activationVectors::AbstractVector{SparseVector{UInt8}}, nfunctions::Integer)::AbstractFloat
    return tsp_coverage(testsactive.nzind, activationVectors, nfunctions)
end

function tsp_coverage(testsactive::AbstractVector{UInt8}, activationVectors::AbstractVector{SparseVector{UInt8}}, nfunctions::Integer)::AbstractFloat
    return tsp_coverage(findall(x -> x == 1, testsactive), activationVectors, nfunctions)
end

using DataFrames

# OBS per contract, all vectors are binary vectors of same length, but actual datatype may vary.
function tsp_criterion_coverage(testsactive::AbstractVector{T}, properties_df::DataFrame) where {T <: Integer}
    activation_per_property = map(n -> properties_df[:,n] .* testsactive, names(properties_df))
    cov_per_property = map(av -> sum(av) > 0, activation_per_property)
    return sum(cov_per_property ./ ncol(properties_df))
end
