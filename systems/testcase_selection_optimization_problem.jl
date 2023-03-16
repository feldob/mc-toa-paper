using JuMP, Cbc, JSON

# We create a simple struct to hold all these intermediate data structures.
# We typically use sparse matrices since they are often not very completely filled.
# But we use dense vectors since they are typically fully filled.
mutable struct TestcaseSelectionOptimizationProblem
    # Following 7 datastructures are used during the optimization
    functionByTestcase::AbstractMatrix{UInt8} # We only need to represent a bit so use UInt8
    selectedTestcases::AbstractVector{UInt8}
    daysSinceLastExecPerTestcase::AbstractVector{Float64}
    testExecutionTimePerTestcase::AbstractVector{Float64}
    testfeatureByTestcase::AbstractMatrix{UInt8}
    hwmorphologyByTestcase::AbstractMatrix{UInt8}
    failRatePerTestcase::AbstractVector{Float64}

    # Following datastructures are needed to map to/from the intermediate
    # data structures above and the outside world (test case names etc).
    func2row::Dict{Any,Int}     # map from function name to its number
    tc2col::Dict{Any,Int}       # map from testcase to its number
    feature2row::Dict{Any,Int} # testfeatures in the order of rows used above
    morphology2row::Dict{Any,Int} # hardware morphologies in the order of rows used above
    testcases::AbstractVector{Any} # testcases in the order they are used in all datastructures

    # Save timing data for later reporting
    timings::Dict{String, Float64}
end

using SparseArrays

# Typically we load the data from two json files, one mapping testcases to
# the functions, the other mapping test cases to info about them.
function TestcaseSelectionOptimizationProblem(tc2funcsfile::String,
    tcinfofile::String; excludeifnochain = true)
    timings = Dict{String, Float64}()

    t0 = time()
    tc2funcs, tcinfo = loadfromfiles(tc2funcsfile, tcinfofile)
    timings["load_data_from_files"] = time() - t0

    t0 = time()
    func2row, tc2col, feature2row, morphology2row, testcases =
        unpackdata(tc2funcs, tcinfo; excludeifnochain = true)
    timings["unpack_data"] = time() - t0

    t0 = time()
    optdata = make_internal_opt_datastructures(tc2funcs, tcinfo,
        func2row, tc2col, feature2row, morphology2row, testcases)
    timings["make_internal_structs"] = time() - t0
    timings["preprocess_data"] = timings["unpack_data"] +
        timings["make_internal_structs"]

    TestcaseSelectionOptimizationProblem(optdata...,
        func2row, tc2col, feature2row, morphology2row, testcases, timings)
end

ntestcases(p::TestcaseSelectionOptimizationProblem) = length(p.testcases)
nfunctions(p::TestcaseSelectionOptimizationProblem) = size(p.functionByTestcase, 1)
ntestfeatures(p::TestcaseSelectionOptimizationProblem) = size(p.testfeatureByTestcase, 1)
nmorphologies(p::TestcaseSelectionOptimizationProblem) = size(p.hwmorphologyByTestcase, 1)

const TSOP = TestcaseSelectionOptimizationProblem

totaldayssincelastexecuted(p::TSOP, sel::AbstractVector{Int}) =
    sum(p.daysSinceLastExecPerTestcase[sel])

nfunccalls(p::TSOP, sel::AbstractVector{Int}) = sum(callsperfunc(p, sel))

callsperfunc(p::TSOP, sel::AbstractVector{Int}) =
    p.functionByTestcase * sparsevec(sel, 1, ntestcases(p))

coveredfunctions(p::TSOP, sel::AbstractVector{Int}) =
    callsperfunc(p, sel) .>= 1

changedfunccoverage(p::TSOP, sel::AbstractVector{Int}) =
    sum(coveredfunctions(p, sel)) / nfunctions(p)

coverageperfeature(p::TSOP, sel::AbstractVector{Int}) =
    p.testfeatureByTestcase * sparsevec(sel, 1, ntestcases(p))

coveredfeatures(p::TSOP, sel::AbstractVector{Int}) =
    coverageperfeature(p, sel) .>= 1

featurecoverage(p::TSOP, sel::AbstractVector{Int}) =
    sum(coveredfeatures(p, sel)) / ntestfeatures(p)

coveragepermorphology(p::TSOP, sel::AbstractVector{Int}) =
    p.hwmorphologyByTestcase * sparsevec(sel, 1, ntestcases(p))

coveredmorphologies(p::TSOP, sel::AbstractVector{Int}) =
    coveragepermorphology(p, sel) .>= 1

morphologycoverage(p::TSOP, sel::AbstractVector{Int}) =
    sum(coveredmorphologies(p, sel)) / nmorphologies(p)

testexecutiontime(p::TSOP, sel::AbstractVector{Int}) =
        sum(p.testExecutionTimePerTestcase[sel])

totalfailrate(p::TSOP, sel::AbstractVector{Int}) =
        sum(p.failRatePerTestcase[sel])

avgfailrate(p::TSOP, sel::AbstractVector{Int}) =
        totalfailrate(p, sel) / length(sel)

r(v) = round(v, digits=2)
aspct(v) = r(100.0 * v)
aspct(n, N) = aspct(n/N)

function objectives(p::TSOP, sel::AbstractVector{Int})
    Dict("O1_ChangedFuncCoverage" => aspct(changedfunccoverage(p, sel)),
        "O2_SubsetSize" => length(sel),
        "O3_ChangedFuncCalls" => nfunccalls(p, sel),
        "O4_DaysSinceLastTestExec" => totaldayssincelastexecuted(p, sel),
        "O5_TestExecutionTime" => testexecutiontime(p, sel),
        "O6_TestFeatureCoverage" => aspct(featurecoverage(p, sel)),
        "O7_MorphologyCoverage" => aspct(morphologycoverage(p, sel)),
        "O8_AvgFailRate" => aspct(avgfailrate(p, sel)),
    )
end

function pctdiff(v2, v1)
    isnan(v1) && isnan(v2) && return ""
    p = aspct((v2-v1)/v1)
    p == 0.0 && return(" 0.0%")
    (p > 0.0) ? ("+" * string(p) * "%") : (string(p) * "%")
end
rightfillto(s::String, len) = (length(s) < len) ? (s * (" "^(len-length(s)))) : s
rightfillto(o, len) = rightfillto(string(o), len)

function make_internal_opt_datastructures(tc2funcs, tcinfo,
        func2row, tc2col, feature2row, morphology2row, testcases)
    nfuncs = length(func2row)
    ntcs = length(testcases)
    @assert ntcs == length(tc2col)
    nfeatures = length(feature2row)
    nmorphologies = length(morphology2row)

    # Create the datastructures first and then we fill them in
    selectedTestcases = sparsevec(Int[], UInt8[], ntcs)
    daysSinceLastExecPerTestcase = zeros(Float64, ntcs)
    testExecutionTimePerTestcase  = zeros(Float64, ntcs)
    failRatePerTestcase = zeros(Float64, ntcs)

    I, J = Int[], Int[]
    for (tc, funcs) in tc2funcs
        for func in funcs
            push!(I, func2row[func])
            push!(J, tc2col[tc])
        end
    end
    functionByTestcase = sparse(I, J, UInt8(1), nfuncs, ntcs)

    tfI, tfJ, mI, mJ = Int[], Int[], Int[], Int[]
    for (tc, info) in tcinfo
        if !haskey(tc2col, tc)
            continue # This was a testcase with no function call chain and it has not been added to tc2col during unpackdata so we skip it also here...
        end

        tci = tc2col[tc]

        daysSinceLastExecPerTestcase[tci] =
            valueforfirstkey(info, AltDaysSinceLastExecKeys, 1.0)
        @assert daysSinceLastExecPerTestcase[tci] >= 0.0

        testExecutionTimePerTestcase[tci] =
            valueforfirstkey(info, AltTestExecTimeKeys, 10.0)
        @assert testExecutionTimePerTestcase[tci] >= 0.0

        failRatePerTestcase[tci] =
            valueforfirstkey(info, AltFailRateKeys, 10.0)
        if failRatePerTestcase[tci] > 1.0
            if 0.0 <= failRatePerTestcase[tci] <= 100.0
                @info "Normalizing fail rate to 0.0-1.0 range"
                failRatePerTestcase[tci] = failRatePerTestcase[tci]/100.0
            else
                error("A failrate outside of the range 0.0-1.0 or 0.0-100.0 was found $(failRatePerTestcase[tci])")
            end
        end

        for feature in getelements(info, AltTestFeatureKeys)
            push!(tfI, feature2row[feature])
            push!(tfJ, tci)
        end

        for morphology in getelements(info, AltMorphologiesKeys)
            push!(mI, morphology2row[morphology])
            push!(mJ, tci)
        end
    end
    testfeatureByTestcase = sparse(tfI, tfJ, UInt8(1), nfeatures, ntcs)
    hwmorphologyByTestcase = sparse(mI, mJ, UInt8(1), nmorphologies, ntcs)

    return functionByTestcase, selectedTestcases, daysSinceLastExecPerTestcase,
        testExecutionTimePerTestcase, testfeatureByTestcase, hwmorphologyByTestcase,
        failRatePerTestcase
end

function loadfromfiles(tc2funcsfile::String, tcinfofile::String)
    t0 = time()
    tc2funcs = open(tc2funcsfile, "r") do fh
        JSON.parse(read(fh, String))
    end
    tcinfo = open(tcinfofile, "r") do fh
        JSON.parse(read(fh, String))
    end
    t = time() - t0
    return tc2funcs, tcinfo, t
end

const AltTestFeatureKeys = [
    "testcasefeatures", "testfeatures", "testcasefeature", "testfeature"
]

const AltMorphologiesKeys = [
    "hwmorphologies", "morphologies", "hwmorphology", "morphology",
]

const AltFailRateKeys = [
    "failrate", "failurerate", "failureprobability"
]

const AltTestExecTimeKeys = [
    "executiontime", "testexecutiontime"
]


const AltDaysSinceLastExecKeys = [
    "dayssincelastexecution", "dayssincelasttestexecution"
]

function unpackdata(tc2funcs::Dict, tcinfo::Dict; excludeifnochain = true)
    func2row = Dict{Any, Int}()
    tc2col = Dict{Any, Int}()
    testcases = Any[]
    lastusedfuncidx = 0
    for (tc, funcs) in tc2funcs
        if !haskey(tc2col, tc)
            push!(testcases, tc)
            tc2col[tc] = length(testcases)
        end
        for func in funcs
            !haskey(func2row, func) && (func2row[func] = (lastusedfuncidx += 1))
        end
    end

    feature2row = Dict{Any,Int}()
    morphology2row = Dict{Any,Int}()
    lastusedmorphologyidx = lastusedfeatureidx = 0
    for (tc, info) in tcinfo
        if excludeifnochain && !haskey(tc2col, tc)
            continue
        end

        if !haskey(tc2col, tc)
            push!(testcases, tc)
            tc2col[tc] = length(testcases)
        end

        for feature in getelements(info, AltTestFeatureKeys)
            if !haskey(feature2row, feature)
                feature2row[feature] = (lastusedfeatureidx += 1)
            end
        end

        for morphology in getelements(info, AltMorphologiesKeys)
            if !haskey(morphology2row, morphology)
                morphology2row[morphology] = (lastusedmorphologyidx += 1)
            end
        end
    end

    return func2row, tc2col, feature2row, morphology2row, testcases
end

function valueforfirstkey(d::Dict{K,V}, keys::AbstractVector{K}, default) where {K,V}
    k = firstkey(d, keys)
    if isnothing(k) && !in(nothing, keys)
        return default
    else
        return d[k]
    end
end

function firstkey(d::Dict{K,V}, keys::AbstractVector{K}) where {K,V}
    for k in keys
        haskey(d, k) && return(k)
    end
    return nothing
end

function getelements(info::Dict, altkeys)
    k = firstkey(info, altkeys)
    return isnothing(k) ? [] : parsecommaseparated(info[k])
end

parsecommaseparated(f::AbstractVector) = f
parsecommaseparated(f::AbstractVector{S}) where {S<:AbstractString} =
    map(s -> strip(string(s)), f)
parsecommaseparated(f::AbstractString) =
    map(s -> strip(string(s)), split(string(f), r"\s*,\s*"))
