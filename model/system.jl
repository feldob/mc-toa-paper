using JSON

function import_systems_json(file_name::String)
    local d::AbstractDict

    open(file_name, "r") do file
        s = read(file, String)
        d = JSON.parse(s)
    end
    return d
end

struct System
    systemid::String
    testCaseCountPerFunction::AbstractVector{<:UInt}
    functionsCoveredPerTestCase::AbstractVector{<:UInt}

    function System(d::Dict, systemid::String)
        sysdict = d[systemid]

        tcc = convert.(UInt64, sysdict["num_testcases_called_per_func"])
        fctc = convert.(UInt64, sysdict["num_funcs_called_per_testcase"])
        return new(systemid, tcc, fctc)
    end

    System(systemid, tcc, fctc) = new(systemid, tcc, fctc)
end

systemid(s::System) = s.systemid
numtests(s::System) = length(s.functionsCoveredPerTestCase)
numfuncs(s::System) = length(s.testCaseCountPerFunction)
empiricalTestCaseCount(s::System) = s.testCaseCountPerFunction
sparsity(s::System) = 1 - (sum(s.testCaseCountPerFunction) / (numtests(s)*numfuncs(s)))
