#!/bin/sh

# install Julia
# -------------------------------------------------------------------

JULIA_VERSION_MAJOR=1.5
JULIA_VERSION=${JULIA_VERSION_MAJOR}.0-rc1
JULIA_TAR=julia-${JULIA_VERSION}-linux-x86_64.tar.gz
JULIA_HOME=/usr/julia

mkdir ${JULIA_HOME}
wget -q https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION_MAJOR}/${JULIA_TAR}
tar -C ${JULIA_HOME} -xzf ${JULIA_TAR}

# put host Julia on global path
ln -s ${JULIA_HOME}/julia-${JULIA_VERSION}/bin/julia /usr/local/bin/julia
rm -f ${JULIA_TAR}

julia -e 'using Pkg; Pkg.add("JSON")'
julia -e 'using Pkg; Pkg.add("DataFrames")'
julia -e 'using Pkg; Pkg.add("DelimitedFiles")'
julia -e 'using Pkg; Pkg.add("CSV")'
julia -e 'using Pkg; Pkg.add("JuMP")'
julia -e 'using Pkg; Pkg.add("Cbc")'
julia -e 'using Pkg; Pkg.add("Clp")'
julia -e 'using Pkg; Pkg.add("BenchmarkTools")'
julia -e 'using Pkg; Pkg.add("Gurobi"); Pkg.build("Gurobi")'
julia -e 'using Pkg; Pkg.add("Plots")'
julia -e 'using Pkg; Pkg.add("SparseArrays")'
julia -e 'using Pkg; Pkg.add("Distributions")'
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/feldob/BlackBoxOptim.jl"))'
julia -e 'using Pkg; Pkg.add("LinearAlgebra")'
julia -e 'using Pkg; Pkg.add("CodecZlib")'
julia -e 'using Pkg; Pkg.add("Lazy")'
julia -e 'using Pkg; Pkg.add("StatsPlots")'
julia -e 'using Pkg; Pkg.add("GR")'

# install Gurobi
# -------------------------------------------------------------------

GUROBI_VERSION="9.1.0"
GUROBI_TAR=gurobi${GUROBI_VERSION}_linux64.tar.gz
GUROBI_HOME_DIR=/opt/gurobi910
GUROBI_HOME=${GUROBI_HOME_DIR}/linux64

wget -q https://packages.gurobi.com/9.1/${GUROBI_TAR}
tar -C /opt/ -xzf ${GUROBI_TAR}
rm -f ${GUROBI_TAR}

export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

# activate license (uncomment line below and run; store the file in default, i.e. users home directory)
#grbgetkey YOURVALIDKEYGOESHERE

# verify gurobi license (uncomment line below)
gurobi_cl

# execute experiment
# -------------------------------------------------------------------
