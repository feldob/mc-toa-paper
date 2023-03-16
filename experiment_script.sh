#!/bin/sh

JULIA_COPY_STACKS=1 JULIA_NUM_THREADS=$(nproc) GKSwstype=100 julia experiment_series_1.jl
JULIA_COPY_STACKS=1 JULIA_NUM_THREADS=$(nproc) GKSwstype=100 julia experiments_results_1.jl
