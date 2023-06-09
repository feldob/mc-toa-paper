#!/bin/sh

JULIA_COPY_STACKS=1 JULIA_NUM_THREADS=$(nproc) GKSwstype=100 julia example.jl
