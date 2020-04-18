# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

import MetabolicEP: MetNet
import SparseArrays: SparseMatrixCSC

function linear_S(M = 3)
    N = M + 1
    # S
    S = zeros(M,N)
    for i in 1:M
        S[i,i] = 1.0
        S[i,i+1] = -1.0
    end
    return SparseMatrixCSC{Float64,Int64}(S)
end

function simple_linear_MetNet(M)
    S = linear_S(M)
    M,N = size(S)
    return MetNet(S,zeros(Float64, M), zeros(Float64,N), ones(Float64, N))
end

