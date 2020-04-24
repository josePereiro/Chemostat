"""
    A dict (json, mat) base data structure to store results
    related with the Chemostat
"""

import JSON: json, parse
import MetabolicEP: MetNet
import SparseArrays: sparse, sparsevec, findnz, SparseVector

include("metnet_json.jl")