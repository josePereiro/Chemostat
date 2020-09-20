# function metnet_to_json(metnet::MetNet)
#     data_dict = Dict()
#     data_dict["N"] = metnet.N
#     data_dict["M"] = metnet.M
#     data_dict["S"] = findnz(sparse(metnet.S))
#     data_dict["b"] = findnz(sparsevec(metnet.b))
#     data_dict["c"] = findnz(sparsevec(metnet.c))
#     data_dict["lb"] = findnz(sparsevec(metnet.lb))
#     data_dict["ub"] = metnet.ub # ub use to be non-dense
#     data_dict["genes"] = metnet.genes
#     data_dict["rxnGeneMat"] = findnz(sparse(metnet.rxnGeneMat))
#     data_dict["grRules"] = metnet.grRules
#     data_dict["mets"] = metnet.mets
#     data_dict["rxns"] = metnet.rxns
#     data_dict["metNames"] = metnet.metNames
#     data_dict["metFormulas"] = metnet.metFormulas
#     data_dict["rxnNames"] = metnet.rxnNames
#     data_dict["rev"] = "Should be reconstructed from bounds"
#     data_dict["subSystems"] = metnet.subSystems
#     return json(data_dict)
# end

# # TODO check this for improvement
# function metnet_from_json(json_str)
#     data_dict = parse(json_str)
#     N = Int(data_dict["N"])
#     M = Int(data_dict["M"])
#     S = sparse(data_dict["S"]..., M, N)
#     b = collect(SparseVector{Float64, Int}(sparsevec(data_dict["b"]..., M)))
#     c = collect(SparseVector{Float64, Int}(sparsevec(data_dict["c"]..., N)))
#     lb = collect(SparseVector{Float64, Int}(sparsevec(data_dict["lb"]..., N)))
#     ub = Vector{Float64}(data_dict["ub"])
#     genes = Vector{String}(data_dict["genes"])
#     rxnGeneMat = sparse(data_dict["rxnGeneMat"]..., M, N)
#     grRules = data_dict["grRules"]
#     mets = data_dict["mets"]
#     rxns = data_dict["rxns"]
#     metNames = data_dict["metNames"]
#     metFormulas = data_dict["metFormulas"]
#     rxnNames = data_dict["rxnNames"]
#     rev = (lb .< 0.0) .& (lb .> 0.0)
#     subSystems = data_dict["subSystems"]
    
#     return MetNet(S,b,c,lb,ub,genes,
#         rxnGeneMat,grRules,mets,rxns,metNames,
#         metFormulas,rxnNames,rev,subSystems)
# end
