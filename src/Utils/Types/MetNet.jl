# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

fake_metsid(M) = ["M$i" for i in 1:M]
fake_rxnsid(N) = ["r$i" for i in 1:N]

struct MetNet{T<:Real}

    # Fields important for modeling
    S::AbstractMatrix{T} # Stoichiometric matrix M x N sparse
    b::AbstractVector{T} # right hand side of equation  S ν = b 
    lb::AbstractVector{T} # fluxes lower bound N elements vector
    ub::AbstractVector{T} # fluxes upper bound N elements vector
    c::AbstractVector{T} # reaction index of biomass
    rxns::Vector{String} # reactions short-name N elements
    mets::Vector{String} # metabolites short-name M elements
    intake_info::Dict # Information required for enforcing the Chemostat bounds
    
    # Mustly for compatibility with COBRA
    metNames # metabolites long-names M elements
    rxnNames # reactions long-names N elements
    metFormulas # metabolites formula M elements
    genes # gene names 
    rxnGeneMat # 
    grRules # gene-reaction rule N elements vector of strings (and / or allowed)
    rev # reversibility of reactions N elements
    subSystems # cellular component of fluxes N elements

    
    function MetNet{T}(;kwargs...) where {T<:Real}
        kwargs = Dict(kwargs)

        # Just for better errors printing
        function checkget(k, T2) 
            !haskey(kwargs, k) && error("Mandatoty field '$k' not found!!!")
            dat = kwargs[k]
            !(dat isa T2) && error("'$k' is a '$(typeof(dat))', expected '$T2'")
            return dat
        end

        # LP fields
        S = checkget(:S, AbstractMatrix{T})
        M, N = size(S)
        b = checkget(:b, AbstractVector{T})
        lb = checkget(:lb, AbstractVector{T})
        ub = checkget(:ub, AbstractVector{T})
        c = checkget(:c, AbstractVector{T})
        rxns = get(kwargs, :rxns, fake_rxnsid(N))
        mets = get(kwargs, :mets, fake_metsid(N))
        intake_info = get(kwargs, :intake_info, Dict())

        # Others
        metNames = get(kwargs, :metNames, [])
        rxnNames = get(kwargs, :rxnNames, [])
        metFormulas = get(kwargs, :metFormulas, [])
        genes = get(kwargs, :genes, [])
        rxnGeneMat = get(kwargs, :rxnGeneMat, [])
        grRules = get(kwargs, :grRules, [])
        rev = get(kwargs, :rev, [])
        subSystems = get(kwargs, :subSystems, [])

        new{T}(S, b, lb, ub, c, rxns, mets, intake_info, 
            metNames, rxnNames, metFormulas, genes, 
            rxnGeneMat, grRules, rev, subSystems)

    end
end

# Helpers
MetNet(;kwargs...) = MetNet{Float64}(;kwargs...)

# For COBRA .mat files use reshape
function MetNet(model_dict::Dict; reshape = false) 
    net = to_symbol_dict(model_dict)
    reshape && (net = reshape_mat_dict(net))
    return MetNet(;net...)
end

"""
    Create a new MetNet from a template but overwriting the fields
    of the template with the given as kwargs
"""
function MetNet(metnet::MetNet; reshape = false, net...)
    net = Dict(net)

    metnet_dict = Dict()
    for field in fieldnames(typeof(metnet))
        metnet_dict[field] = getfield(metnet, field)
    end
    
    for (k, v) in metnet_dict
        if haskey(net, k)
            metnet_dict[k] = net[k]
        end
    end
    return MetNet(metnet_dict; reshape = reshape)
end

# For compatibility with cobra mat models
function reshape_mat_dict(mat_dict::Dict; S::DataType = String, 
    F::DataType = Float64, N::DataType = Int, I::DataType = Int64)

    new_dict::Dict = Dict()

    function _reshape!(rshf::Function, k::Symbol, T::Type)
        !haskey(mat_dict, k) && return
        dat = mat_dict[k]
        new_dict[k] = dat isa T ? deepcopy(dat) : rshf(dat)
    end

    # ----------------- Typed -----------------
    # String vectors
    for k in [:comps, :metNames, :metFormulas, 
            :rxnFrom, :rxnNames, :genes, :inchis, 
            :grRules, :mets, :rxns]
        _reshape!(k, Vector{S}) do dat
            dat |> vec .|> S
        end
    end

    # Float vectors
    for k in [:b, :lb, :ub, :c]
        _reshape!(k, Vector{F}) do dat
            dat |> vec .|> F
        end
    end

    # Integer vectors
    for k in [:metComps]
        _reshape!(k, Vector{I}) do dat
            floor.(I, dat |> vec)
        end
    end

    # S
    _reshape!(:S, Matrix{F}) do dat
        Matrix{F}(dat)
    end

    # ----------------- Type free -----------------
    # Matrices
    for k in [:rxnGeneMat]
        _reshape!(k, Matrix) do dat
            Matrix(dat)
        end
    end

    # Vectors
    for k in [:subSystems, :rev]
        _reshape!(k, Matrix) do dat
            dat |> vec |> collect
        end
    end

    return new_dict
end