# general functions
# indexing using ξ, data_key
function _apply_out_getter(getter::Function, 
        bundle::ChstatBundle, ξ::Real, data_key::Symbol)
    return getter(bundle[ξ, data_key])
end

# indexing using ξ, data_key, ider(s)
function _apply_out_getter(getter::Function, 
        bundle::ChstatBundle, ξ::Real, data_key::Symbol, 
        ider, metnet_data_key::Symbol = metnet_data_key)
    
    metnet = bundle[ξ, metnet_data_key]
    out = bundle[ξ, data_key]

    return getter(metnet, out, ider)
end

# indexing using ξs, data_key, ider
_apply_out_getter(getter::Function, 
        bundle::ChstatBundle, ξs::Vector, data_key::Symbol, 
        ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key) = 
    [_apply_out_getter(getter, bundle, ξ, data_key, ider, metnet_data_key) for ξ in ξs]

# indexing using ξ, β, data_key
function _apply_out_getter(getter::Function, 
        bundle::ChstatBundle, ξ::Real, β::Real, data_key::Symbol)
    return getter(bundle[ξ, β, data_key])
end
# indexing using ξ, β, data_key, ider(s)
function _apply_out_getter(getter::Function, bundle::ChstatBundle, 
        ξ::Real, β::Real, data_key::Symbol, 
        ider, metnet_data_key::Symbol = metnet_data_key)
    
    metnet = bundle[ξ, metnet_data_key]
    out = bundle[ξ, β, data_key]

    return getter(metnet, out, ider)
end

# indexing using ξs, β, data_key, ider
_apply_out_getter(getter::Function, bundle::ChstatBundle, 
        ξs::Vector, β::Real, data_key::Symbol, 
        ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key) = 
    [_apply_out_getter(getter, bundle, ξ, β, data_key, ider, metnet_data_key) for ξ in ξs]

# indexing using ξ, βs, data_key, ider
_apply_out_getter(getter::Function, bundle::ChstatBundle, 
        ξ::Real, βs::Vector, data_key::Symbol, 
        ider::IDER_TYPE, metnet_data_key::Symbol = metnet_data_key) = 
    [_apply_out_getter(getter, bundle, ξ, βs, data_key, ider, metnet_data_key) for βs in βs]

# eval for all Out getters
for fun in [av, va, μ, σ]
    fun_name = string(nameof(fun))

    # indexing using ξ, data_key
    fun_str = 
        """
        $(fun_name)(bundle::ChstatBundle, 
                ξ::Real, data_key::Symbol) = 
            _apply_out_getter($(fun_name), bundle, ξ, data_key)
        """
    eval(Meta.parse(fun_str))

    # indexing using ξ(s), data_key, ider(s)
    fun_str = 
        """
         $(fun_name)(bundle::ChstatBundle, 
                ξ, data_key::Symbol, ider, 
                metnet_data_key::Symbol = metnet_data_key) = 
            _apply_out_getter($(fun_name), bundle, ξ, data_key, ider, metnet_data_key)
        """
    eval(Meta.parse(fun_str))

    # indexing using ξ, β, data_key
    fun_str = 
        """
        $(fun_name)(bundle::ChstatBundle, 
                ξ::Real, β::Real, data_key::Symbol) = 
            _apply_out_getter($(fun_name), bundle, ξ, β, data_key)
        """
    eval(Meta.parse(fun_str))

    # indexing using ξ(s), β(s), data_key, ider(s)
    fun_str = 
        """
        $(fun_name)(bundle::ChstatBundle, 
                ξ, β, data_key::Symbol, 
                ider, metnet_data_key::Symbol = metnet_data_key) = 
            _apply_out_getter($(fun_name), bundle, ξ, β, data_key, ider, metnet_data_key)
        """
    eval(Meta.parse(fun_str))

end