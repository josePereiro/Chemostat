const cost_met_id = "cost"
const cost_met_name = "Cost per flux unit"
const cost_rxn_id = "tot_cost"
const cost_rxn_name = "Total flux cost"
const cost_rxn_lb = 0.0
const cost_rxn_ub = 1.0

"""
    Add a new metabolite simulating the cost penalazing 
    reaction fluxes. 
    A new balance equations is then added:
        Σ(rᵢ*costᵢ) + tot_cost = 0
    Because the cost coefficients (costᵢ) > 0, the system must allocate 
    the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
    are usually bounded [0.0, 1.0]
"""
function add_costs(metnet::MetNet, cost_info::Dict; 
        cost_met_id::AbstractString = cost_met_id,
        cost_met_name::AbstractString = cost_met_name, 
        cost_rxn_id::AbstractString = cost_rxn_id, 
        cost_rxn_name::AbstractString = cost_rxn_name,
        cost_rxn_lb::Real = cost_rxn_lb, 
        cost_rxn_ub::Real = cost_rxn_ub)
    
    # TODO check that internal reactions should be irreversibles
    
    # For being concistence with the current implementation,
    # costs must be deffined negative
    for (rxn, cost) in cost_info
        cost > 0 && error("cost of ($rxn) must be negative, just negate it, the abs value is what matters!!!")
    end

    metnet_with_costs =  add_met(metnet, cost_met_id, 
        rxns = cost_info, metName = cost_met_name)
    return add_rxn(metnet_with_costs, cost_rxn_id,
        mets = Dict(cost_met_id => 1.0),
        lb = cost_rxn_lb, ub = cost_rxn_ub, 
        rxnName = cost_rxn_name)
end