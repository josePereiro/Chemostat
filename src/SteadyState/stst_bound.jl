"""

The chemostat steady state impose a contraint over the exchanges:

``intake <= c/ξ``
    
where ``c`` is the concentration of the metabolite in the medium and ``ξ`` is the cell-specific dilution rate.
Together with the thermodinamic bound we have

stst_bound(c, ξ, b) = sign(b) * min(c/ ξ, abs(b))

"""
stst_bound(c, ξ, b) = sign(b) * min(c/ ξ, abs(b))