
function grad_desc(f;
        target::Vector,
        x0::Vector, x1::Vector,
        C::Vector = abs.(x0 - x1), th = 1e-5,
        maxiters = 1000, 
        toshow::Vector = [],
        Err = (fᵢ) -> abs.(target .- fᵢ) ./ ifelse.(iszero.(target), 1.0, abs.(target)),
        verbose = true
    )

    # initializing
    xᵢ₋₁, xᵢ = x0, x1
    fᵢ₋₁ = f(xᵢ₋₁)
    ϵᵢ₋₁ = Err(fᵢ₋₁)
    Δx = zero(xᵢ)
    sense = ones(length(target))

    verbose && (prog = ProgressThresh(th, "Grad desc: "))
    for it in 1:maxiters
        
        xᵢ += Δx
        fᵢ = f(xᵢ)
        ϵᵢ = Err(fᵢ)
        sense .*= -sign.(ϵᵢ .- ϵᵢ₋₁)
        sense == 0.0 && error("sense == 0. Descend gets stocked, target unreachable!!!")
        Δx = sense .* C .* ϵᵢ

        maxϵᵢ = maximum(ϵᵢ)
        maxϵᵢ < th && break

        xᵢ₋₁ =  xᵢ
        ϵᵢ₋₁ = ϵᵢ

        verbose && update!(prog, maxϵᵢ; showvalues = vcat(
                [
                    ("it", it),
                    ("maxϵᵢ", maxϵᵢ),
                    ("ϵᵢ", ϵᵢ),
                    ("sense", sense),
                    ("xᵢ", xᵢ),
                    ("t", target),
                    ("fᵢ", fᵢ),
                ], toshow
            )
        )
    end
    verbose && finish!(prog)

    return xᵢ
end