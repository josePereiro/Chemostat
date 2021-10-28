using Test
import Pkg


@testset "Chemostat.jl" begin

    for pkg in ["MetNets", "MetLP", "MetEP"]
        println("-"^60)
        @info("Tesing", pkg)
        Pkg.test(pkg)
        println()
    end
end
