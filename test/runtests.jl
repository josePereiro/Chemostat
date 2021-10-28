using Test
import Pkg


@testset "Chemostat.jl" begin

    for pkg in ["MetNets", "MetLP", "MetEP"]
        println("-"^60)
        println("Tesing: ", pkg)
        Pkg.test(pkg)
        println()
    end
end
