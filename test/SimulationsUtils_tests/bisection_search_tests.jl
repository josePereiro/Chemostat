function test_bisection_search() 
        n = 500
    x0 = fill(0.0, n)
    x1 = fill(1.0, n)
    t = clamp.(rand(n), 1e-3, 1.0 - 1e-3)
    for e in [10.0^i for i in -1:0.1:1]
        println("f(x) = x.^$(e)")
        f(x) = x.^e
        x = Chemostat.SimulationUtils.bisection_search(f, x0, x1, t; verbose = false)
        @test !isnothing(x) && all(isapprox.(f(x), t; atol = 1e-4))
    end
end
test_bisection_search() 