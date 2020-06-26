function test_linear_S()
    S =[1.0 -1.0  0.0  0.0;
        0.0  1.0 -1.0  0.0;
        0.0  0.0  1.0 -1.0]
    @test all(Chemostat.Test.lineal_S(size(S,1)) .== S)
end
test_linear_S()

