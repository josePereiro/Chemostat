function test_TestUtils_linear_S()
    S =[1.0 -1.0  0.0  0.0;
        0.0  1.0 -1.0  0.0;
        0.0  0.0  1.0 -1.0]
    @test  all(Chemostat.Utils.lineal_S(size(S,1)) .== S)
end
test_TestUtils_linear_S()

