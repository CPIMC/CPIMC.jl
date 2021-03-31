using Test

tests = [ "test_configuration.jl", "test_planewave.jl", "test_UniformElectronGas.jl" ]


for t in tests
    @testset "$t" begin
        include(t)
    end
end
