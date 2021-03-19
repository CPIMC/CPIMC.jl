using Test

include("../src/Configuration.jl")
include("../src/UEG/model.jl")


tests = [ "test_configuration.jl", "UEG/test_model.jl" ]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
