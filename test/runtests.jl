using Test

include("../src/Configuration.jl")
include("../src/UEG/model.jl")


tests = [ "test_configuration.jl" ]

for t in tests
    @testset "$t" begin
        include(t)
    end
end

