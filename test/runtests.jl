include("../src/Configuration.jl")
include("../src/UEG/model.jl")
include("../src/Updates/Ideal-Updates.jl")
include("../src/Updates/Other-Updates.jl")
include("../src/Updates/Type-A-Updates.jl")
include("../src/Updates/Type-B-Updates.jl")
include("../src/Updates/Type-C-Updates.jl")
include("../src/Updates/Type-D-Updates.jl")
include("../src/Updates/Type-E-Updates.jl")
include("../src/UEG/estimators.jl")
include("../src/CPIMC.jl")

using Test

include("../src/Configuration.jl")
include("../src/UEG/model.jl")


tests = [ "test_configuration.jl", "UEG/test_model.jl", "UEG/test_estim.jl"]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
