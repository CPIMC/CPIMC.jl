include("../../src/UEG/estimators.jl")

e = CEnsemble(2,5.68089,7)

@testset "Atomic_units" begin
    for _ in 1:100
        local E = rand()*100
        @test Et_Ry(E, e) == 2*Et_Ha(E, e)
    end
end
