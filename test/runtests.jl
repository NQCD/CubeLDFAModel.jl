using NonadiabaticModels: Free
using NonadiabaticDynamicsBase: PeriodicCell, Atoms
using CubeLDFAModel
using Test

@testset "CubeLDFAModel.jl" begin

    vecs = [0.111175    0.000000    0.000000
           -0.055588    0.096280    0.000000
            0.000000    0.000000    0.703079]' .* 100
    model = LDFAModel(Free(), "test.cube", Atoms([:H]), PeriodicCell(vecs))
end

include("test_cube.jl")
