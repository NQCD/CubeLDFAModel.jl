using NQCModels: Free, FrictionModels
using PyCall
using NQCBase
using CubeLDFAModel
using Test
using Unitful

include("test_cube.jl")

aseio = pyimport("ase.io")
ase_atoms = aseio.read("start.in")
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)

model = LDFAModel(Free(), "test.cube", atoms, cell; friction_atoms=[1, 2])

@testset "friction!" begin

    F = zeros(6, 6)
    r = @view R[:,1]
    r .= 0
    FrictionModels.friction!(model, F, R)

    for _=1:10
        r .= 0
        r += rand() * cell.vectors[:,1]
        r += rand() * cell.vectors[:,2]
        r += rand() * cell.vectors[:,3]

        FrictionModels.friction!(model, F, R)
    end

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3]
    FrictionModels.friction!(model, F, R)

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3] + rand(3)
    FrictionModels.friction!(model, F, R)
end
