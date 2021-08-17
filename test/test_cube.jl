
using Test
using NonadiabaticDynamicsBase
using CubeLDFAModel
using PyCall

py"""
import sys
sys.path.insert(0, ".")
"""

cube = pyimport("cube")

filename = "test.cube"

c = CubeLDFAModel.Cube(filename)
cube_object = cube.cube()
cube_object.read(filename)

@testset "Compare to cube.py" begin
    @test c.density == cube_object."density".reshape(cube_object.x_len, cube_object.y_len, cube_object.z_len)

    vecs = inv(c.inverse)
    for _=1:10
        r = c.origin
        r += rand() * vecs[1,:]
        r += rand() * vecs[2,:]
        r += rand() * vecs[3,:]

        @test c(r) ≈ cube_object(au_to_ang.(r)...)
    end
end

