
using StaticArrays
using LinearAlgebra

struct Cube{T}
    origin::SVector{3,T}
    shape::Tuple{Int,Int,Int}
    inverse::SMatrix{3,3,T,9}
    density::Array{T,3}
end

function Cube(filename, ::Type{T}=Float64) where {T<:AbstractFloat}

    data = readdlm(filename; skipstart=2)
    natoms = convert(Int, data[1,1])
    origin = SVector{3,T}(data[1,2:4])
    shape = Tuple(data[2:4])
    cell = SMatrix{3,3,T,9}(data[2:4,2:4]) .* shape
    inverse = inv(cell)

    volumetric_data = permutedims(data[5+natoms:end,:], [2,1])
    indices = volumetric_data .!= ""
    density = convert(Vector{T}, volumetric_data[indices])
    density = permutedims(reshape(density, reverse(shape)), 3:-1:1)

    Cube(origin, shape, inverse, density)
end

function (cube::Cube)(r::AbstractVector)
    r = r - cube.origin
    r = cube.inverse' * r

    if any(r .> 1) || any(r .< 0)
        throw(DomainError(r, "You cannot evaluate the density outside of the cell."))
    end
    indices = floor.(Int, r .* cube.shape) .+ 1

    return cube.density[indices...]
end
