"""
This uses a cube file to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""
module CubeLDFAModel

using DataInterpolations: CubicSpline
using DelimitedFiles: readdlm
using UnitfulAtomic: austrip
using Unitful: @u_str
using NQCBase: PeriodicCell, apply_cell_boundaries!
using NQCBase
using NQCModels: NQCModels, FrictionModels, Model


include("cube.jl")

export LDFAModel

"""
    LDFAModel(model::Model, filename, atoms, cell;
              friction_atoms=collect(range(atoms)),
              )

Wrapper for existing models that adds LDFA friction.

This model uses a cube file to evaluate the electron density used to calculate the friction.
This model assumes that the cube file has units of bohr for the grid and cell distances,
but provides the density in ``Å^{-3}``, as is the default in `FHI-aims`.
"""
struct LDFAModel{T,M,S,C,L,D,LM,A} <: FrictionModels.AdiabaticFrictionModel
    "Model that provides the energies and forces."
    model::M
    "Splines fitted to the numerical LDFA data."
    splines::S
    "Cube file reader."
    cube::C
    "Temporary array for storing the electron density."
    ρ::Vector{T}
    "Temporary array for the Wigner-Seitz radii."
    radii::Vector{T}
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
    "Periodic cell containing the electron density."
    cell::PeriodicCell{T}
    "LDFA model type" 
    ldfa_model::L
    "Descriptors used for scikit models"
    descriptors::D
    "Loaded ML model"
    loaded_model::LM
    "Atoms"
    atoms::A
    function LDFAModel(model, splines, cube, ρ, radii, friction_atoms, cell, ldfa_model, descriptors, loaded_model, atoms)
        new{eltype(cell),typeof(model),typeof(splines),typeof(cube),typeof(ldfa_model),typeof(descriptors),typeof(loaded_model),typeof(atoms)}(
            model, splines, cube, ρ, radii, friction_atoms, cell, ldfa_model, descriptors, loaded_model, atoms)
    end
end

function LDFAModel(model::Model, atoms, cell, ldfa_model; cube_filename="", friction_atoms=collect(range(atoms)), cell_matching_rtol=1e-3, descriptors=nothing, loaded_model=nothing)

    if ldfa_model == :cube
        cube = Cube(cube_filename)
        if !isapprox(cell.vectors, cube.cell.vectors, rtol=cell_matching_rtol)
            error("the cube file cell vectors do not match the simulation cell vectors.\n",
                "  Simulation vectors: ", cell.vectors, "\n",
                "  Cube vectors: ", cube.cell.vectors
            )
        end
    elseif ldfa_model == :scikit
        cube = nothing
    else
        error("The LDFA model: ", ldfa_model, " is not recognised.\n",
            "Please, choose a different model."
        )   
    end

    ldfa_data, _ = readdlm(joinpath(@__DIR__, "ldfa.txt"), ',', header=true)
    r = ldfa_data[:,1]
    splines = []
    for i in range(atoms)
        η = ldfa_data[:,atoms.numbers[i].+1]
        indices = η .!= ""
        ri = convert(Vector{Float64}, r[indices])
        η = convert(Vector{Float64}, η[indices])
        push!(ri, 10.0) # Ensure it goes through 0.0 for large r.
        push!(η, 0.0)
        push!(splines, CubicSpline(η, ri))
    end

    ρ = zeros(length(atoms))
    radii = zero(ρ)

    LDFAModel(model, splines, cube, ρ, radii, friction_atoms, cell, ldfa_model, descriptors, loaded_model, atoms)
end
NQCModels.ndofs(model::LDFAModel) = NQCModels.ndofs(model.model)

function NQCModels.potential(model::LDFAModel, R::AbstractMatrix)
    NQCModels.potential(model.model, R)
end

function NQCModels.derivative!(model::LDFAModel, D::AbstractMatrix, R::AbstractMatrix)
    NQCModels.derivative!(model.model, D, R)
end

function density!(model::LDFAModel, ρ::AbstractVector, R::AbstractMatrix)
    if model.ldfa_model == :cube
        for i in model.friction_atoms
            r = R[:,i]
            apply_cell_boundaries!(model.cube.cell, r)
            ρ[i] = austrip(model.cube(r + model.cube.origin) * u"Å^-3")
        end
    elseif model.ldfa_model == :scikit
        for i in model.friction_atoms
            ase_atoms = NQCBase.convert_to_ase_atoms(model.atoms, R, model.cell)
            r_desc = model.descriptors.create(ase_atoms, positions=[i-1], n_jobs=1) #n_threads)
            ρ[i] = model.loaded_model.predict(r_desc)[end]
        end
    end
end

function FrictionModels.friction!(model::LDFAModel, F::AbstractMatrix, R::AbstractMatrix)
    density!(model, model.ρ, R)
    clamp!(model.ρ, 0, Inf)
    @. model.radii = 1 / cbrt(4/3 * π * model.ρ)
    DoFs = size(R, 1)
    for i in 1:length(model.friction_atoms)
        η = model.radii[model.friction_atoms[i]] < 10 ? model.splines[model.friction_atoms[i]](model.radii[model.friction_atoms[i]]) : 0.0
        for j in axes(R, 1)
            F[(i-1)*DoFs+j, (i-1)*DoFs+j] = η
        end
    end
    return F
end

end
