module Subsystems

using LinearAlgebra

abstract type AbstractSubsystem end
abstract type LinearSubsystem{T} <: AbstractSubsystem end
 
export LdtiSubsystem, 
    LinearSubsystem, 
    AbstractSubsystem,
    addNeighbour, 
    removeNeighbour,
    getNeighbourIndices,
    stepSubsystem,
    outputMap,
    OutputMatrix,
    UIO, 
    setF,
    stepObs

OutputMatrix{T} = Union{typeof(I), Matrix{T}}

function _string_mat_with_headers(X::VecOrMat)
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

include("LdtiSubsystem.jl")
include("Observers/UIO.jl")

end
