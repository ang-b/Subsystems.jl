module Subsystems

using LinearAlgebra
using BlockArrays

abstract type AbstractSubsystem{T} end
abstract type LinearSubsystem{T} <: AbstractSubsystem{T} end
 
export LdtiSubsystem, 
    LinearSubsystem, 
    OutputMatrix,
    AbstractSubsystem,
    UIO, 
    DLUE,
    UMVE,
    addNeighbour, 
    setNeighbour,
    setNeighbourCoupling,
    removeNeighbour,
    getNeighbourIndices,
    getNeighbourStates,
    getOutboundWeights,
    outputMap,
    setF,
    updateState,
    setInitialState

OutputMatrix{T} = Union{typeof(I), AbstractMatrix{T}}

function _string_mat_with_headers(X::AbstractVecOrMat)
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

include("LdtiSubsystem.jl")
include("Observers/UIO.jl")
include("Observers/DLuenberger.jl")
include("Observers/UMVE.jl")
include("LSS.jl")

export _getdimarrays,
    tomonolith

end
