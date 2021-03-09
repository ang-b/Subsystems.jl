module Subsystems

using LinearAlgebra

abstract type AbstractSubsystem end
abstract type LinearSubsystem{T} <: AbstractSubsystem end

export OutputMatrix

OutputMatrix{T} = Union{typeof(I), Matrix{T}}

include("LdtiSubsystem.jl")
include("Observers/UIO.jl")

end
