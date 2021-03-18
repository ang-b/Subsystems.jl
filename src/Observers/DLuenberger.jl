mutable struct DLUE{V<:Real}
    F::Matrix{V}
    B::Matrix{V}
    C::OutputMatrix{V}
    L::Matrix{V}
    E::Matrix{V}
    z::Vector{V}

    function DLUE(A::T, B::T, C::OutputMatrix{V}, L::T) where {V<:Real, T<:Matrix{V}}
        new{V}(A - L*C, 
            B, 
            C,
            L, 
            Matrix{V}(undef, size(A,1), 0),
            zeros(V, size(A,1)))
    end
end

function setInitialState(o::DLUE{V}, z0::Vector{V}) where {V<:Real}
    setfield!(o, :z, z0)
end

function setNeighbourCoupling(o::DLUE{V}, coupling::AbstractMatrix{V}) where {V} 
    setfield!(o, :E, isempty(coupling) ? Matrix{V}(undef, o.nx, 0) : coupling)
end

function updateState(o::DLUE{V}, u::Vector{V}, y::Vector{V}, hatbarx::AbstractVector{V}) where {V<:Real}
    if isempty(hatbarx)
        setfield!(o, :z, o.F*o.z + o.B*u + o.L*y)
    else
        setfield!(o, :z, o.F*o.z + o.B*u + o.L*y + o.E*hatbarx)
    end
end

updateState(o::DLUE{V}, u::V, y::Vector{V}, hatbarx::Vector{V}) where {V} = updateState(o, [u], y, hatbarx)

outputMap(o::DLUE{<:Real}) = o.C*o.z

function Base.getproperty(o::DLUE, s::Symbol) 
    if s === :nx
        return size(o.F, 1)
    else
        getfield(o, s)
    end
end

function Base.setproperty!(::DLUE, ::Symbol, val) 
    error("private object")
end