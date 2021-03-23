function state_space_validation(A,B,C,D)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = C == I ? nx : size(C, 1)

    if size(A, 2) != nx && nx != 0
        error("A has dimentions $(size(A)), but must be square")
    elseif size(B, 1) != nx
        error("The number of rows of A ($(size(A,1))) and B ($(size(B,1))) are not equal")
    elseif ny != nx
        error("The number of columns of A ($(size(A,2))) and C ($(ny)) are not equal")
    elseif nu != size(D, 2)
        error("The number of columns of B ($(size(B,2))) and D ($(size(D,2))) are not equal")
    elseif ny != size(D, 1)
        error("The number of rows of C ($(ny)) and D ($(size(D,1))) are not equal")
    end

    nx,nu,ny
end

mutable struct LdtiSubsystem{T<:Real} <:LinearSubsystem{T}
    A::Matrix{T}
    B::Matrix{T}
    C::OutputMatrix{T}
    D::Matrix{T}
    E::Matrix{T}
    index::UInt
    neighbours::Array{LdtiSubsystem{T}}
    Aij::Array{Matrix{T}}
    x::Vector{T}
    xnext::Vector{T}
    function LdtiSubsystem(i, A::Matrix{T}, B::Matrix{T}, C::OutputMatrix{T}, D::Matrix{T}) where {T<:Real} 
        (nx, nu, ny) = state_space_validation(A, B, C, D)
        new{T}(A, B, C, D, Matrix{T}(undef, nx, 0), i, [], [], zeros(T, nx), zeros(T, nx))
    end
    function LdtiSubsystem(i, A::Matrix{T}, B::Matrix{T}, C::OutputMatrix{T}) where {T<:Real}
        ny = C == I ? size(A,1) : size(C,1)
        D = zeros(T, ny, size(B, 2))    
        LdtiSubsystem(i, A, B, C, D)
    end
end

function setInitialState(s::LdtiSubsystem{T}, x::Vector{T}) where {T}
    setfield!(s, :x, x)
    setfield!(s, :xnext, x)
end

function updateState(s::LdtiSubsystem{T}, u::Union{T, Vector{T}}, w::Vector{T}) where {T}
    if isempty(s.neighbours) 
        setfield!(s, :xnext, s.A * s.x + (s.B * u)[:] .+ w)
    else
        barx = getNeighbourStates(s)
        setfield!(s, :xnext, s.A * s.x + s.B * u + s.E * barx + w) 
    end
end
updateState(s::LdtiSubsystem{T}, u::Union{T, Vector{T}}) where {T} = updateState(s, u, zeros(T, s.nx))

function outputMap(s::LinearSubsystem{T}, u::Union{T, Vector{T}}, v::Vector{T}) where {T}
    setfield!(s, :x, s.xnext)
    s.C * s.x + (s.D * u)[:] + v
end
outputMap(s::LinearSubsystem{T}, u::Union{T, Vector{T}}) where {T} = outputMap(s, u, zeros(T, s.ny))

getNeighbourStates(s::AbstractSubsystem) = mapreduce(s -> s.x, vcat, s.neighbours)

function addNeighbour(s::LinearSubsystem{T}, neighbour::LinearSubsystem{T}, weight::Matrix{T}, conservation::Bool = false) where {T<:Real}
    if neighbour.index == s.index 
        error("Self loops are not allowed")
    end
    push!(s.neighbours, neighbour)
    push!(s.Aij, weight)
    updateE(s)
    if conservation
       s.A -= weight 
    end
    map(x-> convert(Int, x), getNeighbourIndices(s))
end

function removeNeighbour(s::AbstractSubsystem, neighbour::Integer, conservation::Bool = false) 
    idx = findfirst(x -> x.index == neighbour, s.neighbours)
    deleteat!(s.neighbours, idx)
    weight = splice!(s.Aij, idx)
    if conservation
        s.A += weight
    end
end

function setNeighbour(s::LinearSubsystem, neighbour::LinearSubsystem, neighbourIdx::Integer, weight::Matrix{T}, conservation::Bool = false) where {T} 
    j = findfirst(x -> x.index == neighbourIdx, s.neighbours)
    s.neighbours[j] = neighbour
    if conservation
        s.A = s.A + s.Aij[j] - weight
    end
    s.Aij[j] = weight
    updateE(s)
end

function getNeighbourIndices(s::AbstractSubsystem) 
    Tuple([s.index for s in s.neighbours])
end

function updateE(s::LinearSubsystem{T}) where {T}
    setfield!(s, :E, isempty(s.Aij) ? 
                            Matrix{T}(undef, s.nx, 0) : 
                            foldl(hcat, s.Aij))
end

ninputs(s::LinearSubsystem) = size(s.D, 2)
noutputs(s::LinearSubsystem) = size(s.D, 1)
nstates(s::LinearSubsystem) = size(s.A, 1)

Base.eltype(::Type{S}) where {S<:LinearSubsystem} = S

function Base.getproperty(sys::LinearSubsystem, s::Symbol)
    if s === :nx
        return nstates(sys)
    elseif s === :nu
        return ninputs(sys)
    elseif s === :ny
        return noutputs(sys)
    elseif s === :Cmat
        return sys.C == I ? Matrix(I, sys.ny, sys.ny) : sys.C
    else
        return getfield(sys, s)
    end
end

function Base.setproperty!(sys::LinearSubsystem, s::Symbol, val)
    if s ∉ [:A, :B, :C, :D, :index] 
        error("Protected property")
    else
        if size(val) != size(getfield(sys, s))
            error("Cannot change size")
        else
            return setfield!(sys, s, val)
        end
    end
end

function Base.show(io::IO, s::LinearSubsystem)
    println(io, "Subsystem $(s.index):")
    if size(s.A, 1) > 0
        println(io, "A = \n", _string_mat_with_headers(s.A))
        println(io, "B = \n", _string_mat_with_headers(s.B))
        println(io, "C = ", s.C == I ? "I" : "\n"*_string_mat_with_headers(s.C))
    end
    println(io, "D = \n", _string_mat_with_headers(s.D), "\n")
    if !isempty(s.neighbours)
        println(io, "Neighbours = $(convert.(Int64, getNeighbourIndices(s)))")
        println(io, "E = \n", _string_mat_with_headers(s.E))
    end
end
