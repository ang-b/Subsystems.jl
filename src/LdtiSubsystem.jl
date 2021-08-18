function state_space_validation(A,B,C,D)
    nx = size(A, 1)
    nu = size(B, 2)
    ny = C == I ? nx : size(C, 1)
    colC = C == I ? nx : size(C, 2)

    if size(A, 2) != nx && nx != 0
        error("A has dimentions $(size(A)), but must be square")
    elseif size(B, 1) != nx
        error("The number of rows of A ($(size(A,1))) and B ($(size(B,1))) are not equal")
    elseif colC != nx
        error("The number of columns of A ($(size(A,2))) and C ($(colC)) are not equal")
    elseif nu != size(D, 2)
        error("The number of columns of B ($(size(B,2))) and D ($(size(D,2))) are not equal")
    elseif ny != size(D, 1)
        error("The number of rows of C ($(ny)) and D ($(size(D,1))) are not equal")
    end

    nx,nu,ny
end

mutable struct LdtiSubsystem{T<:Real} <:LinearSubsystem{T}
    A::AbstractMatrix{T}
    B::AbstractMatrix{T}
    C::OutputMatrix{T}
    D::AbstractMatrix{T}
    E::AbstractMatrix{T}
    index::UInt
    neighbours::Array{LdtiSubsystem{T}}
    Aij::Array{AbstractMatrix{T}}
    x::Vector{T}
    xnext::Vector{T}
    function LdtiSubsystem(i, A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::OutputMatrix{T}, D::AbstractMatrix{T}) where {T<:Real} 
        (nx, nu, ny) = state_space_validation(A, B, C, D)
        new{T}(A, B, C, D, zeros(T, nx, 0), i, [], [], zeros(T, nx), zeros(T, nx))
    end
    function LdtiSubsystem(i, A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::OutputMatrix{T}) where {T<:Real}
        ny = C == I ? size(A,1) : size(C,1)
        D = zeros(T, ny, size(B, 2))    
        LdtiSubsystem(i, A, B, C, D)
    end
end

function Base.zero(s::LdtiSubsystem{T}) where {T}
    LdtiSubsystem(0, zeros(T, s.nx,s.nx), zeros(T, s.nx, 1), I)
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
    s.xnext
end
updateState(s::LdtiSubsystem{T}, u::Union{T, Vector{T}}) where {T} = updateState(s, u, zeros(T, s.nx))

function outputMap(s::LinearSubsystem{T}, u::Union{T, Vector{T}}, v::Vector{T}) where {T}
    setfield!(s, :x, s.xnext)
    s.C * s.x + (s.D * u)[:] + v
end
outputMap(s::LinearSubsystem{T}, u::Union{T, Vector{T}}) where {T} = outputMap(s, u, zeros(T, s.ny))

getNeighbourStates(s::AbstractSubsystem) = isempty(s.neighbours) ? zero(s.x) : foldl(vcat, [s.x for s in s.neighbours])

function addNeighbour(s::LinearSubsystem{T}, neighbour::LinearSubsystem{T}, weight::AbstractMatrix{T}, conservation::Bool = false) where {T}
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
    updateE(s)
end

function setNeighbour(s::LinearSubsystem, neighbour::LinearSubsystem, neighbourIdx::Integer, weight::AbstractMatrix{T}, conservation::Bool = false) where {T} 
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
    setfield!(s, :E, isempty(s.Aij) ? zeros(T, s.nx, 0) : foldl(hcat, s.Aij))
end

function getOutboundWeights(s::LinearSubsystem{T}) where {T}
    outWeights = Vector{AbstractMatrix{T}}()
    for sj in s.neighbours
        ji = findfirst(x -> s.index == x, getNeighbourIndices(sj))
        push!(outWeights, sj.Aij[ji])
    end
    outWeights
end

Base.eltype(::Type{S}) where {S<:LinearSubsystem} = S

function Base.getproperty(sys::LinearSubsystem{T}, s::Symbol) where {T}
    if s === :nx
        return size(sys.A, 1)
    elseif s === :nu
        return size(sys.B, 2)
    elseif s === :ny
        return size(sys.D, 1)
    elseif s === :nx̄
        return size(sys.E,2)
    elseif s === :Cmat
        return sys.C == I ? AbstractMatrix{T}(I, sys.ny, sys.ny) : sys.C
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
