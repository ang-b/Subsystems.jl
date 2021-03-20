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

function setInitialState(sys::LdtiSubsystem{T}, x::Vector{T}) where {T}
    setfield!(sys, :x, x)
    setfield!(sys, :xnext, x)
end

function updateState(sys::LdtiSubsystem{T}, u::Vector{T}, w::Vector{T}) where {T}
    if isempty(sys.neighbours) 
        setfield!(sys, :xnext, sys.A * sys.x + sys.B * u .+ w)
    else
        barx = mapreduce(s -> s.x, vcat, sys.neighbours)
        setfield!(sys, :xnext, sys.A * sys.x + sys.B * u + sys.E * barx + w) 
    end
end
updateState(sys::LdtiSubsystem{T}, u::T, w::Vector{T}) where {T} = updateState(sys, [u], w)
updateState(sys::LdtiSubsystem{T}, u) where {T} = updateState(sys, u, convert(T, 0))

function outputMap(sys::LinearSubsystem{T}, u) where {T<:Real}
    setfield!(sys, :x, sys.xnext)
    sys.C * sys.x + sys.D * u 
end

function addNeighbour(sys::LinearSubsystem{T}, neighbour::LinearSubsystem{T}, weight::Matrix{T}, conservation::Bool = false) where {T<:Real}
    if neighbour.index == sys.index 
        error("Self loops are not allowed")
    end
    push!(sys.neighbours, neighbour)
    push!(sys.Aij, weight)
    updateE(sys)
    if conservation
       sys.A -= weight 
    end
    map(x-> convert(Int, x), getNeighbourIndices(sys))
end

function removeNeighbour(sys::AbstractSubsystem, neighbour::Integer, conservation::Bool = false) 
    idx = findfirst(x -> x.index == neighbour, sys.neighbours)
    deleteat!(sys.neighbours, idx)
    weight = splice!(sys.Aij, idx)
    if conservation
        sys.A += weight
    end
end

function setNeighbour(sys::LinearSubsystem, neighbour::LinearSubsystem, neighbourIdx::Integer, weight::Matrix{T}, conservation::Bool = false) where {T} 
    j = findfirst(x -> x.index == neighbourIdx, sys.neighbours)
    # if conservation 
    #     sys.A -= weight
    # end
    sys.neighbours[j] = neighbour
    if conservation
        sys.A = sys.A + sys.Aij[j] - weight
    end
    sys.Aij[j] = weight
    updateE(sys)
end

function getNeighbourIndices(sys::AbstractSubsystem) 
    Tuple([s.index for s in sys.neighbours])
end

function updateE(sys::LinearSubsystem{T}) where {T}
    setfield!(sys, :E, isempty(sys.Aij) ? 
                            Matrix{T}(undef, sys.nx, 0) : 
                            foldl(hcat, sys.Aij))
end

ninputs(sys::LinearSubsystem) = size(sys.D, 2)
noutputs(sys::LinearSubsystem) = size(sys.D, 1)
nstates(sys::LinearSubsystem) = size(sys.A, 1)

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

function Base.show(io::IO, sys::LinearSubsystem)
    # Compose the name vectors
    #inputs = format_names(s.inputnames, "u", "?")
    #outputs = format_names(s.outputnames, "y", "?")
    println(io, "Subsystem $(sys.index):")
    println(io, "Neighbours = $(convert.(Int64, getNeighbourIndices(sys)))")
    # println(io, typeof(sys))
    if size(sys.A, 1) > 0
        #states = format_names(s.statenames, "x", "?")
        println(io, "A = \n", _string_mat_with_headers(sys.A))
        println(io, "B = \n", _string_mat_with_headers(sys.B))
        println(io, "C = \n", sys.C == I ? "I" : _string_mat_with_headers(sys.C))
    end
    println(io, "D = \n", _string_mat_with_headers(sys.D), "\n")
end
