export LdtiSubsystem, 
    LinearSubsystem, 
    AbstractSubsystem,
    addNeighbour, 
    removeNeighbour,
    getNeighbourIndices,
    stepSubsystem,
    outputMap

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
    index::UInt
    neighbours::Array{LdtiSubsystem{T}}
    Aij::Array{Matrix{T}}
    x::Vector{T}
    function LdtiSubsystem(i, A::Matrix{T}, B::Matrix{T}, C::OutputMatrix{T}, D::Matrix{T}) where {T<:Real} 
        (nx, nu, ny) = state_space_validation(A, B, C, D)
        new{T}(A, B, C, D, i, [], [], zeros(T, nx))
    end
    function LdtiSubsystem(i, A::Matrix{T}, B::Matrix{T}, C::OutputMatrix{T}) where {T<:Real}
        ny = C == I ? size(A,1) : size(C,1)
        D = zeros(T, ny, size(B, 2))    
        LdtiSubsystem(i, A, B, C, D)
    end
end

function stepSubsystemLocal(sys::LdtiSubsystem{T}, u::Vector{T}) where {T<:Real}
    sys.x = sys.A * sys.x + sys.B * u
end

stepSubsystemLocal(sys::LdtiSubsystem{T}, u::T) where {T<:Real} = stepSubsystemLocal(sys, [u])

function stepSubsystem(sys::LdtiSubsystem{T}, u::Vector{T}) where {T<:Real}
    if isempty(sys.neighbours) 
        stepSubsystemLocal(sys, u) 
    else
        barx = mapreduce(s -> s.x, vcat, sys.neighbours)
        stepSubsystemLocal(sys, u)
        sys.x += sys.E * barx 
    end
end

stepSubsystem(sys::LdtiSubsystem{T}, u::T) where {T<:Real} = stepSubsystem(sys, [u])

function outputMap(sys::LinearSubsystem{T}, u) where {T<:Real}
   y = sys.C * sys.x + sys.D * u 
end

function addNeighbour(sys::LinearSubsystem{T}, neighbour::LinearSubsystem{T}, weight::Matrix{T}) where {T<:Real}
    if neighbour.index == sys.index 
        error("Self loops are not allowed")
    end
    push!(sys.neighbours, neighbour)
    push!(sys.Aij, weight)
    map(x-> convert(Int, x), getNeighbourIndices(sys))
end

function removeNeighbour(sys::AbstractSubsystem, neighbour::UInt) 
    idx = findfirst(x -> x.index == neighbour, sys.neighbours)
    deleteat!(sys.neighbours, idx)
    deleteat!(sys.Aij, idx)
end

function getNeighbourIndices(sys::AbstractSubsystem) 
    Tuple([s.index for s in sys.neighbours])
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
    elseif s === :E
        return isempty(sys.Aij) ? zeros(typeof(sys.x), sys.nx, sys.nx) : foldl(hcat, sys.Aij)
    else
        return getfield(sys, s)
    end
end

function _string_mat_with_headers(X::VecOrMat)
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

function Base.show(io::IO, sys::LinearSubsystem)
    # Compose the name vectors
    #inputs = format_names(s.inputnames, "u", "?")
    #outputs = format_names(s.outputnames, "y", "?")
    println(io, "Subsystem $(sys.index):")
    # println(io, typeof(sys))
    if size(sys.A, 1) > 0
        #states = format_names(s.statenames, "x", "?")
        println(io, "A = \n", _string_mat_with_headers(sys.A))
        println(io, "B = \n", _string_mat_with_headers(sys.B))
        println(io, "C = \n", sys.C == I ? "I" : _string_mat_with_headers(sys.C))
    end
    println(io, "D = \n", _string_mat_with_headers(sys.D), "\n")
end
