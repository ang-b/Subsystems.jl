mutable struct UMVE{T<:Real}
    A::Matrix{T}
    B::Matrix{T}
    C::OutputMatrix{T}
    G::Matrix{T}
    W::Matrix{T}
    V::Matrix{T}
    Ebar::Matrix{T}
    z::Vector{T}
    xi::Vector{T}
    Pi::Matrix{T} # state estimation error covariance
    Delta::Union{Matrix{T}, Missing} # input estimation error covariance

    function UMVE(A::S,
        B::S,
        C::OutputMatrix{T},
        E::S,
        W::AbstractMatrix{T},
        V::AbstractMatrix{T},
        decompose::Bool = true) where {T<:Real, S<:Matrix{T}}
        if !isposdef(W)
            throw(DomainError(W, "W must be positive definite"))
        end
        if !isposdef(V)
            throw(DomainError(V, "V must be positive definite"))
        end
        W = Matrix(W)
        V = Matrix(V)
        feasible, G, Ebar = checkUMVEFeasibility(C,E, decompose)
        nx = size(A,1)
        nxi = size(G,2)
        if feasible
            Pk = W
            Rk = R(A, C, Pk, W, V)
            Mk = M(G, C, Rk)
            new{T}(
                A,
                B,
                C,
                G,
                W,
                V,
                Ebar,
                zeros(T, nx),
                zeros(T, nxi),
                W,
                missing
            )
        else
            error("unfeasible")
        end
    end
end

function checkUMVEFeasibility(C::OutputMatrix, E::Matrix, decompose::Bool) 
    rE = rank(E)
    feasible = rank(C*E) == rE
    if decompose
        SVD = svd(E)
        G = SVD.U[:,1:rE]
        Ebar = (Diagonal(SVD.S) * SVD.Vt)[1:rE,:]
    else
        nz = size(E,1)
        G = isempty(E) ? zeros(nz, 0) : Matrix{eltype(E)}(I, nz, nz)
        Ebar = E
    end
    feasible, G, Ebar
end

R(A::T, C::OutputMatrix{S}, Pk::T, W::T, V::T) where {S<:Real, T<:Matrix{S}} = C*(A*Pk*A' + W)*C' + V 
M(G::T, C::OutputMatrix{S}, Rk::T) where {S<:Real, T<:Matrix{S}} = inv( G'C'*inv(Rk)*C*G ) * G'C'inv(Rk) 
K(A::T, C::OutputMatrix{S}, Pk::T, W::T, Rk::T) where {S<:Real, T<:Matrix{S}} = C*(A*Pk*A' + W)*C'inv(Rk) 
barA(Kk::T, C::OutputMatrix{S}, G::T, Mk::T) where {S<:Real, T<:Matrix{S}} = (I - Kk*C)*(I - G*Mk*C)
barL(Kk::T, C::OutputMatrix{S}, G::T, Mk::T) where {S<:Real, T<:Matrix{S}} = Kk + (I - Kk*C)*G*Mk

function updateState(o::UMVE{T}, u::Vector{T}, y::Vector{T}) where {T}
    Rk = R(o.A, o.C, o.Pi, o.W, o.V)
    # setfield!(o, :R, Rk)
    Mk = M(o.G, o.C, Rk)
    setfield!(o, :Delta, Mk*Rk*Mk')
    Kk = K(o.A, o.C, o.Pi, o.W, Rk)
    Abar = barA(Kk, o.C, o.G, Mk)
    Lbar = barL(Kk, o.C, o.G, Mk)
    biased_z = o.A*o.z + o.B*u
    setfield!(o, :z, Abar*biased_z + Lbar*y)
    setfield!(o, :xi, Mk*(y - o.C*biased_z)) 
    setfield!(o, :Pi, Abar*o.A*o.Pi*o.A'Abar' + Abar*o.W*Abar' + Lbar*o.V*Lbar')
end
updateState(o::UMVE{T}, u::T, y::Vector{T}) where {T} = updateState(o, [u], y)

outputMap(o::UMVE) = o.C*o.z

function Base.getproperty(o::UMVE, s::Symbol)
    if s === :nz
        return length(o.z)
    elseif s === :nxi
        return length(o.xi)
    else
        return getfield(o, s)
    end
end