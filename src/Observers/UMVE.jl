mutable struct UMVE{T<:Real}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    G::Matrix{T}
    W::Matrix{T}
    V::Matrix{T}
    # Abar::Matrix{T}
    Ebar::Matrix{T}
    # Lbar::Matrix{T}
    z::Vector{T}
    xi::Vector{T}
    Pi::Matrix{T} # state estimation error covariance
    Delta::Union{Matrix{T}, Missing} # input estimation error covariance

    function UMVE(A::S, B::S, C::OutputMatrix{T}, E::S, W::S, V::S) where {T<:Real, S<:Matrix{T}}
        if !isposdef(W)
            throw(DomainError(W, "W must be positive definite"))
        end
        if !isposdef(V)
            throw(DomainError(V, "V must be positive definite"))
        end
        nx = size(A,1)
        nxi = size(G,2)
        feasible, G, Ebar = checkUMVEFeasibility(C,E)
        if feasible
            Pk = W
            Rk = R(A, C, Pk, W, V)
            Mk = M(G, C, Rk)
            # Kk = K(A, C, Pk, W, Rk) 
            new{T}(
                A,
                B,
                C,
                G,
                W,
                V,
                # barA(Kk, C, G, Mk),
                Ebar,
                # barL(Kk, C, G, Mk),
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

function checkUMVEFeasibility(C::OutputMatrix{T}, E::Matrix{T}) where {T}
    feasible = rank(C*E) == rank(E)
    SVD = svd(E)
    sqrtSigma = sqrt(Diagonal(SVD.s))
    G = SVD.U * sqrtSigma
    Ebar = sqrtSigma * SVD.Vt
    feasible, G, Ebar
end

R(A::T, C::OutputMatrix{S}, Pk::T, W::T, V::T) where {S<:Real, T<:Matrix{S}} = C*(A*Pk*A' + W)*C' + V 
M(G::T, C::OutputMatrix{S}, Rk::T) where {S<:Real, T<:Matrix{S}} = inv( G'C'*inv(Rk)*C*G ) * G'C'inv(Rk) 
K(A::T, C::OutputMatrix{S}, Pk::T, W::T, Rk::T) where {S<:Real, T<:Matrix{S}} = C*(A*Pk*A' + W)*C'inv(Rk) 
barA(Kk::T, C::OutputMatrix{S}, G::T, Mk::T) where {S<:Real, T<:Matrix{S}} = (I - Kk*C)*(I - G*Mk*C)
barL(Kk::T, C::OutputMatrix{S}, G::T, Mk::T) where {S<:Real, T<:Matrix{S}} = Kk + (I - Kk*C)*G*Mk

function updateState(o::UMVE{T}, u::Vector{T}, y::Vector{T}) where {T}
    Rk = R(o.A, o.C, o.Pi, o.W, o.V)
    setfield!(o, :R, Rk)
    Mk = M(o.G, o.C, Rk)
    setfield!(o, :Delta, Mk*Rk*Mk')
    Kk = K(o.A, o.C, Pk, o.W, Rk)
    Abar = barA(Kk, o.C, o.G, Mk)
    Lbar = barL(Kk, o.C, o.G, Mk)
    biased_z = o.A*o.z + o.B*u
    setfield!(o, :z, Abar*biased_z + Lbar*y)
    setfield!(o, :xi, Mk*(y - o.C*biased_z)) 
    setfield!(o, :Pi, Abar*o.A*o.Pi*A'Abar' + Abar*W*Abar' + Lbar*V*Lbar')
end