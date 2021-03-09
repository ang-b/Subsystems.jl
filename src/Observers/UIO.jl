mutable struct UIO{V<:Real}
    F::Matrix{V}
    T::Matrix{V}
    B::Matrix{V}
    H::Matrix{V}
    C::OutputMatrix{V}
    A1::Matrix{V}
    K1::Matrix{V}
    K::Matrix{V}
    z::Vector{V}

    function UIO(A::Matrix{V}, 
                B::Matrix{V}, 
                C::Subsystems.OutputMatrix{V}, 
                E::Matrix{V}) where {V<:Real}
        (ny, feasible, lrE) = checkUIOFeasibility(A,C,E);
        if feasible 
            H = lrE*pinv(C*lrE)
            T = I - H*C
            A1 = T*A
            nz = size(A1, 1)
            # K2 = F*H;
            new{V}(
                zeros(V, nz, nz), 
                T, 
                B, 
                H, 
                C, 
                A1, 
                zeros(V, nz, ny), 
                zeros(V, nz, ny), 
                zeros(V, nz));
        else
            error("unfeasible UIO");
        end
    end
end

function UIO(A::Matrix{V}, 
        B::Matrix{V}, 
        C::Subsystems.OutputMatrix{V}, 
        E::Matrix{V},
        K1::Matrix{V}) where {V<:Real}
    uio = UIO(A, B, C, E)
    setF(uio, K1)
    uio    
end

function checkUIOFeasibility(A::Matrix{V}, 
        C::OutputMatrix{V}, 
        E::Matrix{V}) where {V<:Real}
    rE = rank(E);
    nx = size(A,1);
    if rE < nx
        SVD = svd(E);
        lrE = SVD.U[:,1:rE] * sqrt(SVD.S[1:rE, 1:rE]);
    else
        lrE = E;
    end
    matC = C == I ? 1.0 * Matrix(I, nx, nx) : C
    size(matC, 1), rank(matC*lrE) == rE, lrE 
end

function setF(o::UIO{T}, K1::Matrix{T}) where {T<:Real} 
    o.F = o.A1 - K1*o.C
    o.K1 = K1;
    o.K = K1 + o.F*o.H
end

function stepObs(o::UIO{T}, u::Vector{T}, y::Vector{T}) where {T<:Real}
    o.z = o.F*o.z + o.T*o.B*u + o.K*y
end

stepObs(o::UIO{T}, u::T, y::Vector{T}) where {T<:Real} = stepObs(o, [u], y)

outputMap(o::UIO{T}, y::Vector{T}) where {T<:Real} = o.z + o.H*y

function Base.getproperty(sys::UIO, s::Symbol)
    if s === :K2
        return sys.K - sys.K1
    else
        return getfield(sys, s)
    end
end

function Base.show(io::IO, sys::UIO)
    # Compose the name vectors
    #inputs = format_names(s.inputnames, "u", "?")
    #outputs = format_names(s.outputnames, "y", "?")
    println(io, "UIO:")
    # println(io, typeof(sys))
    println(io, "F = \n", _string_mat_with_headers(sys.F))
    println(io, "T = \n", _string_mat_with_headers(sys.T))
    println(io, "B = \n", _string_mat_with_headers(sys.B))
    println(io, "K = \n", _string_mat_with_headers(sys.K))
    println(io, "H = \n", _string_mat_with_headers(sys.H))
end
