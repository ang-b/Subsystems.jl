const LSS{T} = Vector{<:AbstractSubsystem{T}}

Base.sort(s::LSS) = sort(s, by = x -> x.index) 

function _getdimarrays(subs::LSS)
    N = length(subs)
    n = zeros(Int,N)
    m = zeros(Int,N)
    p = zeros(Int,N)
    for i in 1:N
        n[i] = subs[i].nx
        m[i] = subs[i].nu
        p[i] = subs[i].ny
    end
    (n,m,p)
end

function tomonolith(subs::LSS{T}) where T
    # sorted = sort(subs)
    n,m,p = _getdimarrays(subs)
    A = BlockArray{T}(zeros(sum(n), sum(n)), n, n)
    B = BlockArray{T}(zeros(sum(n), sum(m)), n, m)
    C = BlockArray{T}(zeros(sum(p), sum(n)), p, n)
    for s in subs 
        i = s.index
        view(A, Block(i,i)) .= s.A
        view(B, Block(i,i)) .= s.B
        view(C, Block(i,i)) .= s.C
        for j in 1:length(s.neighbours)
           view(A, Block(i, s.neighbours[j].index)) .= s.Aij[j] 
        end
    end
    (A=A, B=B, C=C)
end

