module UMVETests

using LinearAlgebra
using Subsystems
using Test

@testset "umve estimator" begin
    A = Matrix{Float64}(I, 2, 2)
    B = reshape([1., 0.], :, 1)
    C = I
    Aij = [0.1 0; 0 0]
    
    @testset "builds correctly" begin
        N = 4
        W = Diagonal([0.1, 0.1])
        V = W 
        LSS = [LdtiSubsystem(i, A, B, C) for i in 1:N]
        # make ring topology
        [addNeighbour(LSS[i], LSS[i%N + 1], Aij) for i in 1:N]
        O = [UMVE(s.A, s.B, s.C, s.E, W, V) for s in LSS]

        for o in O
            @test o.G ≈ reshape([1. 0], :, 1)
        end
    end


end

end