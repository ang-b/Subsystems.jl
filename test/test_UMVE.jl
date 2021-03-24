module UMVETests

using LinearAlgebra
using Subsystems
using Test

@testset "umve estimator" begin
    A = Matrix{Float64}(I, 2, 2)
    B = reshape([1., 0.], :, 1)
    C = I
    Aij = [0.1 0; 0 0]
    W = Diagonal([0.1, 0.1])
    V = W 
    
    @testset "handles empty neighbourhoods" begin
        s = LdtiSubsystem(1, A, B, C)
        o = UMVE(s.A, s.B, s.C, s.E, W, V)
        Rk = Subsystems.R(o.A, o.C, o.Pi, o.W, o.V)
        M∅ = Subsystems.M(o.G, o.C, Rk)
        Kk = Subsystems.K(o.A, o.C, o.Pi, o.W, Rk)
        Ā∅ = Subsystems.barA(Kk, o.C, o.G, M∅)
        L̄∅ = Subsystems.barL(Kk, o.C, o.G, M∅)
        
        @test M∅ == zeros(0,2)
        @test Ā∅ == (I - Kk*o.C)        
        @test L̄∅ == Kk
    end
    
    @testset "Constructor" begin
        @testset "builds correctly if exists" begin
            N = 4
            LSS = [LdtiSubsystem(i, A, B, C) for i in 1:N]
            # make ring topology
            [addNeighbour(LSS[i], LSS[i%N + 1], Aij) for i in 1:N]
            O = [UMVE(s.A, s.B, s.C, s.E, W, V) for s in LSS]

            for o in O
                @test o.G ≈ reshape([1. 0], :, 1)
            end
        end
        @testset "can build with no neighbours" begin
            s = LdtiSubsystem(1, A, B, C)
            o = UMVE(s.A, s.B, s.C, s.E, W, V)
            @test o.G == zeros(s.nx, 0)
        end
    end

end

end