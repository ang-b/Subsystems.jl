module ldtiSubTests

using LinearAlgebra
using Subsystems
using Test

@testset "ldti Subsystem" begin
    A = Matrix{Float64}(I, 2, 2);
    B = reshape([1., 0.], :, 1)
    C = I;
    
    @testset "defines correct system" begin
        s = LdtiSubsystem(1, A, B, C)
        @test s.A == A
        @test s.C == I
        @test s.B == B
    end
    
    @testset "implements correct conservation" begin
        s = LdtiSubsystem(1, A, B, C)
        s1 = LdtiSubsystem(2, A, B, C)
        k = [0.1 0; 0 0]

        addNeighbour(s, s1, k, true)
        @test s.A == A - k

        removeNeighbour(s, 2, true)
        @test s.A == A
        addNeighbour(s, s1, k, true)
        setNeighbour(s, s1, 2, 0.2*k, true)
        @test s.A == A - 0.2*k
        
        s = LdtiSubsystem(1, A, B, C)
        s2 = LdtiSubsystem(3, A, B, C)

        addNeighbour(s, s1, k, true)
        addNeighbour(s, s2, k, true)
        @test s.A == A - 2*k
    end

    @testset "getNeighbourStates" begin
        s1 = LdtiSubsystem(1, A, B, C)
        s2 = LdtiSubsystem(2, A, B, C)
        s3 = LdtiSubsystem(3, A, B, C)
        Aij = rand(2,2)
        addNeighbour(s1, s2, Aij)
        addNeighbour(s1, s3, Aij)
        x2 = rand(2)
        x3 = rand(2)
        setInitialState(s2, x2)
        setInitialState(s3, x3)
        x1bar = getNeighbourStates(s1)
        @test all(x1bar .≈ vcat(x2, x3))
    end
end

end