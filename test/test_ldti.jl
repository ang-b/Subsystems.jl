module ldtiSubTests

using LinearAlgebra
using Subsystems
using Test

@testset "ldti Subsystem" begin
    A = Matrix{Float64}(I, 2, 2);
    B = reshape([1., 0.], :, 1)
    C = I;
    s = LdtiSubsystem(1, A, B, C)

    @testset "defines correct system" begin
        @test s.A == A
        @test s.C == I
        @test s.B == B
    end

    @testset "implements correct conservation" begin
        s1 = LdtiSubsystem(2, A, B, C)
        k = [0.1 0; 0 0]
        addNeighbour(s, s1, k, true)
        @test s.A == A - k
        removeNeighbour(s, 2, true)
        @test s.A == A
        addNeighbour(s, s1, k, true)
        setNeighbour(s, s1, 2, 0.2*k, true)
        @test s.A == A - 0.2*k
    end
end

end