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

    @testset "updateE" begin
        s = LdtiSubsystem(1, A, B, C)
        s1 = LdtiSubsystem(2, A, B, C)
        k = [0.1 0; 0 0]

        addNeighbour(s, s1, k)
        @test s.E == k

        removeNeighbour(s, 2)
        @test s.E == zeros(s.nx, 0)
        addNeighbour(s, s1, k)
        setNeighbour(s, s1, 2, 0.2*k)
        @test s.E == 0.2*k
        
        s = LdtiSubsystem(1, A, B, C)
        s2 = LdtiSubsystem(3, A, B, C)

        addNeighbour(s, s1, k)
        addNeighbour(s, s2, k)
        @test s.E == repeat(k, 1, 2)
    end

    @testset "getNeighbourStates" begin
        @testset "gets the correct state" begin
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
        @testset "gets correct zero on no neighbours" begin
            s1 = LdtiSubsystem(1, A, B, C)
            x1bar = getNeighbourStates(s1)
            @test x1bar == zero(s1.x)
        end
    end
    
    @testset "updateState" begin
        @testset "implements correct state equation" begin
            s = LdtiSubsystem(1, A, B, C)
            updateState(s, 1.0)
            @test s.xnext == [1.,0.]
            @test s.x == [0., 0.]
        end
    end
end

end