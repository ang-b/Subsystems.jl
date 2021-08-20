module LssTests

using Revise
using Subsystems
using BlockArrays
using Test
using Infiltrator

@testset "LSS" begin
    @testset "has correct" begin
        N = 5;
        ne = rand(1:5, N)
        me = rand(0:2, N)
        pe = rand(1:5, N)

        As = [rand(ne[i],ne[i]) for i in 1:N]
        Bs = [rand(ne[i],me[i]) for i in 1:N]
        Cs = [rand(pe[i],ne[i]) for i in 1:N]
        K = [[rand(ne[i],ne[j]) for j in 1:N] for i in 1:N]
        Inc = rand(Bool, N, N);

        subs = [LdtiSubsystem(i, As[i], Bs[i], Cs[i]) for i in 1:N]
        
        Abe = BlockArray{Real}(zeros(sum(ne), sum(ne)), ne, ne)
        Bbe = BlockArray{Real}(zeros(sum(ne), sum(me)), ne, me)
        Cbe = BlockArray{Real}(zeros(sum(pe), sum(ne)), pe, ne)
        
        for i in 1:N
            Abe[Block(i,i)] = As[i]
            Bbe[Block(i,i)] = Bs[i]
            Cbe[Block(i,i)] = Cs[i]
            
            for j in vcat(1:i-1, i+1:N)
                if Inc[i,j] 
                    addNeighbour(subs[i], subs[j], K[i][j], false)
                    view(Abe, Block(i,j)) .= K[i][j]
                end
            end
        end

        @testset "dimensions" begin
            (ng, mg, pg) = _getdimarrays(subs)
            @test ng == ne
            @test mg == me
            @test pg == pe
        end

        @testset "Block matrices" begin
            ms = tomonolith(subs)
            @test Abe == ms.A
            @test Bbe == ms.B
            @test Cbe == ms.C
        end
    end
end

end