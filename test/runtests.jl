using Subsystems
using Test
using TimerOutputs

macro timed_testset(name::String, block)
    # copied from https://github.com/KristofferC/Tensors.jl/blob/master/test/runtests.jl#L8
    return quote
        @timeit "$($(esc(name)))" begin
            @testset "$($(esc(name)))" begin
                $(esc(block))
            end
        end
    end
end

@testset "Subsystems.jl" begin
    include("test_UIO.jl")
end
