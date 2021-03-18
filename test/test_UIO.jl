@timed_testset "UIO" begin
    @test_throws ErrorException UIO(rand(2,2), rand(2,1), ones(1,2), rand(2,2))
end

