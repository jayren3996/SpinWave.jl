@testset "local frames" begin
    frame = SpinWave.local_frame([0, 0, 2])
    @test frame.z ≈ [0, 0, 1]
    @test dot(frame.x, frame.y) ≈ 0 atol=1e-12
    @test dot(frame.x, frame.z) ≈ 0 atol=1e-12
    @test dot(frame.y, frame.z) ≈ 0 atol=1e-12
    @test cross(frame.x, frame.y) ≈ frame.z

    tilted = SpinWave.local_frame([1, 2, 3])
    @test norm(tilted.x) ≈ 1
    @test norm(tilted.y) ≈ 1
    @test norm(tilted.z) ≈ 1
    @test cross(tilted.x, tilted.y) ≈ tilted.z

    @test_throws DomainError SpinWave.local_frame([0, 0, 0])
end
