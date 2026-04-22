@testset "findpeaks: single isolated Lorentzian peak" begin
    # A narrow Lorentzian (γ = 0.003) on a flat baseline, uniform σ.
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σ = fill(1.0, n)
    γ = 0.003
    q0 = 0.2
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- q0) ./ γ).^2)

    pk = findpeaks(q, I, σ)
    @test length(pk.indices) == 1
    @test abs(pk.q[1] - q0) < step(range(0.005, 0.4; length = n))
    @test pk.prominence[1] > 0
    @test pk.sharpness[1]  > 0
end
