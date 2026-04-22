@testset "knee" begin
    # L-shaped descent — elbow is where steep drop transitions to shallow tail.
    v = [100.0, 90.0, 80.0, 5.0, 4.0, 3.0, 2.0]
    k = Himalaya.knee(v)
    @test 80.0 <= k <= 100.0

    # Linear descent — no clear elbow; kneedle returns a value near the middle.
    v2 = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0]
    k2 = Himalaya.knee(v2)
    @test 6.0 <= k2 <= 9.0

    # Constant sequence — no elbow exists; conservative: return first value.
    @test Himalaya.knee([5.0, 5.0, 5.0, 5.0]) == 5.0

    # Empty — returns 0.0 so nothing passes a downstream `.>= threshold` test.
    @test Himalaya.knee(Float64[]) == 0.0

    # Single value — returns that value.
    @test Himalaya.knee([7.5]) == 7.5

    # Two values — also degenerate; returns first.
    @test Himalaya.knee([7.0, 3.0]) == 7.0
end
