using Test, HimalayaUI
import HimalayaUI: hash_trace_file, hash_peak_set

@testset "hash_trace_file is stable and content-addressed" begin
    mktempdir() do dir
        p = joinpath(dir, "t.dat")
        write(p, "0.1 1.0 0.1\n0.2 2.0 0.1\n")
        h1 = hash_trace_file(p)
        h2 = hash_trace_file(p)
        @test h1 == h2
        @test length(h1) == 64  # hex SHA-256

        # Different content → different hash
        write(p, "0.1 1.0 0.1\n0.2 99.0 0.1\n")
        h3 = hash_trace_file(p)
        @test h1 != h3
    end
end

@testset "hash_peak_set is order-independent and content-addressed" begin
    a = (q = [0.1, 0.2], sharpness = [1.0, 2.0])
    b = (q = [0.2, 0.1], sharpness = [2.0, 1.0])
    @test hash_peak_set(a) == hash_peak_set(b)

    c = (q = [0.1, 0.2], sharpness = [1.0, 2.0001])
    @test hash_peak_set(a) != hash_peak_set(c)
end
