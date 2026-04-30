using Test, FileIO, TiffImages, ImageCore, ColorTypes
using ImageCore.FixedPointNumbers: Q0f31

@testset "load_and_lognormalize" begin
    # Build a 3x3 with three regions:
    #   (1,1) dead pixel  — negative raw count, must clip to pure black
    #   count=50 cells   — faint background, must collapse to black via p5 low clip
    #   count=100 cells  — mid-range, must fall strictly between black and white
    #   count=1e6 cells  — direct beam / bright ring, must clip to white via p99
    # The Q0f31 bit pattern survives the TIFF round trip; the production path
    # recovers raw counts via `reinterpret(Int32, channelview(raw))`.
    counts = Int32[
        -1     50      100;
         100   100     1_000_000;
         50    100     1_000_000;
    ]
    img = Gray.(reinterpret.(Q0f31, counts))
    path = tempname() * ".tiff"
    save(path, img)
    try
        out = HimalayaUI.load_and_lognormalize(path)
        arr = Float32.(channelview(out))

        @test arr[1, 1] == 0f0           # dead pixel → black (same as background)
        @test arr[1, 2] == 0f0           # faint background (count=50) → black
        @test arr[3, 1] == 0f0           # faint background (count=50) → black
        @test arr[2, 3] == 1f0           # bright (count=1e6) → white (clipped at hi)
        @test arr[3, 3] == 1f0           # bright (count=1e6) → white
        @test 0f0 < arr[1, 3] < 1f0      # mid (count=100) sits strictly between
    finally
        rm(path; force=true)
    end
end
