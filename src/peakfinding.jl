using LsqFit

# findpeaks is being rebuilt — see docs/superpowers/specs/2026-04-22-peakfinding-rewrite-design.md
# Tasks 4-9 fill in helpers; Task 10 implements the public function.

"""
    ricker(t, a)

Ricker (Mexican-hat) wavelet of width `a`, evaluated at `t`. Equal to the
negative normalized second derivative of a Gaussian.

    ψ_a(t) = (2 / (√(3a) π^(1/4))) · (1 − t²/a²) · exp(−t² / 2a²)
"""
function ricker(t, a)
    norm_const = 2 / (sqrt(3a) * pi^0.25)
    norm_const * (1 - (t/a)^2) * exp(-(t^2) / (2a^2))
end
