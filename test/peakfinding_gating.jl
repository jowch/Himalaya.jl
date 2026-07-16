using DelimitedFiles

# Unit tests for the noise-gating kwargs added to findpeaks:
#   q_trim_high       — high-q trim before kneedle
#   prom_ratio_floor  — relative-prominence backstop
# See docs/peak-finding.md.

function _load(name)
    A = readdlm(joinpath(@__DIR__, "data", name))
    A[:, 1], A[:, 2], A[:, 3]
end

@testset "noise gating: disabled kwargs reproduce pre-change behaviour" begin
    # With both gates off, findpeaks reduces to kneedle-only — the
    # form-factor trace returns its known 2 spurious peaks.
    q, I, σ = _load("form-factor_tot.dat")
    pk = Himalaya.findpeaks(q, I, σ; q_trim_high = 0.0, prom_ratio_floor = 0.0)
    @test length(pk.q) == 2

    # Real-peak fixtures are unaffected by enabling/disabling the gates
    # (kneedle threshold is far above the ratio floor on these traces).
    for name in ("example_tot.dat", "cubic_tot.dat")
        q, I, σ = _load(name)
        on  = Himalaya.findpeaks(q, I, σ)
        off = Himalaya.findpeaks(q, I, σ; q_trim_high = 0.0, prom_ratio_floor = 0.0)
        @test on.q == off.q
    end
end

@testset "q_trim_high suppresses peaks at the high-q tail" begin
    q, I, σ = _load("form-factor_tot.dat")
    qmin, qmax = extrema(q)

    # Both spurious peaks sit at qfrac ≥ 0.99 in this fixture; a 0.05 trim
    # excludes them entirely.
    pk_default = Himalaya.findpeaks(q, I, σ)
    @test all(pk_default.q .<= qmin + 0.95 * (qmax - qmin))

    # An aggressive trim is even stricter.
    pk_half = Himalaya.findpeaks(q, I, σ; q_trim_high = 0.5, prom_ratio_floor = 0.0)
    @test all(pk_half.q .<= qmin + 0.5 * (qmax - qmin))
end

@testset "prom_ratio_floor backstop rejects pure-noise traces" begin
    q, I, σ = _load("form-factor_tot.dat")

    # With q_trim_high disabled, the backstop alone must still suppress
    # both spurious peaks.
    pk = Himalaya.findpeaks(q, I, σ; q_trim_high = 0.0, prom_ratio_floor = 30.0)
    @test isempty(pk.q)

    # Extreme floor returns zero peaks on every fixture (sanity).
    for name in ("example_tot.dat", "cubic_tot.dat", "form-factor_tot.dat")
        q, I, σ = _load(name)
        pk = Himalaya.findpeaks(q, I, σ; prom_ratio_floor = 1e6)
        @test isempty(pk.q)
    end
end

@testset "prom_ratio_floor is skipped on low-candidate-count traces" begin
    # A synthetic trace with one Lorentzian peak on a flat baseline produces
    # only a handful of candidates. The candidate-prominence median is then
    # dominated by the peak itself, so applying `30 × median` would suppress
    # the very peak we're trying to find. The carve-out at
    # RATIO_FLOOR_MIN_CANDIDATES = 20 keeps the gate inert in this regime.
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σ = fill(1.0, n)
    γ = 0.003
    q0 = 0.2
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- q0) ./ γ).^2)

    pk = Himalaya.findpeaks(q, I, σ; prom_ratio_floor = 30.0)
    @test length(pk.q) == 1
    @test abs(pk.q[1] - q0) < 2 * step(range(0.005, 0.4; length = n))
end

@testset "prom_floor and prom_ratio_floor compose as a max" begin
    q, I, σ = _load("example_tot.dat")

    # A very high manual prom_floor dominates the ratio floor.
    pk_high = Himalaya.findpeaks(q, I, σ; prom_floor = 1e6, prom_ratio_floor = 0.0)
    @test isempty(pk_high.q)

    # A modest manual prom_floor below the ratio floor does NOT bypass it
    # — the ratio floor still applies on a noise-only trace.
    qf, If, σf = _load("form-factor_tot.dat")
    pk_compose = Himalaya.findpeaks(qf, If, σf;
                                    q_trim_high = 0.0,
                                    prom_floor = 0.1,
                                    prom_ratio_floor = 30.0)
    @test isempty(pk_compose.q)
end
