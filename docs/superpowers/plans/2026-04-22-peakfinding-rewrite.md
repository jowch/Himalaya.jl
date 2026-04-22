# Peak-finding rewrite — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `findpeaks` in [src/peakfinding.jl](../../../src/peakfinding.jl) with a CWT-based algorithm whose threshold is a physical SNR in sigmas, eliminating per-trace tuning and rejecting form-factor false positives by construction.

**Architecture:** Two-stage pipeline. Stage 1: CWT with Ricker wavelets at ~12 geometric scales, then ridge linking across scales. Stage 2: weighted Gaussian+linear fit per candidate using the σ(q) column from `_tot.dat`, accept if `A / σ_A ≥ nσ`. New API: `findpeaks(q, I, σ; nσ=5, scales=nothing, min_ridge_length=3)` returns `(; indices, q, snr, width)`.

**Tech Stack:** Julia 1.x, LsqFit.jl (new dep, replaces unused Distributions.jl), SparseArrays (existing), LinearAlgebra (existing).

**Reference spec:** [docs/superpowers/specs/2026-04-22-peakfinding-rewrite-design.md](../specs/2026-04-22-peakfinding-rewrite-design.md)

---

## Background for the implementer (read this first)

You are working in a Julia package called `Himalaya` for indexing SAXS
diffraction patterns. The directory layout you care about:

```
src/
  Himalaya.jl      # module root, exports public API
  peakfinding.jl   # the file you'll rewrite
  index.jl, phase.jl, util.jl  # leave alone
test/
  runtests.jl      # @testset entry point
  index.jl         # existing tests
  Project.toml     # test environment deps
Project.toml       # package deps
scratch/           # gitignored sandbox; contains real test traces
  workflow/data/integrations/*.dat
```

Each `_tot.dat` file is space-separated `q  I(q)  σ(q)` with ~922 rows.
σ is the propagated 1σ uncertainty from azimuthal integration.

To run all tests:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
To run a single test file standalone (faster while iterating):
```bash
julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'
```

The package uses `using Test` and the standard `@testset` macro. Internal
helpers (not exported) are accessed via `Himalaya.helpername` in tests.

---

## File structure

| File | Action | Responsibility |
|---|---|---|
| `Project.toml` | Modify | Swap `Distributions` → `LsqFit` |
| `src/peakfinding.jl` | Full rewrite | Public `findpeaks` + internal helpers `ricker`, `cwt`, `find_ridges`, `fit_peak` |
| `src/Himalaya.jl` | Modify | (No code changes expected — `findpeaks` is already exported. Verify only.) |
| `src/io.jl` | Move | → `examples/io.jl` with header listing extra deps |
| `examples/io.jl` | Create | Moved content + header comment |
| `README.md` | Modify | Update example block |
| `test/data/example_tot.dat` | Create (copy) | Test fixture |
| `test/data/cubic_tot.dat` | Create (copy) | Test fixture |
| `test/data/form-factor_tot.dat` | Create (copy) | Test fixture |
| `test/peakfinding.jl` | Create | Synthetic ground-truth tests |
| `test/peakfinding_real.jl` | Create | Regression tests on labeled real traces |
| `test/runtests.jl` | Modify | Include the two new test files |
| `test/index.jl` | Modify | Add end-to-end integration test using new findpeaks |

---

## Task 1: Project.toml dependency swap

**Files:**
- Modify: `Project.toml`
- Modify: `test/Project.toml`

- [ ] **Step 1.1: Read current Project.toml**

Run: `cat Project.toml`

Expected: shows `Distributions` in `[deps]`, no `LsqFit`.

- [ ] **Step 1.2: Edit Project.toml**

Replace the `[deps]` block with:

```toml
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
```

(`Distributions` removed; `LsqFit` added with its registered UUID.)

- [ ] **Step 1.2b: Update test/Project.toml**

In Julia 1.9+, `DelimitedFiles` is no longer in the default stdlib bundle and
must be explicitly declared. We also need `Statistics` (for `median` in tests
later in Task 8) and `Random` (for `Random.seed!` in Task 10's spike test).
Edit `test/Project.toml` so its `[deps]` reads:

```toml
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
```

- [ ] **Step 1.3: Resolve & instantiate**

Run: `julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`
Expected: no errors; LsqFit and its transitive deps are downloaded.

- [ ] **Step 1.4: Verify the package still loads**

Run: `julia --project=. -e 'using Himalaya; println("ok")'`
Expected: `ok` (no errors). The current `using Distributions` in
`src/peakfinding.jl` will fail — that's the next task's problem.

If the load fails because of `Distributions` import: that's expected,
proceed to Task 2 which removes it.

- [ ] **Step 1.5: Commit**

```bash
git add Project.toml Manifest.toml test/Project.toml
git commit -m "deps: swap Distributions for LsqFit; declare test env deps"
```

---

## Task 2: Move src/io.jl → examples/io.jl

**Files:**
- Delete: `src/io.jl`
- Create: `examples/io.jl`

- [ ] **Step 2.1: Create examples directory and move file**

Run:
```bash
mkdir -p examples
git mv src/io.jl examples/io.jl
```

- [ ] **Step 2.2: Add header comment to examples/io.jl**

Prepend this header (using Edit tool, prepend before the existing `using DelimitedFiles` line):

```julia
# Example utilities for loading SAXS detector images and integration files.
#
# These helpers are NOT part of the Himalaya package — they require extra
# dependencies that are not in Himalaya's Project.toml:
#   - DelimitedFiles  (stdlib, free)
#   - Images, TiffImages
#
# To use, install those packages in your own environment and `include` this file.

```

- [ ] **Step 2.3: Verify Himalaya module no longer references io.jl**

Run: `grep -n "io.jl" src/Himalaya.jl`
Expected: empty (nothing). The file was never `include`d, confirming this
move is safe.

- [ ] **Step 2.4: Commit**

```bash
git add src/io.jl examples/io.jl
git commit -m "refactor: move io.jl to examples/ (was orphaned, deps not in Project.toml)"
```

---

## Task 3: Strip the broken peakfinding.jl down to a stub

We'll rebuild `findpeaks` from scratch via TDD. First, replace the file
with a minimal stub so the package loads cleanly.

**Files:**
- Modify: `src/peakfinding.jl`

- [ ] **Step 3.1: Replace src/peakfinding.jl with stub**

Replace the entire file contents with:

```julia
using LsqFit

# findpeaks is being rebuilt — see docs/superpowers/specs/2026-04-22-peakfinding-rewrite-design.md
# Tasks 4-9 fill in helpers; Task 10 implements the public function.
```

- [ ] **Step 3.2: Verify package loads**

Run: `julia --project=. -e 'using Himalaya; println("ok")'`
Expected: `ok`. (LsqFit is loaded; no `findpeaks` definition yet — that's
fine, no test calls it.)

- [ ] **Step 3.3: Run existing tests to confirm baseline still passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all tests in `test/index.jl` pass. Findpeaks is not yet defined,
but no existing test references it.

- [ ] **Step 3.4: Commit**

```bash
git add src/peakfinding.jl
git commit -m "refactor: strip findpeaks to stub before rewrite"
```

---

## Task 4: Ricker wavelet function

**Files:**
- Modify: `src/peakfinding.jl`
- Create: `test/peakfinding.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 4.1: Wire up the new test file**

Replace `test/runtests.jl` with:

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("peakfinding.jl")
end
```

- [ ] **Step 4.2: Write the failing test for `ricker`**

Create `test/peakfinding.jl` with:

```julia
using LinearAlgebra: norm

@testset "ricker wavelet" begin
    # Ricker (Mexican hat) is the negative normalized 2nd derivative of a Gaussian.
    # At t=0 it has its maximum; at t=±a it crosses zero; for |t|>a it's negative.
    a = 5.0
    @test Himalaya.ricker(0.0, a) > 0
    @test isapprox(Himalaya.ricker(a, a), 0.0; atol = 1e-12)
    @test isapprox(Himalaya.ricker(-a, a), 0.0; atol = 1e-12)
    @test Himalaya.ricker(2a, a) < 0
    # symmetry
    @test Himalaya.ricker(1.3, a) == Himalaya.ricker(-1.3, a)
    # integral over a wide window should be ~0 (zero-mean wavelet)
    ts = -10a:0.01:10a
    @test abs(sum(Himalaya.ricker.(ts, a)) * 0.01) < 1e-6
end
```

- [ ] **Step 4.3: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: `UndefVarError: ricker not defined` (or similar). The "ricker wavelet" testset fails.

- [ ] **Step 4.4: Implement `ricker`**

Append to `src/peakfinding.jl`:

```julia
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
```

- [ ] **Step 4.5: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for "ricker wavelet" testset.

- [ ] **Step 4.6: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl test/runtests.jl
git commit -m "feat(peakfinding): add ricker wavelet"
```

---

## Task 5: CWT computation across scales

**Files:**
- Modify: `src/peakfinding.jl`
- Modify: `test/peakfinding.jl`

- [ ] **Step 5.1: Write the failing test**

Append to `test/peakfinding.jl`:

```julia
@testset "cwt" begin
    # A pure Gaussian peak should produce its largest CWT response at the
    # scale matching its width.
    n = 300
    σ_true = 4.0          # Gaussian std in points
    centre = 150
    y = [exp(-((i - centre)^2) / (2σ_true^2)) for i in 1:n]

    scales = [1.5, 2.5, 4.0, 6.5, 10.5]
    coeffs = Himalaya.cwt(y, scales)

    @test size(coeffs) == (n, length(scales))
    # the row at the peak centre should peak in the column nearest σ_true
    centre_row = coeffs[centre, :]
    best_scale_idx = argmax(centre_row)
    @test scales[best_scale_idx] == 4.0
end
```

- [ ] **Step 5.2: Run test to verify it fails**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: `UndefVarError: cwt not defined`.

- [ ] **Step 5.3: Implement `cwt`**

Append to `src/peakfinding.jl`:

```julia
"""
    cwt(y, scales)

Continuous wavelet transform of `y` using Ricker wavelets at each width in
`scales`. Returns a `length(y) × length(scales)` matrix of coefficients.

Convolution is direct (no FFT). At array edges, the signal is reflected.
"""
function cwt(y, scales)
    n = length(y)
    coeffs = zeros(n, length(scales))
    for (j, a) in enumerate(scales)
        # truncate the kernel beyond ±5a (wavelet ≈ 0)
        m = ceil(Int, 5a)
        kernel = [ricker(t, a) for t in -m:m]
        for i in 1:n
            s = 0.0
            for (k, w) in enumerate(kernel)
                idx = i + (k - m - 1)
                # reflect at edges
                if idx < 1
                    idx = 2 - idx
                elseif idx > n
                    idx = 2n - idx
                end
                s += w * y[idx]
            end
            coeffs[i, j] = s
        end
    end
    coeffs
end
```

- [ ] **Step 5.4: Run test to verify it passes**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: PASS for "cwt" testset.

- [ ] **Step 5.5: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl
git commit -m "feat(peakfinding): add CWT computation with Ricker wavelets"
```

---

## Task 6: Per-scale local maxima

**Files:**
- Modify: `src/peakfinding.jl`
- Modify: `test/peakfinding.jl`

- [ ] **Step 6.1: Write the failing test**

Append to `test/peakfinding.jl`:

```julia
@testset "local_maxima" begin
    # A vector with two clear peaks at indices 3 and 7
    v = [0.0, 1.0, 3.0, 1.0, 0.0, 2.0, 5.0, 2.0, 0.0]
    @test Himalaya.local_maxima(v) == [3, 7]

    # Negative values should not produce maxima (the CWT filtering wants positives only)
    v2 = [-1.0, -0.5, -1.0, 0.0, 1.0, 0.0]
    @test Himalaya.local_maxima(v2) == [5]

    # Endpoints don't count
    v3 = [3.0, 1.0, 2.0, 1.0, 3.0]
    @test Himalaya.local_maxima(v3) == [3]
end
```

- [ ] **Step 6.2: Run test to verify it fails**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: `UndefVarError: local_maxima`.

- [ ] **Step 6.3: Implement `local_maxima`**

Append to `src/peakfinding.jl`:

```julia
"""
    local_maxima(v)

Return indices `i` where `v[i] > v[i-1]` and `v[i] > v[i+1]` and `v[i] > 0`.
Endpoints are never returned.
"""
function local_maxima(v)
    [i for i in 2:(length(v) - 1) if v[i] > 0 && v[i] > v[i-1] && v[i] > v[i+1]]
end
```

- [ ] **Step 6.4: Run test to verify it passes**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: PASS.

- [ ] **Step 6.5: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl
git commit -m "feat(peakfinding): add per-scale local maxima"
```

---

## Task 7: Ridge linking across scales

**Files:**
- Modify: `src/peakfinding.jl`
- Modify: `test/peakfinding.jl`

- [ ] **Step 7.1: Write the failing test**

Append to `test/peakfinding.jl`:

```julia
@testset "find_ridges" begin
    # Synthesize a CWT coefficient matrix where one Gaussian peak at
    # column-position 50 produces ridges across scales.
    n = 100
    scales = [1.5, 2.5, 4.0, 6.5, 10.5]
    σ_true = 4.0
    y = [exp(-((i - 50)^2) / (2σ_true^2)) for i in 1:n]
    coeffs = Himalaya.cwt(y, scales)

    ridges = Himalaya.find_ridges(coeffs, scales; min_ridge_length = 3)

    # Exactly one ridge survives, centred near i=50
    @test length(ridges) == 1
    r = ridges[1]
    @test abs(r.index - 50) <= 1
    @test r.scale == 4.0     # best-matched scale
    @test r.length >= 3
end

@testset "find_ridges rejects single-scale noise" begin
    # A noisy CWT with maxima only at the smallest scale should produce no ridges.
    n = 100
    scales = [1.5, 2.5, 4.0, 6.5, 10.5]
    coeffs = zeros(n, length(scales))
    coeffs[20, 1] = 1.0   # one isolated max at smallest scale only
    coeffs[60, 1] = 1.0   # another isolated max
    ridges = Himalaya.find_ridges(coeffs, scales; min_ridge_length = 3)
    @test isempty(ridges)
end
```

- [ ] **Step 7.2: Run test to verify it fails**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: `UndefVarError: find_ridges`.

- [ ] **Step 7.3: Implement `find_ridges`**

Append to `src/peakfinding.jl`:

```julia
"""
    find_ridges(coeffs, scales; min_ridge_length = 3, max_gap = 1)

Greedy ridge linking across scales. Walks from largest scale to smallest,
extending a ridge if the next-smaller scale has a local maximum within
`±max(1, scale/4)` points of the current ridge head. Allows up to `max_gap`
missing scales before terminating a ridge.

Returns a vector of NamedTuples `(; index, scale, length, max_coeff)`:
- `index`: position of the ridge at the scale where its CWT coefficient is
  maximal (the "best-matched scale")
- `scale`: that best-matched scale value
- `length`: number of scales the ridge spans
- `max_coeff`: maximum CWT coefficient along the ridge

Ridges shorter than `min_ridge_length` scales are discarded.
"""
function find_ridges(coeffs, scales; min_ridge_length = 3, max_gap = 1)
    nscales = length(scales)
    # Per-scale local maxima
    maxima = [local_maxima(coeffs[:, j]) for j in 1:nscales]

    # Each ridge: vector of (scale_idx, position, coeff). Walk largest → smallest.
    ridges = Vector{Vector{Tuple{Int,Int,Float64}}}()
    # Track which maxima at each scale have been consumed
    consumed = [falses(length(maxima[j])) for j in 1:nscales]

    # Seed ridges at the largest scale
    for (k, p) in enumerate(maxima[nscales])
        push!(ridges, [(nscales, p, coeffs[p, nscales])])
        consumed[nscales][k] = true
    end

    # Extend each ridge downward through scales
    for j in (nscales - 1):-1:1
        tol = max(1, ceil(Int, scales[j] / 4))
        for ridge in ridges
            head_scale_idx, head_pos, _ = ridge[end]
            gap = head_scale_idx - j - 1
            if gap > max_gap
                continue   # ridge already terminated
            end
            # find an unconsumed maximum at scale j within tolerance
            best_k = 0
            best_dist = tol + 1
            for (k, p) in enumerate(maxima[j])
                if consumed[j][k]; continue; end
                d = abs(p - head_pos)
                if d <= tol && d < best_dist
                    best_dist = d
                    best_k = k
                end
            end
            if best_k > 0
                p = maxima[j][best_k]
                push!(ridge, (j, p, coeffs[p, j]))
                consumed[j][best_k] = true
            end
        end
        # Any unconsumed maxima at scale j start new ridges
        for (k, p) in enumerate(maxima[j])
            if !consumed[j][k]
                push!(ridges, [(j, p, coeffs[p, j])])
                consumed[j][k] = true
            end
        end
    end

    # Filter and summarize
    out = NamedTuple{(:index, :scale, :length, :max_coeff),
                     Tuple{Int,Float64,Int,Float64}}[]
    for ridge in ridges
        if length(ridge) < min_ridge_length; continue; end
        # Best-matched scale = where coeff is maximal
        _, idx_in_ridge = findmax(map(t -> t[3], ridge))
        scale_idx, pos, coeff = ridge[idx_in_ridge]
        push!(out, (index = pos, scale = scales[scale_idx],
                    length = length(ridge), max_coeff = coeff))
    end
    out
end
```

- [ ] **Step 7.4: Run tests to verify they pass**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: PASS for both new ridge testsets.

- [ ] **Step 7.5: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl
git commit -m "feat(peakfinding): add greedy ridge linking across scales"
```

---

## Task 8: Stage 2 Gaussian fit + SNR

**Files:**
- Modify: `src/peakfinding.jl`
- Modify: `test/peakfinding.jl`

- [ ] **Step 8.1: Write the failing test**

Append to `test/peakfinding.jl`:

```julia
using Statistics: median

@testset "fit_peak" begin
    # Synthesize a clean Gaussian + linear baseline + Gaussian noise.
    # Then assert fit_peak recovers amplitude and a high SNR.
    n = 60
    q = collect(range(0.0, 1.0; length = n))
    A_true = 100.0
    q0_true = 0.5
    σw_true = 0.05
    baseline = 5.0 .+ 2.0 .* q
    peak = A_true .* exp.(-((q .- q0_true).^2) ./ (2σw_true^2))
    σ_arr = fill(1.0, n)
    I = peak .+ baseline   # noiseless to make assertions tight
    centre_idx = argmin(abs.(q .- q0_true))

    fit = Himalaya.fit_peak(q, I, σ_arr, centre_idx, σw_true / step(q))

    @test fit !== nothing
    @test isapprox(fit.A, A_true; rtol = 0.01)
    @test isapprox(fit.q0, q0_true; atol = step(q))
    @test isapprox(fit.σw, σw_true; rtol = 0.05)
    # SNR should be very high for a noiseless peak with σ=1
    @test fit.snr > 50
end

@testset "fit_peak returns nothing on failure" begin
    # Pathological window: flat data, no peak. LsqFit may not converge or
    # may converge to A≈0 with huge σ_A. fit_peak should return nothing
    # in either case (we treat both as "no detectable peak here").
    n = 60
    q = collect(range(0.0, 1.0; length = n))
    I = fill(10.0, n)
    σ_arr = fill(1.0, n)
    fit = Himalaya.fit_peak(q, I, σ_arr, 30, 5.0)
    @test fit === nothing || fit.snr < 1
end
```

- [ ] **Step 8.2: Run tests to verify they fail**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: `UndefVarError: fit_peak`.

- [ ] **Step 8.3: Implement `fit_peak`**

Append to `src/peakfinding.jl`:

```julia
"""
    fit_peak(q, I, σ, centre_idx, width_pts)

Locally fit a Gaussian + linear baseline to the trace `(q, I)` with
per-point uncertainties `σ`, in a window of `±2.5 · width_pts` points
around `centre_idx`. Uses weighted Levenberg-Marquardt (`LsqFit.curve_fit`)
with weights `1/σᵢ²`.

Returns a NamedTuple `(; A, q0, σw, snr)` on success:
- `A`: fitted peak amplitude (above baseline)
- `q0`: fitted peak centre in q-units
- `σw`: fitted Gaussian σ in q-units
- `snr`: `A / σ_A` where `σ_A` is the standard error on `A` from the
  covariance matrix

Returns `nothing` on convergence failure or non-physical fits (A ≤ 0,
σw ≤ 0, or window too small).
"""
function fit_peak(q, I, σ, centre_idx, width_pts)
    n = length(q)
    half = max(5, ceil(Int, 2.5 * width_pts))
    lo = max(1, centre_idx - half)
    hi = min(n, centre_idx + half)
    if hi - lo < 8
        return nothing
    end

    qw = q[lo:hi]
    Iw = I[lo:hi]
    σw_arr = σ[lo:hi]

    Δq = q[2] - q[1]
    A0 = I[centre_idx] - median(Iw)
    q0_init = q[centre_idx]
    σw_init = max(width_pts * Δq / 2.355, Δq)
    b0 = (Iw[1] + Iw[end]) / 2
    m0 = (Iw[end] - Iw[1]) / (qw[end] - qw[1])

    model(x, p) = p[1] .* exp.(-((x .- p[2]).^2) ./ (2 * p[3]^2)) .+ p[4] .+ p[5] .* x
    p0 = [A0, q0_init, σw_init, b0, m0]
    weights = 1.0 ./ (σw_arr .^ 2)

    fit = try
        curve_fit(model, qw, Iw, weights, p0)
    catch
        return nothing
    end

    fit.converged || return nothing
    A_fit, q0_fit, σw_fit = fit.param[1], fit.param[2], fit.param[3]
    if A_fit <= 0 || σw_fit <= 0
        return nothing
    end

    σ_A = try
        stderror(fit)[1]
    catch
        return nothing
    end
    if !isfinite(σ_A) || σ_A <= 0
        return nothing
    end

    (A = A_fit, q0 = q0_fit, σw = σw_fit, snr = A_fit / σ_A)
end
```

- [ ] **Step 8.4: Run tests to verify they pass**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: PASS for both fit_peak testsets.

- [ ] **Step 8.5: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl
git commit -m "feat(peakfinding): add weighted Gaussian+linear local fit with SNR"
```

---

## Task 9: Top-level `findpeaks` tying it all together

**Files:**
- Modify: `src/peakfinding.jl`
- Modify: `test/peakfinding.jl`

- [ ] **Step 9.1: Write the failing test (single isolated peak)**

Append to `test/peakfinding.jl`:

```julia
@testset "findpeaks: single isolated peak" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    A = 100.0
    q0 = 0.2
    σw = 0.003
    σ_arr = fill(5.0, n)   # ⇒ A/σ = 20 at the peak
    I = 50.0 .+ A .* exp.(-((q .- q0).^2) ./ (2σw^2))

    pk = findpeaks(q, I, σ_arr)
    @test length(pk.q) == 1
    @test abs(pk.q[1] - q0) < step(q)
    @test pk.snr[1] > 10
end
```

- [ ] **Step 9.2: Run test to verify it fails**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: `UndefVarError: findpeaks` or `MethodError`.

- [ ] **Step 9.3: Implement `findpeaks`**

Append to `src/peakfinding.jl`:

```julia
const DEFAULT_SCALES = [2.0, 2.6, 3.5, 4.5, 6.0, 8.0, 10.0, 13.0, 17.0, 22.0, 28.0, 35.0]

"""
    findpeaks(q, I, σ; nσ = 5, scales = nothing, min_ridge_length = 3)

Detect Bragg peaks in a SAXS integration trace.

# Arguments
- `q`, `I`, `σ`: equal-length vectors from a `_tot.dat` file. `σ` is the
  per-point intensity uncertainty (third column).
- `nσ`: SNR threshold in sigmas. A candidate is accepted iff `A / σ_A ≥ nσ`,
  where `A` is the fitted peak amplitude and `σ_A` its standard error.
  Default 5.
- `scales`: vector of CWT widths in *points*. `nothing` ⇒ a default
  geometric series from 2 to 35 points covering typical Bragg widths.
- `min_ridge_length`: a CWT ridge must persist across this many adjacent
  scales to be considered. Default 3.

# Returns
A NamedTuple `(; indices, q, snr, width)` of equal-length vectors:
- `indices`: nearest grid indices into the input arrays
- `q`: sub-pixel peak positions from the local Gaussian fit
- `snr`: physical SNR (`A / σ_A`)
- `width`: fitted FWHM in q-units

Peaks are returned sorted by ascending `q`.
"""
function findpeaks(q, I, σ; nσ = 5, scales = nothing, min_ridge_length = 3)
    @assert length(q) == length(I) == length(σ) "q, I, σ must be the same length"
    sc = scales === nothing ? DEFAULT_SCALES : collect(float.(scales))

    # Stage 1: CWT on log10(I) (compress dynamic range so high-q peaks register)
    y = log10.(max.(I, eps()))
    coeffs = cwt(y, sc)
    ridges = find_ridges(coeffs, sc; min_ridge_length = min_ridge_length)

    # Stage 2: local fit + σ-based SNR validation
    indices = Int[]
    qs = Float64[]
    snrs = Float64[]
    widths = Float64[]

    for r in ridges
        width_pts = r.scale
        fit = fit_peak(q, I, σ, r.index, width_pts)
        fit === nothing && continue
        fit.snr < nσ && continue

        # Snap fitted q0 back to nearest grid index
        idx = argmin(abs.(q .- fit.q0))
        push!(indices, idx)
        push!(qs, fit.q0)
        push!(snrs, fit.snr)
        push!(widths, 2.355 * fit.σw)   # FWHM in q-units
    end

    # Sort by q
    perm = sortperm(qs)
    (indices = indices[perm], q = qs[perm], snr = snrs[perm], width = widths[perm])
end
```

- [ ] **Step 9.4: Run test to verify it passes**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: PASS for "findpeaks: single isolated peak" testset.

- [ ] **Step 9.5: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl
git commit -m "feat(peakfinding): implement public findpeaks(q, I, σ)"
```

---

## Task 10: Synthetic test scenarios from the spec

Six more tests covering the full scenario list from the spec.

**Files:**
- Modify: `test/peakfinding.jl`

- [ ] **Step 10.1: Write the multi-peak test**

Append to `test/peakfinding.jl`:

```julia
@testset "findpeaks: 3 peaks at Pn3m ratios on power-law background" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    # Pn3m first three normalized ratios: √2, √3, √4 → divide by √2 → 1, √(3/2), √2
    base_q = 0.05
    peak_qs = base_q .* [1.0, sqrt(3/2), sqrt(2)]
    σw = 0.003
    background = 1000.0 .* (q ./ q[1]).^(-2)   # power law
    I = copy(background)
    for q0 in peak_qs
        I .+= 100.0 .* exp.(-((q .- q0).^2) ./ (2σw^2))
    end
    σ_arr = fill(10.0, n)   # A/σ ≈ 10 at peak

    pk = findpeaks(q, I, σ_arr)
    @test length(pk.q) == 3
    for q0 in peak_qs
        @test any(abs.(pk.q .- q0) .< step(q))
    end
end
```

- [ ] **Step 10.2: Write below-threshold rejection test**

Append:

```julia
@testset "findpeaks: below-threshold peak rejected" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σw = 0.003
    σ_arr = fill(50.0, n)
    I = 100.0 .+ 100.0 .* exp.(-((q .- 0.2).^2) ./ (2σw^2))   # A/σ ≈ 2

    pk = findpeaks(q, I, σ_arr; nσ = 5)
    @test isempty(pk.q)
end
```

- [ ] **Step 10.3: Write above-threshold detection test**

Append:

```julia
@testset "findpeaks: just-above-threshold peak detected" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σw = 0.003
    σ_arr = fill(15.0, n)
    I = 100.0 .+ 100.0 .* exp.(-((q .- 0.2).^2) ./ (2σw^2))   # A/σ ≈ 6 (margin over nσ=5)

    pk = findpeaks(q, I, σ_arr; nσ = 5)
    @test length(pk.q) == 1
end
```

- [ ] **Step 10.4: Write form-factor false-positive test**

Append:

```julia
@testset "findpeaks: smooth Lorentzian background produces no peaks" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    # Wide Lorentzian (width ~100 points in q-units)
    width_q = 100 * step(q)
    background = 1.0 ./ (1 .+ ((q .- 0.2) ./ width_q).^2)
    I = 1000.0 .* background
    σ_arr = sqrt.(I)   # roughly Poisson-scaled

    pk = findpeaks(q, I, σ_arr)
    @test isempty(pk.q)
end
```

- [ ] **Step 10.5: Write single-pixel-spike test**

Append:

```julia
@testset "findpeaks: single-pixel spikes do not trigger" begin
    using Random
    Random.seed!(0)
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    I = fill(1000.0, n)
    σ_arr = fill(30.0, n)
    # Three single-pixel outliers, each ~10σ above baseline
    for idx in (200, 500, 800)
        I[idx] = 1300.0
    end

    pk = findpeaks(q, I, σ_arr)
    @test isempty(pk.q)
end
```

- [ ] **Step 10.6: Write close-but-resolved peaks test**

Append:

```julia
@testset "findpeaks: two close-but-resolved peaks" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σw = 0.003
    Δ = 8 * step(q)   # 8 q-grid points apart
    q01, q02 = 0.2 - Δ/2, 0.2 + Δ/2
    σ_arr = fill(10.0, n)
    I = 100.0 .+ 100.0 .* exp.(-((q .- q01).^2) ./ (2σw^2)) .+
                  100.0 .* exp.(-((q .- q02).^2) ./ (2σw^2))

    pk = findpeaks(q, I, σ_arr)
    @test length(pk.q) == 2
end
```

- [ ] **Step 10.7: Run all synthetic tests**

Run: `julia --project=. -e 'using Himalaya, Test; include("test/peakfinding.jl")'`
Expected: ALL synthetic testsets PASS.

If any fail, iterate on the algorithm. Most likely culprits:
- Single-pixel-spikes leaking through → tighten `min_ridge_length` or
  raise scales lower bound to 2.5.
- Form-factor leaking through → confirm log-space transform is active and
  that the Lorentzian width is broader than `max(scales)`.
- Adjacent peaks merging → loosen ridge linking tolerance from `scale/4`
  to `scale/3`.

Document any algorithm tweak in the commit message.

- [ ] **Step 10.8: Commit**

```bash
git add test/peakfinding.jl src/peakfinding.jl
git commit -m "test(peakfinding): add synthetic ground-truth scenarios from spec"
```

---

## Task 11: Real-trace regression tests

**Files:**
- Create: `test/data/example_tot.dat`
- Create: `test/data/cubic_tot.dat`
- Create: `test/data/form-factor_tot.dat`
- Create: `test/peakfinding_real.jl`
- Modify: `test/runtests.jl`

(`test/Project.toml` was already updated in Task 1; `DelimitedFiles` is available.)

- [ ] **Step 11.1: Copy test fixtures**

Run:
```bash
mkdir -p test/data
cp scratch/workflow/data/integrations/example_tot.dat test/data/
cp scratch/workflow/data/integrations/cubic_tot.dat test/data/
cp scratch/workflow/data/integrations/form-factor_tot.dat test/data/
ls -la test/data/
```
Expected: three `.dat` files listed.

- [ ] **Step 11.2: Create test/peakfinding_real.jl**

Create with this exact contents:

```julia
using DelimitedFiles
using Statistics: median

const LABELED_PEAKS = Dict(
    "example_tot.dat" => [
        0.053439, 0.065558, 0.075081, 0.091962, 0.107111,
        0.113171, 0.124858, 0.130052, 0.139142, 0.155590,
    ],
    "cubic_tot.dat" => [
        0.046115, 0.058678, 0.065177, 0.071675, 0.079906,
        0.082939, 0.092037, 0.101568, 0.112832, 0.117164,
        0.121496, 0.124529, 0.130594, 0.137959, 0.143591,
        0.145324, 0.152689, 0.159620, 0.171318, 0.177816,
        0.189513, 0.194712, 0.200777,
    ],
    "form-factor_tot.dat" => Float64[],
)

function load_trace(name)
    A = readdlm(joinpath(@__DIR__, "data", name))
    A[:, 1], A[:, 2], A[:, 3]
end

@testset "findpeaks on $name" for name in keys(LABELED_PEAKS)
    q, I, σ = load_trace(name)
    pk = findpeaks(q, I, σ)
    expected = LABELED_PEAKS[name]
    tol = 2 * median(diff(q))

    if isempty(expected)
        # Diagnostic test: form-factor trace should yield zero peaks.
        @test isempty(pk.q)
    else
        # Each labeled peak should have a returned peak nearby.
        for q_label in expected
            @test any(abs.(pk.q .- q_label) .< tol)
        end
        # And the count should match (no major false positives).
        @test length(pk.q) == length(expected)
    end
end
```

Note: the Test environment needs `DelimitedFiles`. It's a stdlib so it doesn't
need to be added to `test/Project.toml` — but verify in the next step.

- [ ] **Step 11.3: Add the file to runtests.jl**

Modify `test/runtests.jl`:

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("peakfinding.jl")
    include("peakfinding_real.jl")
end
```

- [ ] **Step 11.4: Run regression tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for all three regression testsets.

If `cubic_tot.dat` fails on the strict `length == length(expected)`
assertion (it has 23 dense peaks), do not weaken the test silently. Instead:
- Print the actual `pk.q` values returned vs the labeled set.
- Investigate: are returned peaks plausibly at real maxima the user missed,
  or are they spurious (e.g. ridge linking grabbing a shoulder)?
- Decide between (a) re-labeling, (b) tweaking algorithm parameters,
  (c) loosening to `≥ 0.9 * length(expected)` recovered + `≤ 1.1 * length(expected)` returned.
- Document the decision in the commit message.

- [ ] **Step 11.5: Commit**

```bash
git add test/data test/peakfinding_real.jl test/runtests.jl
git commit -m "test(peakfinding): add labeled-trace regression tests"
```

---

## Task 12: Update README example

**Files:**
- Modify: `README.md`

- [ ] **Step 12.1: Replace the Usage example**

Find this block (currently around lines 24-51 of `README.md`):

```julia
using DelimitedFiles
using Himalaya

load_integration(path) = readdlm(path, ' ', Float64, '\n')

# `integration` contains the values in a tot_file
integration = load_integration("my-high-impact-sample_tot.dat")
qs, logIs = integration[:, 1], log10.(integration[:, 2])

# indices of peaks in the integration array
peak_locations, peak_proms = findpeaks(logIs)
peak_qs = qs[peak_locations]
```

Replace it with:

```julia
using DelimitedFiles
using Himalaya

# `_tot.dat` files are space-separated `q  I(q)  σ(q)` (third column is the
# propagated intensity uncertainty from azimuthal integration).
A = readdlm("my-high-impact-sample_tot.dat", ' ', Float64, '\n')
qs, Is, σs = A[:, 1], A[:, 2], A[:, 3]

# Detect Bragg peaks. Threshold `nσ` is in physical sigmas (default 5);
# the same value works across scatters — no per-file tuning.
pk = findpeaks(qs, Is, σs)        # NamedTuple: (indices, q, snr, width)
peak_qs = pk.q
peak_snrs = pk.snr

# Index the detected peaks against known phases.
indices = indexpeaks(peak_qs, peak_snrs)
```

(Leave the rest of the README — phase-specific indexing example, etc. — unchanged.)

- [ ] **Step 12.2: Verify the example is internally consistent**

Eyeball the change. The downstream `indexpeaks(peak_qs, peak_snrs)` call
takes `peaks, proms` positional args — `peak_snrs` plays the role of
`peak_proms`. This works because `indexpeaks` only sums them in `score`;
SNR is a strictly better quality metric.

- [ ] **Step 12.3: Commit**

```bash
git add README.md
git commit -m "docs: update README example for new findpeaks API"
```

---

## Task 13: Final integration test

**Files:**
- Modify: `test/index.jl`

- [ ] **Step 13.1: Append end-to-end test to test/index.jl**

Append this to the existing `@testset "Indexing" begin … end` block (just
before the final `end`):

```julia
    # End-to-end: load a real trace, run findpeaks, run indexpeaks, assert
    # the top-scoring index is a sensible phase. This catches API contract
    # drift between findpeaks and indexpeaks.
    A = readdlm(joinpath(@__DIR__, "data", "example_tot.dat"))
    q, I, σ = A[:, 1], A[:, 2], A[:, 3]
    pk = findpeaks(q, I, σ)
    @test !isempty(pk.q)
    indices = indexpeaks(pk.q, pk.snr)
    @test !isempty(indices)
    top = first(sort(indices; by = score, rev = true))
    @test phase(top) in (Pn3m, Im3m, Ia3d, Lamellar, Hexagonal, Square, Fd3m)
```

Add `using DelimitedFiles` at the top of `test/index.jl` if it's not already there.

- [ ] **Step 13.2: Run all tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: ALL testsets PASS — index, peakfinding (synthetic), peakfinding_real, integration.

- [ ] **Step 13.3: Commit**

```bash
git add test/index.jl
git commit -m "test: add end-to-end findpeaks → indexpeaks integration test"
```

---

## Task 14: Clean up & version bump

**Files:**
- Modify: `Project.toml`

- [ ] **Step 14.1: Bump patch/minor version**

Read the current `version = "0.4.5"` in `Project.toml`. This is a breaking
API change (`findpeaks` signature changed), so bump to `0.5.0`.

Edit `Project.toml`:
```toml
version = "0.5.0"
```

- [ ] **Step 14.2: Verify final state**

Run all tests one more time:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: full green.

Run: `git status`
Expected: clean working tree (everything committed).

Run: `git log --oneline -20`
Expected: clean commit sequence from this rewrite.

- [ ] **Step 14.3: Commit version bump**

```bash
git add Project.toml
git commit -m "chore: bump to 0.5.0 (breaking findpeaks API change)"
```

---

## Done

The new `findpeaks` is in place with synthetic and real-trace test
coverage. The threshold parameter `nσ` is in physical sigmas — no more
per-trace tuning. Form-factor and single-pixel-noise false positives are
rejected by construction.

Follow-up work (out of scope for this plan, suggested as separate spec):
- Replace magic constants in `score` ([src/index.jl:264](../../../src/index.jl)) with
  an SNR-aware formula now that `findpeaks` returns physical SNRs.
- Add a `return_cwt = false` debugging kwarg to `findpeaks` for
  visualization of the wavelet coefficients.
- Reconsider whether `Square` and `Fm3m` should be in the `indexpeaks`
  loop ([src/index.jl:103](../../../src/index.jl)).
