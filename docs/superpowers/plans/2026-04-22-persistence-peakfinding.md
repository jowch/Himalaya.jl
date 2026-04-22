# Persistence peak-finding — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the v1 CWT + Lorentzian-fit implementation of `findpeaks` with a shape-agnostic persistence + sharpness detector, using kneedle thresholding instead of hand-set or label-trained thresholds.

**Architecture:** Three small modules — `persistence.jl` (topological prominence, wrapper over `Peaks.peakproms`), `sharpness.jl` (per-sample curvature by `:savgol` or `:cwt`, sharing `convolve_reflect`), `threshold.jl` (kneedle elbow finder). Top-level `findpeaks` is ~20 lines composing them via an AND-gate on adaptively-thresholded prominence and sharpness.

**Tech Stack:** Julia 1.x. Package deps: `Peaks` (kept, used by `persistence`), `LinearAlgebra` (for `savitzky_golay`'s `I`), `Statistics` (unchanged, used by `index.jl`), `SparseArrays` (unchanged). Dropped: `LsqFit` (no more fitting).

**Reference spec:** [docs/superpowers/specs/2026-04-22-persistence-peakfinding-design.md](../specs/2026-04-22-persistence-peakfinding-design.md)

---

## Background for the implementer

You are on the `peakfinding-rewrite` branch (worktree at `/Users/me/projects/Himalaya-peakfinding`). The v1 CWT+Lorentzian implementation lives in `src/peakfinding.jl` and has comprehensive tests in `test/peakfinding.jl` and `test/peakfinding_real.jl`. We are **replacing** the v1 implementation wholesale — not patching it.

Project layout after this plan:

```
src/
  Himalaya.jl       # module entry (small edits: adjust includes + exports)
  peakfinding.jl    # REWRITTEN — top-level findpeaks (~30 lines)
  persistence.jl    # NEW — wrapper over Peaks.peakproms
  sharpness.jl      # NEW — :savgol + :cwt + helpers
  threshold.jl      # NEW — knee(v) kneedle algorithm
  phase.jl, index.jl, util.jl   # unchanged
test/
  runtests.jl       # include the new module tests
  persistence.jl    # NEW
  sharpness.jl      # NEW
  threshold.jl      # NEW
  peakfinding.jl    # REWRITTEN — synthetic tests against v2 pipeline
  peakfinding_real.jl  # KEPT — re-baseline RECALL_FLOOR / SPURIOUS_CEILING
  index.jl          # KEPT — update end-to-end test for new return shape
```

To run the full test suite:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

To run a single test file standalone:
```bash
julia --project=. -e 'using Himalaya, Test; include("test/<file>.jl")'
```

Tests use `@testset`, `@test`, `@test_throws`. Internal helpers (everything except exported names) are accessed via `Himalaya.<name>`. `findpeaks` is already exported from `src/Himalaya.jl`.

---

## File responsibilities

| File | Responsibility |
|---|---|
| `src/threshold.jl` | `knee(v::AbstractVector) -> Float64` — kneedle elbow finder on a sorted descending sequence |
| `src/persistence.jl` | `persistence(y) -> (; indices, prominence)` — wrapper over `Peaks.peakproms` |
| `src/sharpness.jl` | `sharpness(y; method, scales, m)`, `sharpness_savgol`, `sharpness_cwt`, `ricker`, `convolve_reflect`, `savitzky_golay` (restored verbatim from v0.4.5), `DEFAULT_SCALES` |
| `src/peakfinding.jl` | Top-level `findpeaks(q, I, σ; ...)` and `findpeaks(q, I; ...)` composing the above |

Each module file has one clear responsibility and exposes a minimal surface. The rewrite intentionally avoids cross-imports within the modules — each leaf file depends only on the Julia standard library and explicit deps (e.g., `Peaks`).

---

## Task 1: Strip v1 and drop `LsqFit`

**Files:**
- Modify: `src/peakfinding.jl` (full replacement → stub)
- Modify: `Project.toml` (drop `LsqFit`)

- [ ] **Step 1.1: Stub out `src/peakfinding.jl`**

Replace the entire contents of `src/peakfinding.jl` with:

```julia
# Rewritten in v2: see docs/superpowers/specs/2026-04-22-persistence-peakfinding-design.md
# This file will be filled in by Task 5 after the leaf modules are in place.
```

- [ ] **Step 1.2: Drop `LsqFit` from `Project.toml`**

Replace the `[deps]` block of `Project.toml` with:

```toml
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
```

- [ ] **Step 1.3: Resolve**

Run: `julia --project=. -e 'using Pkg; Pkg.resolve()'`
Expected: succeeds; `LsqFit` and its transitive deps are removed from the resolved set.

- [ ] **Step 1.4: Verify module loads**

Run: `julia --project=. -e 'using Himalaya; println("ok")'`
Expected: `ok`. (No callers of the old `findpeaks` remain in `src/`; the stub is a no-op.)

- [ ] **Step 1.5: Commit**

Note that `Manifest.toml` is gitignored, so we do not stage it.

```bash
git add src/peakfinding.jl Project.toml
git commit -m "refactor: strip v1 peakfinding; drop LsqFit dep

v1 used CWT + Lorentzian fit + sigma-SNR; superseded by v2 (persistence +
sharpness + kneedle). Leaf modules added in subsequent tasks; this commit
sets up the empty slot."
```

---

## Task 2: Kneedle threshold (`threshold.jl`)

TDD: write the failing tests first, then the implementation.

**Files:**
- Create: `src/threshold.jl`
- Create: `test/threshold.jl`
- Modify: `src/Himalaya.jl` (add `include("threshold.jl")`)
- Modify: `test/runtests.jl` (add `include("threshold.jl")`)

- [ ] **Step 2.1: Wire up the test file**

Overwrite `test/runtests.jl` with:

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("threshold.jl")
end
```

The other test files (`persistence.jl`, `sharpness.jl`, `peakfinding.jl`, `peakfinding_real.jl`) are added to `runtests.jl` in later tasks as they become runnable.

- [ ] **Step 2.2: Create the failing test file**

Create `test/threshold.jl`:

```julia
@testset "knee" begin
    # L-shaped descent — elbow is where steep drop transitions to shallow tail.
    # The descent 100→80 is steep; 5→2 is shallow. Kneedle should return a
    # value between 80 and 5 (the elbow region), landing on one of the
    # "signal" values (100, 90, 80).
    v = [100.0, 90.0, 80.0, 5.0, 4.0, 3.0, 2.0]
    k = Himalaya.knee(v)
    @test 80.0 <= k <= 100.0

    # Linear descent — no clear elbow; kneedle still returns a value near
    # the middle by symmetry.
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
```

- [ ] **Step 2.3: Run test to see it fail**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: `UndefVarError: knee not defined in Himalaya` (or similar).

- [ ] **Step 2.4: Create `src/threshold.jl`**

Create `src/threshold.jl` with:

```julia
"""
    knee(v) -> Float64

Find the elbow of a sorted descending sequence `v` using the kneedle
algorithm: the value at the point of maximum perpendicular distance
from the chord connecting the first and last entries, after normalising
both axes to [0, 1] so the distance is scale-invariant.

Used as an adaptive threshold: values ≥ knee(v) are "signal", values
< knee(v) are "noise". No magic numbers, no labels, no per-trace tuning —
the threshold emerges from the shape of the data's own distribution.

Edge cases:
- Empty input → 0.0 (nothing passes).
- Single or two values → v[1].
- Constant sequence → v[1] (no elbow; nothing passes).
- No clear elbow (convex-up curve) → v[1] (conservative; nothing passes).

Precondition: `v` must be sorted descending.
"""
function knee(v)
    n = length(v)
    n == 0 && return 0.0
    n == 1 && return float(v[1])
    n == 2 && return float(v[1])

    y_min   = v[end]
    y_range = v[1] - y_min
    y_range == 0 && return float(v[1])

    best_score = 0.0
    best_idx   = 1
    for i in 2:(n - 1)
        x = (i - 1) / (n - 1)
        y = (v[i] - y_min) / y_range
        score = 1 - x - y                 # > 0 ⇔ point below the (0,1)→(1,0) chord
        if score > best_score
            best_score = score
            best_idx   = i
        end
    end

    float(v[best_idx])
end
```

- [ ] **Step 2.5: Wire into `Himalaya.jl`**

Edit `src/Himalaya.jl`. Find the block of `include(...)` lines near the bottom of the module (they currently are `include("util.jl"); include("phase.jl"); include("peakfinding.jl"); include("index.jl")`). Replace with:

```julia
include("util.jl")
include("phase.jl")
include("threshold.jl")
include("peakfinding.jl")
include("index.jl")
```

(The other `include`s for `persistence.jl` and `sharpness.jl` will be added in Tasks 3 and 4; we do it incrementally so each task leaves the package loadable.)

- [ ] **Step 2.6: Run tests to confirm passing**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all existing `Indexing` tests plus the new `knee` testset pass.

- [ ] **Step 2.7: Commit**

```bash
git add src/threshold.jl src/Himalaya.jl test/threshold.jl test/runtests.jl
git commit -m "feat(peakfinding): add knee() kneedle threshold finder"
```

---

## Task 3: Persistence (`persistence.jl`)

**Files:**
- Create: `src/persistence.jl`
- Create: `test/persistence.jl`
- Modify: `src/Himalaya.jl` (add `include`)
- Modify: `test/runtests.jl` (add `include`)

- [ ] **Step 3.1: Add to `test/runtests.jl`**

Update `test/runtests.jl` to:

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("threshold.jl")
    include("persistence.jl")
end
```

- [ ] **Step 3.2: Write failing tests**

Create `test/persistence.jl`:

```julia
@testset "persistence" begin
    # Single Gaussian peak on a flat baseline (baseline = 1, peak height = 10).
    n = 101
    xs = 1:n
    y = [1.0 + 9.0 * exp(-((x - 50)^2) / (2 * 5^2)) for x in xs]

    p = Himalaya.persistence(y)

    @test length(p.indices) == 1
    @test p.indices[1] == 50
    # Prominence = peak height above flat baseline = 9.0, within numerical tolerance.
    @test isapprox(only(p.prominence), 9.0; atol = 1e-6)

    # Two peaks: heights 10 and 6 at indices 30 and 70, saddle around 3 between them.
    y2 = zeros(Float64, 100)
    y2 .= 0.0
    y2[30] = 10.0; y2[29] = 6.0; y2[31] = 6.0
    y2[70] =  6.0; y2[69] = 3.0; y2[71] = 3.0
    y2[45:55] .= 3.0   # plateau between the peaks

    p2 = Himalaya.persistence(y2)
    @test length(p2.indices) == 2
    @test 30 in p2.indices
    @test 70 in p2.indices

    # Monotonic signal has no local maxima (strict definition excludes endpoints).
    p3 = Himalaya.persistence(collect(1.0:100.0))
    @test isempty(p3.indices)
    @test isempty(p3.prominence)
end
```

- [ ] **Step 3.3: Run to see it fail**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: `UndefVarError: persistence`.

- [ ] **Step 3.4: Create `src/persistence.jl`**

```julia
using Peaks: findmaxima, peakproms

"""
    persistence(y) -> (; indices, prominence)

Compute the topological persistence — equivalently, the prominence — of
every local maximum in the 1D signal `y`.

Intuition: imagine flooding the signal from above. Each local max is
born when the water level passes it and dies when its hill merges with
a taller neighbour (or hits the array boundary). The vertical distance
between birth and death is its prominence.

This is the cleanest mathematical formalisation of criterion (A): a
peak is real if it stands clearly above the local floor. Internally
delegates to `Peaks.peakproms`, which implements the same computation
with careful handling of edges and plateaus.
"""
function persistence(y)
    maxima = findmaxima(y)
    with_proms = peakproms(maxima)
    (indices = with_proms.indices, prominence = with_proms.proms)
end
```

- [ ] **Step 3.5: Wire into `Himalaya.jl`**

Edit `src/Himalaya.jl`'s include block:

```julia
include("util.jl")
include("phase.jl")
include("threshold.jl")
include("persistence.jl")
include("peakfinding.jl")
include("index.jl")
```

- [ ] **Step 3.6: Run tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: `persistence` testset passes.

- [ ] **Step 3.7: Commit**

```bash
git add src/persistence.jl src/Himalaya.jl test/persistence.jl test/runtests.jl
git commit -m "feat(peakfinding): add persistence() over Peaks.peakproms"
```

---

## Task 4: Sharpness (`sharpness.jl`)

The biggest leaf module — bundles five functions and a constant. Still only ~90 lines.

**Files:**
- Create: `src/sharpness.jl`
- Create: `test/sharpness.jl`
- Modify: `src/Himalaya.jl` (add `include`)
- Modify: `test/runtests.jl` (add `include`)

- [ ] **Step 4.1: Add to `test/runtests.jl`**

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("threshold.jl")
    include("persistence.jl")
    include("sharpness.jl")
end
```

- [ ] **Step 4.2: Write failing tests**

Create `test/sharpness.jl`:

```julia
@testset "ricker wavelet" begin
    # Ricker (Mexican hat) is the negative normalized 2nd derivative of a Gaussian.
    a = 5.0
    @test Himalaya.ricker(0.0, a) > 0                      # positive at centre
    @test isapprox(Himalaya.ricker(a, a), 0.0; atol = 1e-12)  # zero-crossing at ±a
    @test Himalaya.ricker(2a, a) < 0                       # negative beyond ±a
    @test Himalaya.ricker(1.3, a) == Himalaya.ricker(-1.3, a)  # symmetric
end

@testset "convolve_reflect" begin
    # Reflecting convolution of a delta with a known kernel should reproduce the
    # kernel at the delta's position.
    y = zeros(21); y[11] = 1.0
    k = [1.0, 2.0, 3.0, 2.0, 1.0]      # length 5 (m = 2)
    out = Himalaya.convolve_reflect(y, k)
    @test out[11] == 3.0                # kernel centre at the delta
    @test out[10] == 2.0
    @test out[12] == 2.0
    @test out[9]  == 1.0
    @test out[13] == 1.0
end

@testset "savitzky_golay reproduces a polynomial's derivative" begin
    # A degree-3 polynomial fitted by SG of order 4 should be recovered
    # exactly within the window's interior. Its 2nd derivative 6·x at
    # interior points should match.
    xs = collect(-10.0:10.0)                   # 21 points, indices 1..21
    y = xs.^3                                  # degree-3 polynomial
    d2 = Himalaya.savitzky_golay(5, 4, y; order = 2)
    # Interior points (indices 6..16 correspond to xs = -5..5); at index i,
    # xs[i] = i - 11 and d²/dx² x^3 = 6x.
    for i in 6:16
        @test isapprox(d2[i], 6 * xs[i]; atol = 1e-6)
    end
end

@testset "sharpness: :savgol picks up a Gaussian peak" begin
    xs = 1:201
    y = [exp(-((x - 101)^2) / (2 * 5^2)) for x in xs]
    s = Himalaya.sharpness(y; method = :savgol, m = 5)
    @test length(s) == length(y)
    # The sharpness is maximised at the peak centre.
    @test argmax(s) == 101
    @test s[101] > 0
    # Far from the peak, sharpness is close to zero.
    @test abs(s[1])   < 0.01
    @test abs(s[end]) < 0.01
end

@testset "sharpness: :cwt picks up a Gaussian peak" begin
    xs = 1:201
    y = [exp(-((x - 101)^2) / (2 * 5^2)) for x in xs]
    s = Himalaya.sharpness(y; method = :cwt)
    @test length(s) == length(y)
    @test argmax(s) == 101
    @test s[101] > 0
end

@testset "sharpness: invalid method throws" begin
    y = rand(100)
    @test_throws ArgumentError Himalaya.sharpness(y; method = :bogus)
end
```

- [ ] **Step 4.3: Run to see failures**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: `UndefVarError` for `ricker`, `convolve_reflect`, `savitzky_golay`, `sharpness`.

- [ ] **Step 4.4: Create `src/sharpness.jl`**

```julia
using LinearAlgebra: I

const DEFAULT_SCALES = [2.0, 3.0, 5.0, 8.0, 13.0]

# ---------------------------------------------------------------------------
# Public: top-level sharpness dispatcher
# ---------------------------------------------------------------------------

"""
    sharpness(y; method = :savgol, scales = DEFAULT_SCALES, m = 5) -> Vector{Float64}

Compute a per-sample sharpness measure for the 1D signal `y`. High value =
sharply curved local feature. Two interchangeable methods:

- `:savgol` (default) — magnitude of the smoothed second derivative via a
  Savitzky-Golay filter (window `2m+1`, polynomial order 4). Single-scale,
  fast, O(n).
- `:cwt` — maximum Ricker (Mexican-hat) wavelet response across `scales`,
  normalised by 1/a so coefficients are scale-comparable. Multi-scale,
  width-aware. O(n × length(scales)).
"""
function sharpness(y; method = :savgol, scales = DEFAULT_SCALES, m = 5)
    method === :savgol && return sharpness_savgol(y, m)
    method === :cwt    && return sharpness_cwt(y, scales)
    throw(ArgumentError("Unknown sharpness method: $(repr(method)) (use :savgol or :cwt)"))
end

# ---------------------------------------------------------------------------
# Savitzky-Golay second derivative (sharpness_savgol)
# ---------------------------------------------------------------------------

"""
    sharpness_savgol(y, m) -> Vector{Float64}

Return `-d²y/dx²` computed via Savitzky-Golay smoothing (window `2m+1`,
polynomial order 4). A sharp local maximum has strongly negative d²y, so
we flip sign — large positive output = sharply peaked.
"""
sharpness_savgol(y, m) = -savitzky_golay(m, 4, y; order = 2)

# `savitzky_golay(m, n, y; order)` — preserved verbatim from the v0.4.5
# src/peakfinding.jl (commit fddd611). General-purpose SG: any window,
# polynomial order, derivative order.
function savitzky_golay(m, n, y; order = 0)
    num_y = length(y)
    z = -m:m
    J = zeros(2m + 1, n + 1)

    for i = 0:n
        @inbounds J[:, i + 1] .= z .^ i
    end

    # The convolution term matrix
    C = J' \ I(n .+ 1)[:, order .+ 1]   # = pinv(J) picking out the requested order(s)
    Y = zeros(num_y, length(order))

    for i in 1:num_y
        if i <= m
            window_indices = abs.(z .+ i) .+ 1
        elseif i > num_y - m
            window_indices = -abs.(z .+ i .- num_y) .+ num_y
        else
            window_indices = z .+ i
        end

        for j in eachindex(order)
            @inbounds Y[i, j] = C[:, j]' * y[window_indices]
        end
    end

    length(order) == 1 ? Y[:, 1] : Y
end

# ---------------------------------------------------------------------------
# CWT-based sharpness (sharpness_cwt)
# ---------------------------------------------------------------------------

"""
    sharpness_cwt(y, scales) -> Vector{Float64}

For each sample of `y`, return the maximum Ricker-wavelet response across
`scales`, normalised by 1/a so coefficients are scale-comparable. Peak-like
features produce their maximum response at the scale matching their width.
"""
function sharpness_cwt(y, scales)
    n = length(y)
    out = fill(-Inf, n)
    for a in scales
        m = ceil(Int, 5a)
        kernel = [ricker(t, a) for t in -m:m]
        response = convolve_reflect(y, kernel) ./ a
        out .= max.(out, response)
    end
    out
end

"""
    ricker(t, a)

Ricker (Mexican-hat) wavelet of width `a`, evaluated at `t`. Equal to the
negative normalized second derivative of a Gaussian:

    ψ_a(t) = (2 / (√(3a) π^(1/4))) · (1 − t²/a²) · exp(−t² / (2a²))
"""
ricker(t, a) = (2 / (sqrt(3a) * π^0.25)) * (1 - (t/a)^2) * exp(-(t^2) / (2a^2))

# ---------------------------------------------------------------------------
# Shared: direct convolution with edge reflection
# ---------------------------------------------------------------------------

"""
    convolve_reflect(y, kernel) -> Vector{Float64}

Direct 1D convolution of `y` with a symmetric `kernel` of length `2m+1`,
reflecting at edges (no wrap-around, no zero-padding artefacts).
"""
function convolve_reflect(y, kernel)
    n = length(y)
    m = (length(kernel) - 1) ÷ 2
    out = zeros(n)
    for i in 1:n, k in 1:length(kernel)
        idx = i + (k - m - 1)
        idx = idx < 1 ? 2 - idx : idx > n ? 2n - idx : idx
        out[i] += kernel[k] * y[idx]
    end
    out
end
```

- [ ] **Step 4.5: Wire into `Himalaya.jl`**

```julia
include("util.jl")
include("phase.jl")
include("threshold.jl")
include("persistence.jl")
include("sharpness.jl")
include("peakfinding.jl")
include("index.jl")
```

- [ ] **Step 4.6: Run tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all five new sharpness testsets pass, in addition to the prior ones.

- [ ] **Step 4.7: Commit**

```bash
git add src/sharpness.jl src/Himalaya.jl test/sharpness.jl test/runtests.jl
git commit -m "feat(peakfinding): add sharpness() with :savgol and :cwt methods

Restores the v0.4.5 savitzky_golay function verbatim as the engine under
sharpness_savgol; sharpness_cwt uses a small bank of Ricker scales with
1/a normalization. Shared convolve_reflect helper."
```

---

## Task 5: Top-level `findpeaks` (`peakfinding.jl`)

**Files:**
- Modify: `src/peakfinding.jl` (replace stub)
- Create: `test/peakfinding.jl` (replaces the v1 test file; the old content was v1-specific and is being scrapped)
- Modify: `test/runtests.jl` (add `include`)

- [ ] **Step 5.1: Replace v1 synthetic test file**

Overwrite `test/peakfinding.jl` with a single minimal integration test for this task. We add the full stress scenarios in Task 6.

```julia
@testset "findpeaks: single isolated Lorentzian peak" begin
    # A narrow Lorentzian (γ = 0.003) on a flat baseline, uniform σ so
    # normalize_by_σ has no preferential effect.
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
```

- [ ] **Step 5.2: Add to `test/runtests.jl`**

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("threshold.jl")
    include("persistence.jl")
    include("sharpness.jl")
    include("peakfinding.jl")
end
```

- [ ] **Step 5.3: Run to see failure**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: `MethodError` for `findpeaks(q, I, σ)` with the new signature — the stub has no methods.

- [ ] **Step 5.4: Fill in `src/peakfinding.jl`**

Replace the stub's entire contents with:

```julia
"""
    findpeaks(q, I, σ; normalize_by_σ = true,
                        sharpness_method = :savgol,
                        prom_floor  = nothing,
                        sharp_floor = nothing) -> NamedTuple
    findpeaks(q, I;    normalize_by_σ = false, kwargs...) -> NamedTuple

Detect peaks in a 1D signal using topological prominence (criterion A) and
curvature/sharpness (criterion B), each adaptively thresholded via the
kneedle elbow finder. Shape-agnostic; no per-trace tuning.

# Arguments
- `q`, `I`, `σ`: equal-length vectors. `σ` is per-point intensity
  uncertainty (e.g., the third column of a SAXS `_tot.dat` file).
- `normalize_by_σ`: if `true`, internally work on `I ./ σ` so prominence
  is in noise-relative units. Defaults to `true` when σ is given.
- `sharpness_method`: `:savgol` (default, single-scale 2nd derivative) or
  `:cwt` (multi-scale Ricker max response).
- `prom_floor`, `sharp_floor`: optional manual thresholds. When `nothing`
  (default), the kneedle algorithm chooses each from the data's own
  distribution.

# Returns
`(; indices, q, prominence, sharpness)` — four equal-length vectors,
sorted by ascending q.
"""
function findpeaks(q, I, σ; normalize_by_σ    = true,
                              sharpness_method = :savgol,
                              prom_floor       = nothing,
                              sharp_floor      = nothing)
    @assert length(q) == length(I) == length(σ) "q, I, σ must be equal length"

    y = normalize_by_σ ? I ./ σ : I

    cands           = persistence(y)
    sharps_full     = sharpness(y; method = sharpness_method)
    sharps_at_peaks = sharps_full[cands.indices]

    pf = something(prom_floor,  knee(sort(cands.prominence; rev = true)))
    sf = something(sharp_floor, knee(sort(sharps_at_peaks;  rev = true)))

    keep = (cands.prominence .>= pf) .& (sharps_at_peaks .>= sf)

    # Sort surviving candidates by ascending q (== ascending index, since q
    # is monotonic in our convention).
    idx  = cands.indices[keep]
    perm = sortperm(idx)

    (indices    = idx[perm],
     q          = q[idx[perm]],
     prominence = cands.prominence[keep][perm],
     sharpness  = sharps_at_peaks[keep][perm])
end

findpeaks(q, I; normalize_by_σ = false, kwargs...) =
    findpeaks(q, I, ones(length(I)); normalize_by_σ = normalize_by_σ, kwargs...)
```

- [ ] **Step 5.5: Run tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: the new `findpeaks: single isolated Lorentzian peak` testset passes, along with all prior testsets.

- [ ] **Step 5.6: Commit**

```bash
git add src/peakfinding.jl test/peakfinding.jl test/runtests.jl
git commit -m "feat(peakfinding): implement v2 findpeaks via persistence + sharpness AND-gate"
```

---

## Task 6: Full synthetic stress suite

Now add the six scenarios from the spec's Testing section. Adding one new `@testset` at a time with its own commit would be over-granular; we batch them into one task but keep each testset self-contained.

**Files:**
- Modify: `test/peakfinding.jl` (append testsets)

- [ ] **Step 6.1: Append the full scenario suite**

Append these six `@testset` blocks to `test/peakfinding.jl`:

```julia
@testset "findpeaks: 3 peaks at Pn3m ratios on power-law background" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    σ = fill(10.0, n)
    # Pn3m ratios √2, √3, √4 → normalized to 1, √(3/2), √2
    base = 0.05
    peak_qs = base .* [1.0, sqrt(3/2), sqrt(2)]
    γ = 0.003
    I = 1000.0 .* (q ./ q[1]).^(-2)              # power-law background
    for q0 in peak_qs
        I .+= 100.0 ./ (1 .+ ((q .- q0) ./ γ).^2)
    end

    pk = findpeaks(q, I, σ)
    @test length(pk.q) == 3
    for q0 in peak_qs
        @test any(abs.(pk.q .- q0) .< step(range(0.005, 0.4; length = n)))
    end
end

@testset "findpeaks: smooth Lorentzian background produces no peaks" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    Δq = step(range(0.005, 0.4; length = n))
    width_q = 100 * Δq                           # ≫ any Bragg width
    I = 1000.0 ./ (1 .+ ((q .- 0.2) ./ width_q).^2)
    σ = sqrt.(I)                                 # roughly Poisson-scaled

    pk = findpeaks(q, I, σ)
    @test isempty(pk.q)
end

@testset "findpeaks: single-pixel spikes do not trigger :cwt" begin
    using Random
    Random.seed!(0)
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    I = fill(1000.0, n)
    σ = fill(30.0, n)
    for idx in (200, 500, 800)
        I[idx] = 1300.0
    end

    # :cwt is scale-aware; single-pixel spikes produce a large response only
    # at the smallest scale (width 1 ≪ 2), so they fail the AND gate.
    pk = findpeaks(q, I, σ; sharpness_method = :cwt)
    @test isempty(pk.q)
end

@testset "findpeaks: two close-but-resolved peaks" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    Δq = step(range(0.005, 0.4; length = n))
    γ = 0.003
    Δ = 20 * Δq                                  # ~3 HWHM separation
    q01, q02 = 0.2 - Δ/2, 0.2 + Δ/2
    σ = fill(10.0, n)
    I = 100.0 .+ 100.0 ./ (1 .+ ((q .- q01) ./ γ).^2) .+
                  100.0 ./ (1 .+ ((q .- q02) ./ γ).^2)

    pk = findpeaks(q, I, σ)
    @test length(pk.q) == 2
end

@testset "findpeaks: works without σ" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- 0.2) ./ γ).^2)

    pk = findpeaks(q, I)
    @test length(pk.q) >= 1
    @test any(abs.(pk.q .- 0.2) .< 2 * step(range(0.005, 0.4; length = n)))
end

@testset "findpeaks: normalize_by_σ on vs off produces sensible results" begin
    n = 922
    q = collect(range(0.005, 0.4; length = n))
    γ = 0.003
    I = 10.0 .+ 100.0 ./ (1 .+ ((q .- 0.2) ./ γ).^2)
    σ = fill(1.0, n)

    pk_on  = findpeaks(q, I, σ; normalize_by_σ = true)
    pk_off = findpeaks(q, I, σ; normalize_by_σ = false)

    # Uniform σ means both modes should find the same peak.
    @test !isempty(pk_on.q)
    @test !isempty(pk_off.q)
    @test abs(first(pk_on.q) - first(pk_off.q)) < 2 * step(range(0.005, 0.4; length = n))
end
```

- [ ] **Step 6.2: Run tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all new testsets pass.

If any test fails, do not silently weaken the test. Report DONE_WITH_CONCERNS with diagnosis:
- Which scenario failed? Quote the actual `pk.q` returned.
- Is the algorithm genuinely struggling with this case, or is the test expectation off?
- If a test is clearly wrong (e.g., Lorentzians too close to resolve), widen the separation or adjust amplitudes — and document the change in the commit message. Do not loosen `@test` assertions to make numerics pass.

- [ ] **Step 6.3: Commit**

```bash
git add test/peakfinding.jl
git commit -m "test(peakfinding): add v2 synthetic stress suite"
```

---

## Task 7: Re-baseline real-trace regression

The `test/peakfinding_real.jl` from v1 has `RECALL_FLOOR` / `SPURIOUS_CEILING` tuned to the v1 algorithm. v2 will produce different numbers; we run, observe, and update the floors to reflect v2's actual behaviour.

**Files:**
- Modify: `test/runtests.jl` (add `include`)
- Modify: `test/peakfinding_real.jl` (update floors based on observed v2 behaviour)

- [ ] **Step 7.1: Add to `test/runtests.jl`**

```julia
using Himalaya
using Test

@testset "Himalaya" begin
    include("index.jl")
    include("threshold.jl")
    include("persistence.jl")
    include("sharpness.jl")
    include("peakfinding.jl")
    include("peakfinding_real.jl")
end
```

- [ ] **Step 7.2: Run to see what happens**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

The v1 floors were: example `>= 6/10`, cubic `>= 10/23`, form-factor allows `<= 1` spurious. v2 may pass, fail, or over-achieve. Capture the actual numbers.

- [ ] **Step 7.3: Diagnose and re-baseline**

For each of the three traces (`example_tot.dat`, `cubic_tot.dat`, `form-factor_tot.dat`), run a quick one-off to observe the actual v2 behaviour and derive new floors:

```bash
julia --project=. -e '
using Himalaya, DelimitedFiles, Statistics
for name in ("example_tot.dat", "cubic_tot.dat", "form-factor_tot.dat")
    A = readdlm(joinpath("test/data", name))
    q, I, σ = A[:,1], A[:,2], A[:,3]
    pk = findpeaks(q, I, σ)
    labels = get(Dict(
        "example_tot.dat" => [0.053439, 0.065558, 0.075081, 0.091962, 0.107111,
                               0.113171, 0.124858, 0.130052, 0.139142, 0.155590],
        "cubic_tot.dat"   => [0.046115, 0.058678, 0.065177, 0.071675, 0.079906,
                               0.082939, 0.092037, 0.101568, 0.112832, 0.117164,
                               0.121496, 0.124529, 0.130594, 0.137959, 0.143591,
                               0.145324, 0.152689, 0.159620, 0.171318, 0.177816,
                               0.189513, 0.194712, 0.200777],
        "form-factor_tot.dat" => Float64[],
    ), name, Float64[])
    tol = 2 * median(diff(q))
    recov = sum(any(abs.(pk.q .- L) .< tol) for L in labels; init = 0)
    spur  = sum(!any(abs.(P - L) < tol for L in labels) for P in pk.q; init = 0)
    println(name, "  returned=", length(pk.q), "  recovered=", recov,
            "/", length(labels), "  spurious=", spur)
end'
```

Capture the output. Update the `RECALL_FLOOR` and `SPURIOUS_CEILING` dicts in `test/peakfinding_real.jl` to match observed behaviour. These floors are **regression floors** — they document current behaviour so future changes either preserve or improve them, never regress silently.

Edit `test/peakfinding_real.jl`. Keep the `LABELED_PEAKS` dict unchanged. Update only:

```julia
const RECALL_FLOOR = Dict(
    "example_tot.dat"     => <observed_recall>,
    "cubic_tot.dat"       => <observed_recall>,
    "form-factor_tot.dat" => 0,
)
const SPURIOUS_CEILING = Dict(
    "example_tot.dat"     => <observed_spurious>,
    "cubic_tot.dat"       => <observed_spurious>,
    "form-factor_tot.dat" => <observed_spurious>,
)
```

- [ ] **Step 7.4: Run tests to confirm passing**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all testsets pass.

- [ ] **Step 7.5: Commit**

Include the observed numbers in the commit message for clarity.

```bash
git add test/peakfinding_real.jl test/runtests.jl
git commit -m "test(peakfinding): re-baseline real-trace floors for v2 algorithm

Observed behaviour on hand-labeled traces with v2 defaults:
  example_tot.dat:    <n>/<10>  recovered, <n> spurious
  cubic_tot.dat:      <n>/<23>  recovered, <n> spurious
  form-factor_tot.dat: 0 peaks, <n> spurious

Floors document current behaviour; tightening is a deliberate choice in
the same commit that improves the algorithm."
```

---

## Task 8: Update README and end-to-end test

**Files:**
- Modify: `README.md`
- Modify: `test/index.jl`

- [ ] **Step 8.1: Update the README usage example**

In `README.md`, find the v1 example block (which uses `pk.snr`). Replace the Julia code block with:

```julia
using DelimitedFiles
using Himalaya

# `_tot.dat` files are space-separated `q  I(q)  σ(q)` (third column is the
# per-point intensity uncertainty from azimuthal integration).
A = readdlm("my-high-impact-sample_tot.dat", ' ', Float64, '\n')
qs, Is, σs = A[:, 1], A[:, 2], A[:, 3]

# Detect peaks using topological prominence + curvature, with adaptive
# (kneedle) thresholds on both. Shape-agnostic; no per-trace tuning.
pk = findpeaks(qs, Is, σs)        # NamedTuple: (indices, q, prominence, sharpness)
peak_qs    = pk.q
peak_proms = pk.prominence

# Index the detected peaks against known phases.
indices = indexpeaks(peak_qs, peak_proms)
```

(Leave the phase-specific example and other sections of the README unchanged.)

- [ ] **Step 8.2: Update end-to-end test**

Open `test/index.jl`. Find the end-to-end block added in v1 (it calls `findpeaks` on `example_tot.dat` and uses `pk.snr`). Replace the relevant lines with:

```julia
    # End-to-end: load a real trace, run findpeaks (v2), run indexpeaks,
    # assert the top-scoring index is in the supported phase set. This
    # catches API contract drift between findpeaks and indexpeaks.
    A = readdlm(joinpath(@__DIR__, "data", "example_tot.dat"))
    q, I, σ = A[:, 1], A[:, 2], A[:, 3]
    pk = findpeaks(q, I, σ)
    @test !isempty(pk.q)
    indices = indexpeaks(pk.q, pk.prominence)
    @test !isempty(indices)
    top = first(sort(indices; by = score, rev = true))
    @test phase(top) in (Pn3m, Im3m, Ia3d, Lamellar, Hexagonal, Square, Fm3m, Fd3m)
```

If `using DelimitedFiles` is not already at the top of `test/index.jl`, add it.

- [ ] **Step 8.3: Run tests**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all tests pass.

- [ ] **Step 8.4: Commit**

```bash
git add README.md test/index.jl
git commit -m "docs,test: update README example + end-to-end test for v2 API"
```

---

## Task 9: Version bump

**Files:**
- Modify: `Project.toml`

- [ ] **Step 9.1: Bump version to 0.5.0**

Find `version = "0.4.5"` in `Project.toml` and change to:

```toml
version = "0.5.0"
```

(Breaking API change from 0.4.5: `findpeaks` signature and return shape; previous `findpeaks(y)` no longer exists.)

- [ ] **Step 9.2: Final verification**

Run all tests one more time:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: full green.

Run: `git status`
Expected: clean working tree.

Run: `git log --oneline | head -12`
Expected: a clean task-by-task commit sequence from this plan.

- [ ] **Step 9.3: Commit version bump**

```bash
git add Project.toml
git commit -m "chore: bump to 0.5.0 (breaking findpeaks API v2)"
```

---

## Done

The v2 pipeline is in place:

- `persistence(y)` gives topological prominence for every local max.
- `sharpness(y; method)` gives a per-sample sharpness measure (`:savgol` or `:cwt`).
- `knee(v)` finds the elbow of a sorted descending sequence — used for adaptive thresholding.
- `findpeaks(q, I, σ)` composes these via an AND-gate, with optional σ-normalisation.

The threshold is **data-adaptive** at every stage — no magic numbers, no labels, no per-trace θ. `nothing`-defaulted `prom_floor` and `sharp_floor` kwargs let callers lock down specific thresholds if needed.

### Follow-up work (out of scope for this plan)

- **Sub-pixel position refinement** (parabolic interpolation around each surviving peak), if callers need better-than-grid position accuracy.
- **`score` function in `src/index.jl`** — its magic constants could be revisited now that peaks come with `prominence` and `sharpness` instead of a fit SNR.
- **Pre-filtering for single-pixel spikes** if `:savgol` trips on them in real data (the `:cwt` path already handles them naturally).
