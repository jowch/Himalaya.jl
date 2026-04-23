# Index Scoring Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the ad hoc `score` formula with `coverage × sharpness_consistency`, thread `sharpness` through `indexpeaks`, and add an R² hard gate to the PhasePanel frontend.

**Architecture:** The `Index` struct drops `prom::Real` and gains `sharpness::SparseVector` (same sparsity structure as `peaks`). `indexpeaks` already accepts a second positional array (`proms`) — rename it to `sharpness` and store it per-peak instead of summing. `score` uses the new fields. The pipeline passes `peaks_result.sharpness` instead of `peaks_result.prominence`. The frontend adds a client-side threshold filter on `r_squared`.

**Tech Stack:** Julia (SparseArrays, Statistics already imported), React/TypeScript (Tailwind v4, existing PhasePanel component)

---

## File map

| File | Change |
|------|--------|
| `src/index.jl` | Remove `prom` field; add `sharpness::SparseVector` to `Index`; remove `totalprom`; rename `proms`→`sharpness` in `indexpeaks` internals; store per-peak sharpness vector; rewrite `score` |
| `src/Himalaya.jl` | Remove `totalprom` from exports |
| `test/index.jl` | Update struct assertions; update `score` test; add sharpness-consistency test |
| `packages/HimalayaUI/src/pipeline.jl` | Pass `peaks_result.sharpness` instead of `peaks_result.prominence` |
| `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx` | Add `R2_THRESHOLD` constant; dim alternatives below threshold |
| `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx` | Add test asserting low-R² alternative is dimmed |

---

## Task 1: Update `test/index.jl` to reflect the new `Index` API

**Files:**
- Modify: `test/index.jl`

- [ ] **Step 1: Replace the existing Index test file content**

Open `test/index.jl` and replace the entire file with:

```julia

@testset "Indexing" begin
    @test all(hasfield.(Index, [:basis, :peaks, :sharpness]))
    @test !hasfield(Index, :prom)

    test_peaks    = [1.0, √3, √3 + eps(), 2.0]
    test_sharpness = ones(length(test_peaks))
    indices = indexpeaks(test_peaks, test_sharpness, 0:0.1:2; gaps = false)
    @test length(indices) == 1

    index = only(indices)
    @test phase(index) == Hexagonal
    @test basis(index) == 1
    @test peaks(index) == test_peaks[[1, 2, 4]]
    @test numpeaks(index) == 3

    @test predictpeaks(index) == phaseratios(Hexagonal)

    (; d, R²) = fit(index)
    @test round(d; digits = 3) == 7.255
    @test R² == 1

    # score: coverage × consistency
    # Hexagonal has 14 ratios; peaks at ranks [1,2,3], uniform sharpness → consistency=1
    # coverage = (1/1 + 1/2 + 1/3) / sum(1/r for r in 1:14) ≈ 1.8333 / 3.2515 ≈ 0.564
    @test isapprox(score(index), 0.564; atol = 0.01)
    @test 0.0 ≤ score(index) ≤ 1.0

    # sharpness consistency: heterogeneous peaks score lower
    test_peaks2     = [1.0, √3, 2.0]
    test_sharpness2 = [1.0, 10.0, 1.0]   # mixed wide/sharp
    indices2 = indexpeaks(test_peaks2, test_sharpness2, 0:0.1:2; gaps = false)
    @test length(indices2) == 1
    @test score(only(indices2)) < score(index)

    test_peaks3 = 1:5
    indices3 = indexpeaks(test_peaks3)
    @test length(indices3) == 2

    a, b = indices3
    @test a != b
    @test issubset(a, a)
    @test !issubset(a, b)

    # End-to-end: top-scoring index is a supported phase, score in [0,1]
    A = readdlm(joinpath(@__DIR__, "data", "example_tot.dat"))
    q, I, σ = A[:, 1], A[:, 2], A[:, 3]
    pk = findpeaks(q, I, σ)
    @test !isempty(pk.q)
    indices4 = indexpeaks(pk.q, pk.sharpness)
    @test !isempty(indices4)
    top = first(sort(indices4; by = score, rev = true))
    @test phase(top) in (Pn3m, Im3m, Ia3d, Lamellar, Hexagonal, Square, Fm3m, Fd3m)
    @test 0.0 ≤ score(top) ≤ 1.0
end
```

- [ ] **Step 2: Run the tests to confirm they fail**

```bash
cd /path/to/Himalaya.jl
julia --project=. -e 'using Himalaya, Test, DelimitedFiles; include("test/index.jl")'
```

Expected: error — either `UndefVarError: totalprom` removed from exports, struct field mismatch, or `score` returning values outside `[0,1]`. The point is the tests do not pass yet.

---

## Task 2: Update `Index` struct, remove `totalprom`, update exports

**Files:**
- Modify: `src/index.jl`
- Modify: `src/Himalaya.jl`

- [ ] **Step 1: Replace the `Index` struct and remove `totalprom` in `src/index.jl`**

Find and replace this block (lines 1–31 of `src/index.jl`):

```julia
# OLD — remove this entire block:
struct Index{P<:Phase}
    basis::Real
    peaks::SparseVector{<:Real, <:Integer}
    prom::Real
end

function show(io::IO, index::Index{P}) where P
    idx, xs = findnz(index.peaks)
    peak_str = fill("⋅", length(index.peaks))
    peak_str[idx] = string.(round.(xs; digits = 5))

    print(io, "Index(::$P, $(round(basis(index); digits = 5)), [$(join(peak_str, ' ')...)])")
end

# getters
phase(::Index{P}) where P = P
basis(index::Index) = index.basis
peaks(index::Index) = nonzeros(index.peaks)
numpeaks(index::Index) = nnz(index.peaks)

"""
    totalprom(index)

The total prominence of an index is the sum of the prominences of its peaks.
"""
totalprom(index::Index) = index.prom
```

Replace with:

```julia
struct Index{P<:Phase}
    basis::Real
    peaks::SparseVector{<:Real, <:Integer}
    sharpness::SparseVector{<:Real, <:Integer}
end

function show(io::IO, index::Index{P}) where P
    idx, xs = findnz(index.peaks)
    peak_str = fill("⋅", length(index.peaks))
    peak_str[idx] = string.(round.(xs; digits = 5))

    print(io, "Index(::$P, $(round(basis(index); digits = 5)), [$(join(peak_str, ' ')...)])")
end

# getters
phase(::Index{P}) where P = P
basis(index::Index) = index.basis
peaks(index::Index) = nonzeros(index.peaks)
numpeaks(index::Index) = nnz(index.peaks)
```

- [ ] **Step 2: Remove `totalprom` from the export list in `src/Himalaya.jl`**

Find this line in `src/Himalaya.jl`:

```julia
Index, phase, basis, peaks, numpeaks, totalprom, predictpeaks, missingpeaks,
```

Replace with:

```julia
Index, phase, basis, peaks, numpeaks, predictpeaks, missingpeaks,
```

---

## Task 3: Update `indexpeaks` internals to thread sharpness as a per-peak `SparseVector`

**Files:**
- Modify: `src/index.jl`

The change affects two functions: the public 4-arg dispatcher and the private 5-arg core. Rename `proms` → `sharpness` throughout and change the `Index` constructor to store a `SparseVector` instead of a summed scalar.

- [ ] **Step 1: Update the multi-phase `indexpeaks` dispatcher**

Find:

```julia
function indexpeaks(peaks, proms, domain; removesubsets = true, kwargs...)
    indices = Index[]

    for phase in (Lamellar, Hexagonal, Square, Pn3m, Im3m, Ia3d, Fd3m)
        push!(indices, indexpeaks(phase, peaks, proms, domain; kwargs...)...)
    end
```

Replace with:

```julia
function indexpeaks(peaks, sharpness, domain; removesubsets = true, kwargs...)
    indices = Index[]

    for phase in (Lamellar, Hexagonal, Square, Pn3m, Im3m, Ia3d, Fd3m)
        push!(indices, indexpeaks(phase, peaks, sharpness, domain; kwargs...)...)
    end
```

- [ ] **Step 2: Update the per-phase 4-arg public dispatcher**

Find:

```julia
function indexpeaks(::Type{P}, peaks, proms, domain; gaps = true, tol = 0.0025, requiremin = true) where {P<:Phase}
    indices = indexpeaks(P, peaks, proms, domain, tol, requiremin)
```

Replace with:

```julia
function indexpeaks(::Type{P}, peaks, sharpness, domain; gaps = true, tol = 0.0025, requiremin = true) where {P<:Phase}
    indices = indexpeaks(P, peaks, sharpness, domain, tol, requiremin)
```

- [ ] **Step 3: Update the convenience dispatch methods**

Find:

```julia
indexpeaks(peaks; kwargs...) = indexpeaks(peaks, ones(length(peaks)), peaks; kwargs...)
indexpeaks(peaks, proms; kwargs...) = indexpeaks(peaks, proms, peaks; kwargs...)
indexpeaks(phase::Type{<:Phase}, peaks; kwargs...) = indexpeaks(phase, peaks, ones(length(peaks)), peaks; kwargs...)
indexpeaks(phase::Type{<:Phase}, peaks, proms; kwargs...) = indexpeaks(phase, peaks, proms, peaks; kwargs...)
```

Replace with:

```julia
indexpeaks(peaks; kwargs...) = indexpeaks(peaks, ones(length(peaks)), peaks; kwargs...)
indexpeaks(peaks, sharpness; kwargs...) = indexpeaks(peaks, sharpness, peaks; kwargs...)
indexpeaks(phase::Type{<:Phase}, peaks; kwargs...) = indexpeaks(phase, peaks, ones(length(peaks)), peaks; kwargs...)
indexpeaks(phase::Type{<:Phase}, peaks, sharpness; kwargs...) = indexpeaks(phase, peaks, sharpness, peaks; kwargs...)
```

- [ ] **Step 4: Update the private 5-arg core — signature and `Index` constructor**

Find:

```julia
function indexpeaks(::Type{P}, peaks, proms, domain, tol, requiremin) where {P<:Phase}
```

Replace with:

```julia
function indexpeaks(::Type{P}, peaks, sharpness, domain, tol, requiremin) where {P<:Phase}
```

Then find the `Index` constructor call inside that function (near the bottom):

```julia
        push!(
            indices,
            Index{P}(
                domain[i],
                SparseVector{eltype(peaks),UInt8}(
                    length(ratios), ratio_idx, peaks[peak_idx]
                ),
                sum(proms[peak_idx])
            )
        )
```

Replace with:

```julia
        push!(
            indices,
            Index{P}(
                domain[i],
                SparseVector{eltype(peaks),UInt8}(
                    length(ratios), ratio_idx, peaks[peak_idx]
                ),
                SparseVector{eltype(sharpness),UInt8}(
                    length(ratios), ratio_idx, sharpness[peak_idx]
                )
            )
        )
```

---

## Task 4: Rewrite `score`

**Files:**
- Modify: `src/index.jl`

- [ ] **Step 1: Replace the `score` function**

Find the entire current `score` function:

```julia
function score(index::Index)
    if numpeaks(index) > 1
        _, rsquared = fit(index)
    else
        rsquared = 1
    end

    num_gaps = let
        peak_idx, _ = findnz(index.peaks)
        count(==(0), view(index.peaks, first(peak_idx):last(peak_idx)))
    end
    
    (numpeaks(index) + totalprom(index)) * (1 - 0.25 * num_gaps) * rsquared
end
```

Replace with:

```julia
function score(index::Index{P}) where P
    n = length(phaseratios(P))
    found_idx, _ = findnz(index.peaks)
    denom = sum(1/r for r in 1:n)
    coverage = sum(1/r for r in found_idx) / denom

    _, sharps = findnz(index.sharpness)
    consistency = length(sharps) > 1 ? clamp(1 - std(sharps) / mean(sharps), 0, 1) : 1.0

    coverage * consistency
end
```

---

## Task 5: Run core tests and commit

**Files:** none (verification only)

- [ ] **Step 1: Run the core test suite**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected output: all tests pass, including the updated `@testset "Indexing"`.

- [ ] **Step 2: Commit**

```bash
git add src/index.jl src/Himalaya.jl test/index.jl
git commit -m "feat(score): coverage × sharpness-consistency replaces ad hoc score

- Index struct: prom field → sharpness SparseVector
- indexpeaks: threads sharpness per-peak instead of summing prominence
- score: harmonic-weighted coverage × CV-based sharpness consistency
- totalprom removed from exports"
```

---

## Task 6: Update the HimalayaUI pipeline call site

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl`

- [ ] **Step 1: Update `analyze_exposure!` to pass sharpness**

Find in `packages/HimalayaUI/src/pipeline.jl`:

```julia
    candidates   = Himalaya.indexpeaks(peaks_result.q, peaks_result.prominence)
```

Replace with:

```julia
    candidates   = Himalaya.indexpeaks(peaks_result.q, peaks_result.sharpness)
```

- [ ] **Step 2: Run the HimalayaUI backend tests**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: all tests pass.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl
git commit -m "feat(pipeline): pass sharpness (not prominence) to indexpeaks"
```

---

## Task 7: Frontend — R² hard gate in PhasePanel

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx`
- Modify: `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx`

Indices below `R2_THRESHOLD` remain visible in the Alternatives list but are visually dimmed (reduced opacity) and show a `low R²` label so users can still manually promote them if needed.

- [ ] **Step 1: Add the threshold constant and dim low-R² alternatives in `PhasePanel.tsx`**

Find the constant at the top of the file (after the imports):

```tsx
export interface PhasePanelProps {
```

Insert before it:

```tsx
const R2_THRESHOLD = 0.98;
```

Then find the alternatives list item:

```tsx
            {alternatives.map((ix) => (
              <li
                key={ix.id}
                data-alternative-id={ix.id}
                className="flex items-center gap-2 px-2 py-1 rounded-md hover:bg-bg-hover cursor-default"
                onMouseEnter={() => setHoveredIndex(ix.id)}
                onMouseLeave={() => setHoveredIndex(undefined)}
              >
                <span
                  className="w-2 h-2 rounded-full"
                  style={{ background: phaseColor(ix.phase) }}
                  aria-hidden
                />
                <span className="font-medium">{ix.phase}</span>
                <span className="text-fg-muted text-[13px]">
                  a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                </span>
                <button
                  className="ml-auto text-fg-muted hover:text-success focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                  aria-label={`Add index ${ix.id}`}
                  onClick={() => { if (active) addMember.mutate(ix.id); }}
                  disabled={active === undefined}
                >
                  +
                </button>
              </li>
            ))}
```

Replace with:

```tsx
            {alternatives.map((ix) => {
              const lowR2 = ix.r_squared != null && ix.r_squared < R2_THRESHOLD;
              return (
                <li
                  key={ix.id}
                  data-alternative-id={ix.id}
                  data-low-r2={lowR2 ? "true" : undefined}
                  className={`flex items-center gap-2 px-2 py-1 rounded-md hover:bg-bg-hover cursor-default${lowR2 ? " opacity-40" : ""}`}
                  onMouseEnter={() => setHoveredIndex(ix.id)}
                  onMouseLeave={() => setHoveredIndex(undefined)}
                >
                  <span
                    className="w-2 h-2 rounded-full"
                    style={{ background: phaseColor(ix.phase) }}
                    aria-hidden
                  />
                  <span className="font-medium">{ix.phase}</span>
                  <span className="text-fg-muted text-[13px]">
                    a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                    {lowR2 && <span className="ml-1 text-warning text-[11px]">low R²</span>}
                  </span>
                  <button
                    className="ml-auto text-fg-muted hover:text-success focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                    aria-label={`Add index ${ix.id}`}
                    onClick={() => { if (active) addMember.mutate(ix.id); }}
                    disabled={active === undefined}
                  >
                    +
                  </button>
                </li>
              );
            })}
```

- [ ] **Step 2: Add a test for the R² gate in `PhasePanel.test.tsx`**

Add a new `it` block inside `describe("<PhasePanel> — alternatives", ...)`, after the existing `"renders alternative indices with a + button"` test:

```tsx
  it("dims alternatives with r_squared below 0.98 and shows 'low R²' label", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 0.8,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4, 0.6], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await waitFor(() => expect(screen.getByText("Im3m")).toBeInTheDocument());

    const lowR2Row = document.querySelector<HTMLElement>('[data-low-r2="true"]');
    expect(lowR2Row).not.toBeNull();
    expect(lowR2Row!.className).toMatch(/opacity-40/);
    expect(screen.getByText(/low R²/i)).toBeInTheDocument();

    // High-R² alternative should NOT be dimmed
    const goodRow = document.querySelector<HTMLElement>('[data-alternative-id="10"]');
    expect(goodRow).toBeNull(); // id:10 is in the active group, not alternatives
  });
```

Note: index 10 is in the active group (members: [10]) so it won't appear in alternatives. Index 11 (r_squared: 0.71 < 0.98) is the alternative and should be dimmed.

- [ ] **Step 3: Run frontend tests**

```bash
cd packages/HimalayaUI/frontend
npm test
```

Expected: all tests pass including the new `"dims alternatives with r_squared below 0.98"` test.

- [ ] **Step 4: Run type check and build**

```bash
npm run build
```

Expected: exits 0 with no TypeScript errors.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx \
        packages/HimalayaUI/frontend/test/PhasePanel.test.tsx
git commit -m "feat(ui): dim low-R² alternatives in PhasePanel (threshold 0.98)"
```
