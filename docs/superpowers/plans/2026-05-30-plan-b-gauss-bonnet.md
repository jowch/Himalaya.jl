# Plan B — Gauss–Bonnet Coexistence Predictor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** When a bicontinuous cubic (Pn3m / Im3m / Ia3d) is in an exposure's assignment, predict the lattice parameter a coexisting cubic *should* have (via the Bonnet lattice ratios), and flag candidates whose measured lattice matches — so the human sees "this Im3m candidate is exactly where Bonnet says it should be sitting next to your Pn3m."

**Architecture:** A pure predictor kernel in core `Himalaya` (`src/bonnet.jl`) computes the predicted coexisting-cubic lattice and a consistency predicate — no scoring change, no schema change (physics integrity: `score()` stays coverage×consistency per `docs/scoring.md`). The HimalayaUI `GET /api/exposures/{id}/indices` route consumes the assignment (Plan A) to attach a per-candidate `bonnet` field. The frontend type gains an optional `bonnet`. The visual badge + optional reorder is Plan D's job; Plan B delivers the data.

**Tech Stack:** Julia (core `Himalaya` + `HimalayaUI`), SQLite/DBInterface, JSON3, stdlib `Test`. TypeScript type addition only on the frontend.

---

## Dependencies & scope

- **Tasks 1–3 (kernel + predicate + unit tests)** are self-contained in core `Himalaya` — **no dependency on Plan A**, can land first/in parallel.
- **Tasks 4–6 (backend wiring + frontend type)** depend on **Plan A** (they read the assignment's member indices). Sequence after A.
- **Out of scope:** the `⭙ Bonnet` badge, the contextual note, and any candidate reordering live in **Plan D** (Focus). Plan B stops at the `bonnet` field on the wire.

## Physics constants — confirmed

The Bonnet lattice-parameter ratios among the three bicontinuous cubics are the canonical lipid-mesophase triple `a_D(Pn3m) : a_P(Im3m) : a_G(Ia3d) = 1.000 : 1.279 : 1.576` (Diamond/Primitive/Gyroid; Hyde; Squires/Templer/Seddon). **Both `1.279` and `1.576` are confirmed by `saxs-physics-reviewer` (2026-05-30) and can ship** — they track the NGC constants already hard-coded in `index.jl:311-313` (A₀ = 1.919/2.345/3.091), which derive from the same IPMS geometry. The kernel still returns `nothing` for any non-bicontinuous-cubic pairing (correct by construction — the Bonnet relation is cubic↔cubic only).

## File structure

| File | Change | Responsibility |
|---|---|---|
| `src/bonnet.jl` | **Create** | Pure Bonnet predictor kernel + consistency predicate |
| `src/Himalaya.jl` | Modify | `include("bonnet.jl")` + export `bonnet_lattice`, `bonnet_consistent` |
| `test/bonnet.jl` | **Create**; include in `test/runtests.jl` | Kernel unit tests |
| `packages/HimalayaUI/src/routes_analysis.jl` | Modify | `_bonnet_for_index` helper; attach `:bonnet` in the `/indices` serialization (line 118-133) |
| `packages/HimalayaUI/test/test_assignments.jl` | Modify | `_bonnet_for_index` route-logic test |
| `packages/HimalayaUI/frontend/src/api.ts` | Modify | Add `bonnet?: BonnetFlag \| null` to `IndexEntry` (line 268-282) |
| `docs/scoring.md` | Modify | Short note: Bonnet is a display flag, never folded into `score()` |

---

### Task 1: Bonnet predictor kernel (`bonnet_lattice`)

**Files:**
- Create: `src/bonnet.jl`
- Modify: `src/Himalaya.jl`
- Test: `test/bonnet.jl` (create), `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Create `test/bonnet.jl`:

```julia
using Test
using Himalaya
using Himalaya: Pn3m, Im3m, Ia3d, Lamellar, bonnet_lattice

@testset "bonnet_lattice predicts the coexisting cubic" begin
    # a_Im3m ≈ 1.279 · a_Pn3m  (the canonical ratio).
    @test bonnet_lattice(Pn3m, 100.0, Im3m) ≈ 127.9 atol=1e-6
    @test bonnet_lattice(Im3m, 127.9, Pn3m) ≈ 100.0 atol=1e-3

    # Ia3d uses the 1.576 reference (confirm with domain expert before shipping).
    @test bonnet_lattice(Pn3m, 100.0, Ia3d) ≈ 157.6 atol=1e-6

    # Non-bicontinuous phase → nothing (no Bonnet relation defined).
    @test bonnet_lattice(Pn3m, 100.0, Lamellar) === nothing
    @test bonnet_lattice(Lamellar, 100.0, Im3m) === nothing

    # Same phase → identity.
    @test bonnet_lattice(Pn3m, 100.0, Pn3m) ≈ 100.0 atol=1e-9
end
```

- [ ] **Step 2: Register + run to verify it fails**

In `test/runtests.jl`, add `include("bonnet.jl")` alongside the other `include(...)` test lines.
Run: `julia --project=. test/bonnet.jl`
Expected: FAIL — `UndefVarError: bonnet_lattice not defined`.

- [ ] **Step 3: Write the kernel**

Create `src/bonnet.jl`:

```julia
# Gauss–Bonnet lattice-parameter ratios among the three bicontinuous cubic
# phases. Reference scale (Pn3m = 1.0): a_Pn3m : a_Im3m : a_Ia3d.
#
# ⚠ Ia3d = 1.576 must be confirmed by the domain expert (saxs-physics-reviewer)
# before this is wired into the API (Plan B Task 3). Only Pn3m↔Im3m (1.279) is
# committed by the design spec. To restrict to the confirmed pair, delete the
# Ia3d entry — `bonnet_lattice` then returns `nothing` for any Ia3d pairing.
const BONNET_SCALE = Dict{DataType,Float64}(
    Pn3m => 1.000,
    Im3m => 1.279,
    Ia3d => 1.576,
)

"""
    bonnet_lattice(from::Type{<:Phase}, a_from::Real, to::Type{<:Phase}) -> Union{Float64,Nothing}

Predict the lattice parameter a coexisting cubic of phase `to` should have,
given an assigned cubic of phase `from` with lattice `a_from`, via the
Gauss–Bonnet ratios. Returns `nothing` when either phase is not one of the
three bicontinuous cubics (no Bonnet relation defined).
"""
function bonnet_lattice(from::Type{<:Phase}, a_from::Real, to::Type{<:Phase})
    (haskey(BONNET_SCALE, from) && haskey(BONNET_SCALE, to)) || return nothing
    Float64(a_from) * (BONNET_SCALE[to] / BONNET_SCALE[from])
end
```

In `src/Himalaya.jl`, add `include("bonnet.jl")` after the other includes, and add `bonnet_lattice` to the `export` list (or leave unexported and access via `Himalaya.bonnet_lattice` — match the surrounding export convention; `phaseratios`/`fit` patterns show what's exported).

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=. test/bonnet.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bonnet.jl src/Himalaya.jl test/bonnet.jl test/runtests.jl
git commit -m "feat(core): bonnet_lattice — predict coexisting bicontinuous-cubic lattice"
```

---

### Task 2: Consistency predicate (`bonnet_consistent`)

**Files:**
- Modify: `src/bonnet.jl`, `src/Himalaya.jl` (export)
- Test: `test/bonnet.jl`

- [ ] **Step 1: Write the failing test**

Append to `test/bonnet.jl`:

```julia
@testset "bonnet_consistent within relative tolerance" begin
    # Observed Im3m at 128 vs predicted 127.9 from a Pn3m at 100 → consistent.
    @test Himalaya.bonnet_consistent(Pn3m, 100.0, Im3m, 128.0; reltol=0.02)
    # Observed Im3m at 140 → off by ~9% → not consistent at 2%.
    @test !Himalaya.bonnet_consistent(Pn3m, 100.0, Im3m, 140.0; reltol=0.02)
    # Undefined pair → false (no relation to be consistent with).
    @test !Himalaya.bonnet_consistent(Pn3m, 100.0, Lamellar, 100.0; reltol=0.02)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=. test/bonnet.jl`
Expected: FAIL — `bonnet_consistent not defined`.

- [ ] **Step 3: Implement**

Append to `src/bonnet.jl`:

```julia
"""
    bonnet_consistent(from, a_from, to, a_to_observed; reltol=0.02) -> Bool

True when an observed coexisting-cubic lattice `a_to_observed` matches the
Bonnet-predicted lattice within relative tolerance `reltol`. False for phase
pairs with no Bonnet relation.
"""
function bonnet_consistent(from::Type{<:Phase}, a_from::Real,
                           to::Type{<:Phase}, a_to_observed::Real; reltol::Real=0.02)
    pred = bonnet_lattice(from, a_from, to)
    pred === nothing && return false
    abs(a_to_observed - pred) <= reltol * pred
end
```

Export `bonnet_consistent` in `src/Himalaya.jl` (match Task 1's choice).

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=. test/bonnet.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/bonnet.jl src/Himalaya.jl test/bonnet.jl
git commit -m "feat(core): bonnet_consistent predicate"
```

---

### Task 3: `_bonnet_for_index` route helper (depends on Plan A)

**Files:**
- Modify: `packages/HimalayaUI/src/routes_analysis.jl`
- Test: `packages/HimalayaUI/test/test_assignments.jl`

> **Depends on Plan A** — reads `assignment_members` to find the assigned bicontinuous cubic(s). The helper computes, for a given candidate index, whether its `(phase, lattice_d)` is Bonnet-consistent with any *assigned* bicontinuous cubic of a different phase, and what the predicted lattice is.

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_assignments.jl`:

```julia
@testset "_bonnet_for_index flags a coexisting cubic" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # Assigned Pn3m at a=100; a candidate Im3m at a=128 (Bonnet-consistent).
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis, lattice_d) VALUES (10, ?, 'Pn3m', 0.1, 100.0)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis, lattice_d) VALUES (11, ?, 'Im3m', 0.1, 128.0)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis, lattice_d) VALUES (12, ?, 'Im3m', 0.1, 200.0)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
        DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

        # The assigned Pn3m (10) itself: no flag (it's the anchor, same phase).
        @test HimalayaUI._bonnet_for_index(db, e_id, "Pn3m", 100.0, 10) === nothing
        # Candidate Im3m at 128 → consistent, predicted ≈ 127.9.
        b = HimalayaUI._bonnet_for_index(db, e_id, "Im3m", 128.0, 11)
        @test b !== nothing && b[:consistent] == true
        @test b[:predicted_a] ≈ 127.9 atol=1.0
        # Candidate Im3m at 200 → predicted still 127.9, not consistent.
        b2 = HimalayaUI._bonnet_for_index(db, e_id, "Im3m", 200.0, 12)
        @test b2 !== nothing && b2[:consistent] == false
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: FAIL — `_bonnet_for_index not defined`.

- [ ] **Step 3: Implement the helper**

In `packages/HimalayaUI/src/routes_analysis.jl`, add near `_ngc_for_phase` (which ends ~line 109):

```julia
"""
    _bonnet_for_index(db, exposure_id, phase_name, lattice_d, index_id) -> Union{Dict,Nothing}

For a candidate index, return `{predicted_a, consistent}` when the exposure's
assignment contains a bicontinuous cubic of a DIFFERENT phase that predicts this
candidate's lattice via the Bonnet ratio. Returns `nothing` when there is no
applicable anchor (no assigned cubic, same phase, or non-bicontinuous phase).
The anchor index itself is never flagged.
"""
function _bonnet_for_index(db::SQLite.DB, exposure_id::Integer,
                           phase_name::AbstractString, lattice_d, index_id::Integer)
    (lattice_d === nothing || ismissing(lattice_d)) && return nothing
    P = resolve_phase(phase_name)
    P === nothing && return nothing
    a = Float64(lattice_d)
    a > 0 || return nothing

    # Assigned bicontinuous-cubic anchors (phase + lattice), excluding this index.
    anchors = Tables.rowtable(DBInterface.execute(db,
        """SELECT i.phase, i.lattice_d
           FROM assignment_members m JOIN indices i ON i.id = m.index_id
           WHERE m.exposure_id = ? AND i.id != ?
             AND i.lattice_d IS NOT NULL
             AND i.phase IN ('Pn3m', 'Im3m', 'Ia3d')""",
        [Int(exposure_id), Int(index_id)]))

    for anc in anchors
        Pa = resolve_phase(String(anc.phase))
        Pa === nothing && continue
        Pa === P && continue   # same phase: not a coexisting pair
        pred = Himalaya.bonnet_lattice(Pa, Float64(anc.lattice_d), P)
        pred === nothing && continue
        return Dict(:predicted_a => pred,
                    :consistent  => abs(a - pred) <= 0.02 * pred)
    end
    nothing
end
```

- [ ] **Step 4: Attach `:bonnet` in the `/indices` serialization**

In the `GET "/api/exposures/{id}/indices"` map block (lines 118-133), after the `d[:ngc] = ...` line, add:

```julia
            d[:bonnet] = _bonnet_for_index(db, id, String(ix.phase), ix.lattice_d, Int(ix.id))
```

- [ ] **Step 5: Run the test + smoke load**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_assignments.jl`
Expected: PASS.
Run: `julia --project=packages/HimalayaUI -e 'using HimalayaUI'`
Expected: no error.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_analysis.jl packages/HimalayaUI/test/test_assignments.jl
git commit -m "feat(backend): _bonnet_for_index flag on the indices response"
```

---

### Task 4: Frontend `IndexEntry.bonnet` type

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Test: `packages/HimalayaUI/frontend/test/api.test.ts` (or the nearest type/contract test — add a parse assertion)

- [ ] **Step 1: Add the type**

In `packages/HimalayaUI/frontend/src/api.ts`, above `IndexEntry` (line 268), add:

```typescript
export interface BonnetFlag {
  predicted_a: number;
  consistent: boolean;
}
```

Then add the field to `IndexEntry` (after `predicted_q`):

```typescript
  /** Gauss–Bonnet coexistence flag vs the current assignment (null when N/A).
   *  Display-only — never folded into `score`. Rendered as the ⭙ Bonnet badge
   *  in Plan D. */
  bonnet?: BonnetFlag | null;
```

- [ ] **Step 2: Verify the type compiles**

Run: `(cd packages/HimalayaUI/frontend && npx tsc --noEmit)`
Expected: no new type errors.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts
git commit -m "feat(frontend): IndexEntry.bonnet flag type"
```

---

### Task 5: Docs + full-suite green

**Files:**
- Modify: `docs/scoring.md`

- [ ] **Step 1: Note the display-only contract**

In `docs/scoring.md`, add a short paragraph: the Gauss–Bonnet flag (`bonnet` on the indices response) is a **display-and-ranking affordance computed from the assignment**, never folded into `score()` (which stays coverage×consistency). It is recomputed per request, never persisted. Cross-reference `src/bonnet.jl`.

- [ ] **Step 2: Run the core suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: all core tests pass, including the new `bonnet.jl` testsets.

- [ ] **Step 3: Run the HimalayaUI suite (capture once)**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test-planB.out 2>&1`
Then: `grep -E "Test Summary|FAIL|Error|bonnet" /tmp/jl-test-planB.out | head -40`
Expected: all pass; `_bonnet_for_index` test green; no regressions.

- [ ] **Step 4: Commit**

```bash
git add docs/scoring.md
git commit -m "docs(scoring): Gauss–Bonnet is a display flag, never scored"
```

---

## Self-Review

**1. Spec coverage** (survey item B3 + design §"the single smart automation"):
- Predict coexisting cubic lattice from assigned cubic → Task 1. ✓
- Consistency check → Task 2. ✓
- Per-candidate flag on the response (no scoring change, no persistence) → Tasks 3–5. ✓
- Badge/note/reorder explicitly deferred to Plan D. ✓

**2. Placeholder scan:** none. The one physics value flagged for confirmation (`Ia3d = 1.576`) is called out explicitly with a degrade path, not left vague.

**3. Type/name consistency:** `bonnet_lattice`, `bonnet_consistent`, `BONNET_SCALE`, `_bonnet_for_index`, `BonnetFlag`/`{predicted_a, consistent}` are consistent across kernel, route helper, and frontend type. The route returns `Dict(:predicted_a, :consistent)`; the TS type mirrors `{predicted_a, consistent}`.

**Open question for the user / saxs-physics-reviewer:** confirm `a_Ia3d/a_Pn3m = 1.576` before Task 3 ships; restrict `BONNET_SCALE` to Pn3m↔Im3m if not.

---

## Execution Handoff

1. **Subagent-Driven (recommended)** — fresh subagent per task, two-stage review. Tasks 1–2 can run before Plan A; 3–5 after.
2. **Inline Execution** — batch with checkpoints.
