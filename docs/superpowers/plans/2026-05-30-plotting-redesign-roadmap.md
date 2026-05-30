# Plotting Redesign — Implementation Roadmap

The plotting redesign is decomposed into **five sequenced plans**, each producing working, testable software on its own. This index is the map; each plan is self-contained.

## Source documents

- **Design spec (approved):** [2026-05-29-plotting-redesign-design.md](../specs/2026-05-29-plotting-redesign-design.md)
- **Implementation survey + settled architecture:** [2026-05-30-plotting-implementation-survey.md](../specs/2026-05-30-plotting-implementation-survey.md)
- **Mockups (visual/interaction contract):** `docs/redesign-mockups/2026-05-29-focus-plot.html`, `2026-05-29-series-plot.html`
- **Mockup audit (a11y/responsive checklist):** `docs/redesign-mockups/2026-05-30-plotting-mockup-audit.md`

## Settled architecture decisions

1. **One durable assignment per exposure** — retire the multi-group machinery.
2. **Three distinct states** — `indexed` / `form_factor` / `null` as an explicit enum, never inferred.
3. **Cart is durable backend state** — persisted in SQLite, no URL serialization.
4. **Plot layer is a greenfield rebuild** — preserve behaviors (q-link feel, drag/snap), not the custom-SVG code.

## The five plans

| Plan | Scope | Depends on | File |
|---|---|---|---|
| **A** | Backend assignment model — tables, 3-state enum, 3 event kinds, migration, auto-seed, dual-write | — | [plan-a-backend-assignment-model.md](2026-05-30-plan-a-backend-assignment-model.md) |
| **B** | Gauss–Bonnet predictor — kernel (self-contained) + per-candidate `bonnet` flag | kernel: none · flag: A | [plan-b-gauss-bonnet.md](2026-05-30-plan-b-gauss-bonnet.md) |
| **C** | Plot spine — `peakMark()` (retires magenta) + `<PlotSurface>` greenfield | — | [plan-c-plot-spine.md](2026-05-30-plan-c-plot-spine.md) |
| **D** | Focus surface — cart, combs+residual, rings, q-link, preview, custom-index (+B4) | A, B, C | [plan-d-focus-surface.md](2026-05-30-plan-d-focus-surface.md) |
| **E** | Series surface — waterfall/heatmap, migration tracking, phase-strip, reading, member rows, form-factor, clean export | C, A(B2) | [plan-e-series-surface.md](2026-05-30-plan-e-series-surface.md) |

## Dependency graph & suggested order

```
A (backend assignment) ───────────────┐
                                       ├──> D (Focus)
B-kernel ──> B-flag (needs A) ─────────┤
                                       │
C (plot spine) ────────────────────────┴──> E (Series, also needs A's 3-state)
```

**Critical path:** A → C → D. B's kernel and C can start immediately (no backend dep). A gates everything that touches the assignment. C gates every plot feature in D and E. Recommended landing order: **A, B-kernel, C** (parallelizable) → **B-flag** (after A) → **D** → **E**. D-10 (legacy retirement) is the last gate before the old `/groups` machinery disappears; hold it until migration confidence is high.

## Cross-plan invariants (carried through every plan)

- **Queue contract:** `client_op_id` minted inside `mutationFn`; negative optimistic ids; reverse-rollback/insertion-replay; per-iteration try/catch.
- **SSE post_state:** assignment events use a distinct `{assignment:{state,members}}` shape (never the curation `{indices}` shape) — or fall back to `invalidate`.
- **Design guard:** new phase-colour files thread `phaseColor()` through JS mark options / `ui/` primitives, never inline `style`/Tailwind colour literals.
- **Physics integrity:** `score()` stays coverage×consistency; Bonnet is a display flag, never scored or persisted.
- **Colour-parameterized marks:** `peakMark()` takes resolved colour (no CSS-var reads) so the canvas export renderer works.

## Open questions for the user / domain expert

- **Gauss–Bonnet `Ia3d = 1.576`** (Plan B) — confirm with `saxs-physics-reviewer` before B-flag ships; restrict to Pn3m↔Im3m (1.279) if unsure.
- **Dual-write vs hard-cutover** (Plan A) — defaulted to dual-write (keeps `main` green); hard-cutover lands A+D together and deletes `index_groups` outright.
- **S5 q₁ source** (Plan E) — client-derive `min(effective_peaks.q)` (default) vs add `q1` to the snapshot if too costly.
