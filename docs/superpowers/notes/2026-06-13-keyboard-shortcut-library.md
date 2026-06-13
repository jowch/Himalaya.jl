# Unified keyboard shortcut library — design (2026-06-13)

Jonathan: "design a unified keyboard shortcut library of gestures so that we are
using the same shortcuts for the same types of actions." This is the durable
source of truth for the gesture vocabulary + architecture. Supersedes the ad-hoc
per-page key handlers (Loupe / Samples / Scoping each hand-roll one today, and
prev/next-sample is bound TWO ways — `[ ]` on Loupe, `,` `.` globally).

## Principle

One gesture = one *type* of action, on every surface. Shortcuts are declared once
in a registry; **handlers, on-screen hints, the legend, and `aria-keyshortcuts`
all derive from that one source** so they cannot drift apart.

## Locked decisions (Jonathan, 2026-06-13)

1. **Prev/next sample = `[` / `]` ONLY.** Retire `,` / `.` as the sample stepper.
2. **Focus arrows = a two-axis model** (Jonathan's refinement): `←`/`→` step the
   exposure, `↑`/`↓` move a previewed ("hovered") assignment candidate, `Esc`
   un-focuses the preview. Nests under `[`/`]` (sample).
3. **Undo = `⌘Z` / `⌘⇧Z`, extended to the Builder** (not just Scoping), each with
   a visible Undo control derived from the registry.

## The gesture vocabulary (by action type)

| Action type | Gesture | Applies on |
|---|---|---|
| Step primary entity — **sample** | `[` / `]` | Loupe, Focus |
| Step sub-entity — **exposure/frame** | `←` / `→` | Loupe (frames), Focus (detector exposures) |
| Step sub-entity — **candidate** | `↑` / `↓` | Focus (assignment-rail preview) |
| Screen exposure | `x` drop · `k` keep · `r` representative | Loupe (active), Samples (selection) |
| Dismiss / up a level | `Esc` (ladder) | innermost popover → modal → clear preview/selection → back to parent |
| Undo / Redo | `⌘Z` / `⌘⇧Z` | Scoping, Builder |
| Find / jump | `⌘K` (alias `/`) | global |
| Grid traversal | arrows · `Home`/`End` · `Enter` | Samples sheet (ARIA grid, focus-scoped) |
| Reorder list item | `Alt+↑` / `Alt+↓` (mirrors `▲▼` buttons) | Builder, Scoping |

### The Esc ladder (consistent "go up a level")
Innermost-first: open popover closes → open modal closes → active preview/selection
clears → back to the parent surface (Loupe→sheet, Focus→sheet). Each surface
supplies only the steps it has; the gesture is always `Esc`.

### Focus two-axis model (detail)
- `[`/`]` = which sample · `←`/`→` = which exposure · `↑`/`↓` = which candidate.
- `↑`/`↓` drive a **page-level `previewedCandidate` index** (NOT DOM focus), so it
  mirrors mouse-hover (lights that candidate's comb/rings on the plot) and does
  not fight the rail's Tab/roving. `Esc` clears the preview first, then backs out.

## Architecture

- **`src/print/shell/shortcuts.ts`** — registry. Each entry:
  `{ id, keys: string[], label, group }`. `keys` use a normalized form
  (`"["`, `"ArrowLeft"`, `"Mod+z"`, `"Mod+Shift+z"`; `Mod` = ⌘ on mac / Ctrl else).
- **`useShortcuts(bindings)`** — one hook binding a window `keydown`, matching the
  registry, honoring `suppressGlobalKeys`, dispatching to the surface's handler.
  Replaces every ad-hoc per-page keydown effect.
- **`<KbdLegend group=…>`** — renders legend rows from the registry for a scope, so
  the printed legend can never disagree with the live bindings.
- `aria-keyshortcuts` on the matching buttons (stepper, loupe nav, undo) is read
  from the registry too.

## Per-surface application

| Action | Samples | Loupe | Focus | Scoping | Builder |
|---|---|---|---|---|---|
| `[ ]` sample | – (grid) | ✓ | ✓ NEW | – | – |
| `← →` sub-entity | grid | frames ✓ | exposures NEW | – | – |
| `↑ ↓` candidate | grid | – | preview NEW | – | – |
| `x/k/r` screen | sel ✓ | ✓ | – | – | – |
| `Esc` ladder | clear sel ✓ | →sheet ✓ | clear preview→sheet NEW | – | – |
| `⌘Z` undo | – | – | – | ✓ + button NEW | ✓ NEW |
| `⌘K` `/` find | ✓ | ✓ | ✓ | ✓ | ✓ |
| `Alt+↑/↓` reorder | – | – | – | NEW (match ▲▼) | ✓ |

## Build order (each: TDD → gate → live-verify → commit)

1. **Registry + `useShortcuts` + `<KbdLegend>`** (foundation, no behavior change yet).
2. **Retire `,`/`.`; make `[`/`]` the global/Loupe sample step** via the library.
3. **Focus**: `[ ]` + `← →` exposures + `↑ ↓` candidate preview + `Esc` ladder
   (the seed bug + the two-axis model). Add stepper `aria-keyshortcuts`.
4. **Scoping**: visible Undo button + `Alt+↑/↓` reorder (match Builder `▲▼`).
5. **Builder**: `⌘Z` undo.
6. **Samples**: legend hint for `x`/`k` when nothing is selected.
7. Migrate Loupe/Samples ad-hoc handlers onto `useShortcuts` (consolidation).

## Progress / DONE (2026-06-13, branch `greenfield-ui-rebuild`)

The unified keyboard shortcut library is **COMPLETE**. Every surface now keys off
the one `shortcuts.ts` registry via `useShortcuts`; no ad-hoc window keydown
handler remains on Loupe/Samples/Focus. Branch STAYS UNMERGED (awaiting Jonathan).

Foundation + earlier steps (prior session): `6a97005` (registry + `useShortcuts`
with the `false`-return DECLINE convention), `329c93c` (Focus two-axis model),
`7eec3f2` (retire `,`/`.`; `[`/`]` is the only sample stepper).

This session (each its own commit, full gate green at each — build + vitest):
- **`c017d5d` — Focus exposure-axis fix (FO-EXPSKIP).** Live-verify of `329c93c`
  surfaced a real defect: `←`/`→` walked the FULL exposure list, but
  `useAutoPickExposure` only holds an *acceptable* exposure, so stepping onto a
  dropped frame bounced the active exposure to the representative (axis looked
  dead when the rep was flanked by dropped frames, e.g. sample 10). Fix: the
  stepper traverses `acceptableExposures(exposures)`. (Loupe `←`/`→` deliberately
  still walks ALL frames — it's the screening surface.)
- **`01a6c8b` — Alt+↑/↓ reorder power-gesture (Scoping + Builder).** New shared
  `useReorderShortcuts` hook resolves the focused row via a
  `data-reorder-row`/`data-reorder-index` attribute and calls the surface's
  existing move (Scoping `moveRow`; Builder routes through the row's own ▲▼
  buttons, reusing their BU-MOVEFOCUS boundary dance). Declines when focus isn't
  in a row.
- **`b5c7299` — Builder ⌘Z/⌘⇧Z undo + visible Undo control (BU-UNDO).**
  Snapshot-based history of the structural recipe edits (reorder/add/remove);
  title keeps native field undo. New `restoreSeriesDraft(slot)` store primitive.
  Hooks live in the top hooks block (above every early return — a misplacement
  here first tripped "rendered fewer hooks"). Visible "↺ Undo" beside Cancel,
  key hint read from the registry.
- **`16de6aa` — Samples X/K cull hint + registry-driven `<KbdLegend>`
  (SA-CULLHINT).** `<KbdLegend ids=[…]>` (shell/) derives caps + descriptions
  from `SHORTCUTS`, delegating appearance to the `KbLegend` ui primitive (now
  with an optional `testId` so legends can coexist). A "Select frames, then
  [X] Drop [K] Keep" hint shows only while nothing is selected.
- **`1f2a8ec` — Loupe + Samples keydown migrated onto `useShortcuts`; Loupe gets
  a registry-driven legend.** Curated to screen verbs (X/K/R) + sample steps
  ([ ]); frame arrows omitted because their registry label says "exposure"
  (LO-TERM keeps the loupe speaking "frame").

All five commits live-verified in a real browser (Playwright MCP, dev-DB copy on
backend :8091 / vite :5182 — NEVER prod :8080), and corpus-culling e2e (12 specs)
re-run green in real Chromium after the keydown migration. Last gate: build +
2942 vitest green at `1f2a8ec`.

Residual / not done (intentional): the cross-surface redo inconsistency — Scoping
still maps `⌘⇧Z` to undo (no redo) while Builder now has real redo. Unifying
Scoping's redo was out of scope for this pass.
