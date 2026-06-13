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

## Progress / RESUME HERE (2026-06-13, branch `greenfield-ui-rebuild`)

Implemented so far (each its own commit, full gate green at each):
- **Step 1 — foundation**: `6a97005`. `shortcuts.ts` registry + `useShortcuts`
  (incl. the `false`-return DECLINE convention for the Esc ladder).
- **Step 3 — Focus two-axis model**: `329c93c` (built before steps 4-6 per
  Jonathan's ordering). `[`/`]` sample · `←`/`→` exposure · `↑`/`↓` candidate
  preview (page-level `previewIndexId` + a neutral `previewed` ring on Card) ·
  `Esc` ladder (modal → TracePlate disarm → clear preview → back to /samples) ·
  stepper `aria-keyshortcuts`. **LIVE-VERIFY STILL PENDING** — Playwright MCP was
  down when this landed; feel the whole model in the browser first thing.
- **Step 2 — retire `,`/`.`**: `7eec3f2`. `useGlobalShortcuts` is now find/jump
  (`/`, `⌘K`) ONLY; sample step is page-owned (`[`/`]` on Loupe/Focus); the
  contact sheet uses ARIA grid roving. Param + dead route machinery removed.

Key finding that resizes the remaining work — the impeccable polish loop
(`cac7438`) already built most of the old "step 4 Scoping" item:
- Scoping ALREADY has the visible **Undo** control (plate `onUndo` +
  `Step back: <label>`), `⌘Z` (pure StrictMode-safe history), AND
  keyboard-accessible reorder via **▲▼ Move up/Move down** buttons
  (`onMoveBy → moveRow`, with live-region announcements). So **B6 is done**.
- Builder likewise has ▲▼ Move up/Move down buttons (keyboard-accessible) but
  NO `⌘Z` and NO `Alt+↑/↓`.

So the genuinely-remaining keyboard work is just:
- **`Alt+↑/↓` reorder power-gesture** mirroring the ▲▼ buttons on BOTH Scoping
  (`ScopeSampleRow`/`moveRow`) and Builder (`reorderVisual`) — neither has it.
  Bind via the registry's `reorderUp`/`reorderDown` (`Alt+ArrowUp/Down`).
- **Builder `⌘Z`/`⌘⇧Z` undo** + a visible Undo control (registry-driven), to
  match Scoping. Locked decision #3: undo extends to Builder.
- **Samples**: a `<KbdLegend>` hint for `x`/`k` when nothing is selected.
- **Consolidation (step 7)**: migrate Loupe + Samples ad-hoc keydown effects onto
  `useShortcuts`, and add a registry-driven `<KbdLegend>` where each surface shows
  its keys. (`<KbdLegend>` from step 1 is NOT built yet — defer to here.)

Process for the restarted session: Playwright MCP died mid-session and only a
session RESTART repopulates the tool index (see memory
`feedback_playwright_mcp_session_desync`). After restart, verify Playwright is
callable, then: live-verify the Focus model (`329c93c`) FIRST, then do the four
items above each as TDD → gate → live-verify → commit. Live-audit harness:
writable dev-DB copy on backend :8091 / vite :5182 (NEVER prod :8080).
