# App Shell Unification — Design

> Status: **DRAFT (in review)** · 2026-06-20 · branch `ingestion-redesign`
> Scope: the app **shell, the navigation between surfaces, and one unified keyboard set** — NOT the refined surfaces' layouts.

## 1. Problem

The ingestion redesign added its new experiment-scoped chrome (`ExperimentShell` + `ExperimentTopNav`) *in parallel* to the legacy corpus chrome (`CorpusShell` + `CorpusTopbar`) instead of replacing it. The result is three coexisting navigation grammars:

1. `CorpusShell` / `CorpusTopbar` — wraps `/samples`, `/sample/:id` (Focus), `/series/*`, the `*` catch-all. Nav = **Experiments · Samples · Series**, wordmark → `/samples`, plus a beamtime *filter* chip.
2. `ExperimentShell` / `ExperimentTopNav` — wraps `/experiments/:id/{corpus,config,grouping}`. Nav = **Experiments · Series**, wordmark → `/experiments`, plus the experiment header + `Corpus · Configuration` tabs.
3. **Chrome-less** — `/experiments` (home gallery) and `/experiments/new` render with no nav at all.

Symptoms: two near-identical nav components with divergent item lists; two wordmarks pointing at two different "homes"; `ExperimentCorpusPage` is a blank stub (its `corpus-sheet-slot` was never wired), so the only working contact sheet is still the legacy global `/samples`. Two organizing models (a global cross-experiment corpus filtered by beamtime vs. experiment == directory) run at once.

This is a half-finished cutover. The 2026-06-15 ingestion spec already calls for the collapse (retire the flat `/samples` root + beamtime chip; the corpus becomes a per-experiment tab; one `Experiments | Series` nav). This design completes that and resolves the seams the earlier spec left open.

## 2. Decisions (this session)

- **Experiment-first, fully unified.** One shell. The experiment is home base; the corpus, Configuration, grouping-review, and Focus are all experiment-scoped. The global `/samples` flat root + beamtime chip are retired.
- **Series is the one cross-experiment surface.** It genuinely spans experiments, so it stays a top-level destination (peer to Experiments) rather than a per-experiment view.
- **Two-tier chrome.** A lean top bar for *high-level* navigation; a light, contextual bottom **dock** for *in-context* navigation + actions.
- **No in-context container switcher.** You switch experiments via the Experiments landing gallery, and series via the Series folio. The dock never changes which container you are in.
- **One unified keyboard set**, designed holistically (existing bindings not held sacred). See §4.
- **Create-on-approval, permanent.** A new experiment is an ephemeral draft until the Configuration approval gate; Approve creates the DB row and starts the full parse. Cancel (pre-approval only) discards it; after approval there is **no in-app cancel or delete** (§6).
- **Surface layouts are out of scope.** The contact sheet, Focus workspace, Loupe, Series folio/builder/scoping, and Configuration body keep their refined layouts. This design touches the chrome around them, the navigation, and the keyboard wiring.

## 3. Architecture — the two-tier shell

### 3.1 Top tier (high-level nav, ~stable)

A single lean bar: **wordmark** (`Himalaya · SAXS`) · **section tabs** (`Experiments` / `Series`) · **⚙ Configuration**. No experiment switcher (removed — switching goes through the gallery). The current experiment's identity lives in the *page header* (on the surface), not the bar.

### 3.2 Bottom tier — the dock (in-context, contextual)

A persistent, **light** bar (plate + hairline + soft upward shadow; not the old dark CullBar). It is the evolution of the floating `CullBar` container, generalized: it re-populates with the current surface's in-context navigation and actions. It owns **only** connective/contextual verbs — never a refined surface's own command controls (those stay on the surface).

**Grammar (left → right), consistent on every surface:**

```
‹ up-link   |   cursor (steppers)   |   cull (when frames exist)   ……   destinations
```

- **up-link** (left anchor): one level up. Present on *every* surface so the controls never shift position. Corpus → `‹ Experiments`; Focus/Loupe → `‹ Corpus`; Series → `‹ All series`.
- **cursor**: the sample/frame steppers (mirrors the keyboard cursor, §4).
- **cull**: the verdict verbs — only where exposures exist (Corpus, Loupe). Absent on Focus and Series.
- **destinations** (right, primary rightmost): where you can go *down*/across — `Loupe`, `Focus`.

### 3.3 Per-surface dock

| Surface | up-link | cursor | cull | destinations |
|---|---|---|---|---|
| **Corpus** | `‹ Experiments` | Sample ↑↓ · Frame ‹› | Drop · Keep · Restore | Loupe · **Focus** |
| **Focus** | `‹ Corpus` | Sample ↑↓ | — | Loupe |
| **Loupe** | `‹ Corpus` | Sample ↑↓ · Frame ‹› | Drop · Keep · Set representative | **Focus** |
| **Series** | `‹ All series` | Sample ↑↓ | — | **Focus** |

Notes: Corpus keeps both destination buttons (Loupe and Focus) — they are the dock homes for the `L` and `Enter` keys, so a verbose "Open &lt;sample&gt; in Focus" label is unnecessary (clicking a row still opens it too). Series names the current sample (swatch + label + source experiment) so "which sample Focus opens" is unambiguous, and highlights that sample's trace on the overlay.

## 4. Keyboard model (one unified set)

The whole app shares a single keymap. It was designed holistically; existing bindings were **not** held sacred. Four rules:

1. **Bare arrows navigate** — `↑/↓` = which **sample**; `←/→` = which **detail within it** (a *frame* on Corpus/Loupe, a *candidate* on Focus). First-class — active without tabbing into a grid.
2. **Bare letters = verbs** on the current sample/frame (mnemonic).
3. **`Enter` = drill into Focus** (the primary "go"); **`Esc` = back up one level**, a *ladder* that dismisses the innermost transient state first.
4. **Modifiers are consistent**: `Alt` = move/reorder · `⌘/Ctrl` = meta (undo, find, select-all, confirm) · `Shift` = extend/range · `Space` = toggle-select. All bare keys suppress while typing in a text field; **`?`** shows the keymap overlay.

### 4.1 The map

**Navigate** — every surface where the unit exists (Corpus, Focus, Loupe, Series):

| Key | Action |
|---|---|
| `↑` / `↓` | Prev / next **sample** |
| `←` / `→` | Prev / next **detail** — frame (Corpus, Loupe) · candidate (Focus) |
| `Enter` | Open the current sample in **Focus** |
| `L` | Open the current sample in **Loupe** |
| `Esc` | Back / up one level (ladder) |
| `/` , `⌘K` | Find / jump (NavModal) |
| `?` | Keyboard-map overlay |

**Select & cull** — Corpus, Loupe (wherever a frame cursor exists):

| Key | Action |
|---|---|
| `Space` | Toggle the current frame in the selection |
| `Shift+←/→` | Extend the frame selection · `⌘A` select all |
| `X` / `K` | Drop / Keep (the selection, else the current frame) |
| `R` | Set the current frame as representative |
| `Backspace` | Restore — clear the verdict (selection / current) |
| `Esc` | Clear selection (first rung of Back) |

**Focus** — indexing:

| Key | Action |
|---|---|
| `←` / `→` | Prev / next candidate |
| `P` | Toggle add-peak mode |
| `Enter` | Apply the focused candidate / toggle the focused peak |
| `⌘Z` / `⌘⇧Z` | Undo / redo |

**Series edit** — builder, scoping:

| Key | Action |
|---|---|
| `Alt+↑` / `Alt+↓` | Reorder the sample up / down |
| `A` | Add a sample |
| `⌘Enter` | Confirm / build |
| `⌘Z` / `⌘⇧Z` | Undo / redo |

### 4.2 What changes vs. today (from the inventory)

- **Drop `[` / `]`** (sample-step) — `↑/↓` replaces it everywhere; no alias.
- **Candidate nav moves** `↑/↓` → `←/→` on Focus, freeing `↑/↓` for sample-step.
- **`Enter` now opens Focus** on Corpus (was: open loupe / sort header / cell-interaction). `L` (previously unbound) opens Loupe; the old `Shift+Enter` loupe-at-frame role folds into `L` + the frame cursor. Header-sort and cell-interaction `Enter` stay contextual (fire only on header / multi-widget cells).
- **Arrow nav no longer requires tab-in** — see §4.3.
- **New**: the `Space` / `Shift+arrows` selection model, `⌘A` select-all, `Backspace` restore, `?` help overlay.
- **Unchanged**: `X` / `K` (Drop/Keep, preserving the empty-selection no-op), `R` (set representative), `⌘Z` / `⌘⇧Z`, `/` + `⌘K`, `Alt+arrow` reorder.

### 4.3 Implementation notes

- **One registry.** Express the whole set through `src/print/shell/shortcuts.ts` + per-surface action wiring — not scattered listeners.
- **Tab-in / roving grid is the central change.** The contact sheet's roving-tabindex grid (`src/lib/grid/useRovingGrid.ts` + `rovingGrid.ts`) makes the table a single tab stop where arrows fire only after focus is inside (it stops propagation when focused). The new model promotes bare arrows to a **window-level** sample/frame stepper that drives the roving *active cell*; the grid's internal roving still functions once focused (a11y / screen-reader), but the common path no longer needs tab-in. Prevent double-handling.
- **Pointer parity.** Keyboard steppers must update the same active-cell state the pointer path uses (`requestActivate`), or mouse and keyboard diverge.
- **`Alt` guard.** Because bare `↑/↓` is now sample-step, the handler must ignore events when `Alt` is held so `Alt+arrow` reorder is never shadowed.
- **Input suppression.** Bare-letter / `Space` bindings stay suppressed inside text inputs (NavModal search, Series value fields, CustomIndexModal), matching the existing global gate.
- **Esc ladders preserved** per surface (Corpus: exit-interaction → clear-selection → up; Focus: clear-preview → disarm +Peak → up). "Clear selection" is one rung.

## 5. Visual treatment

- **Light dock** throughout (quiet furniture; the dark, attention-grabbing treatment was the *transient* CullBar's, which it earned by appearing only on selection).
- **Key chips**: a `currentColor`-based keycap, centered (flex + `line-height:1`). A **frosted** fill (translucent `currentColor`) is used **only on filled colored buttons** (the `Focus` primary, white frosted ↵). Neutral and colored-outline buttons get a plain chip whose key matches the text colour.
- **Focus button** = filled accent + frosted ↵ (the unambiguous primary). `Loupe` is a neutral pill; `Drop`/`Keep` are colored outlines; `Restore`/`Set representative` are neutral.
- **Uniform text size** across dock controls; hierarchy is carried by weight + colour, not size (e.g. the `Sample`/`Frame` labels recede via lighter colour at the same size).
- **Vocabulary**: the unit is **Sample** on every surface, including Series (members are samples).
- Tokens are "The Print" `@theme` set; appearance must live in `src/print/ui` primitives (the closed-look / open-placement contract; `lint:design` enforced).

## 6. Screen flow & states

### 6.1 First ingest — the funnel (once per experiment)

`/experiments` (gallery) → **New experiment** (directory picker) → **Configuration** → **Scan & review groups** → **Corpus**. The ingest is **two-phase**:

- **Picker — commit a path only.** The directory picker just commits a path (autocomplete + two cheap checks: the directory exists AND isn't already an experiment — backend 409 + inline guard). No indexing here.
- **Phase 1 — index (on Configuration; ephemeral, read-only).** Arriving at Configuration runs a cheap server-side index — file count, paths, patterns matched across image/metadata/integration, plus one sidecar for geometry. It persists nothing (a *manifest preview* — a richer form of the **planned** `/api/fs/validate` probe; **the `/api/fs/*` backend endpoints + the phase-1 index are not built yet — only the `validatePath`/`suggestPaths` frontend fetchers exist as E1 stubs**, see §8). Because it is **seconds-not-instant over SMB**, Configuration shows an **indexing progress** strip + a live "latest files scanned" stream while it runs, and **Approve sits disabled ("Indexing…") until it completes**. Editing a pattern restarts the index.
- **Approve** is the **creation boundary**: it persists the experiment row (config included) and starts **phase 2 — full parse** (PRP, geometry, grouping).
- **Phase 2 runs on the combined Scan & review groups surface** (§6.2): loads unfold live and each settled load is immediately reviewable; **Confirm groups** commits when the scan finishes → **Corpus** (or Corpus directly if there were no flags).

**Cancel** is available any time before Approve and discards the in-memory draft (aborting an in-flight phase-1 index); nothing is written. **After Approve the experiment is permanent — no in-app cancel or delete.** Existing experiments skip the funnel (gallery → their Corpus).

### 6.2 Scan + grouping (combined) and Corpus states

**First ingest** uses one continuous **Scan & review groups** surface: phase-2 loads unfold live; a settled load is immediately reviewable (its flags → Keep separate / Merge / Split via the existing structural-edit events) while later loads still parse. Safe because ingest is **additive** (a settled load is never re-touched) and edits go through the **event log** exactly as post-ingest. The final **Confirm groups** stays disabled until the scan completes (a later load can still raise a flag), then commits → Corpus. The **same component is the standalone `/grouping` route** reached from the Corpus banner after a rescan.

Once ingested, the **Corpus route is state-driven** — one route, body switches on state:

| State | Body |
|---|---|
| scanning (first ingest) | the combined **Scan & review groups** surface (full screen; loads unfold, flags reviewable, Confirm at end) |
| scan failed (first ingest) | **Scan-failed** surface — Open Configuration (primary); unmatched files + nearest existing file per type; adaptive pattern test; two-stage "Ingest N that parsed" |
| has grouping flags | the sheet + a sticky "N samples need grouping review →" banner → the standalone grouping route |
| clean | the sheet |
| rescanning | progress **inline** on Corpus (§6.3) — never the full-screen surface |

### 6.3 Rescan

Routine rescan (manual "Rescan" on Corpus, or auto on the tiered backoff) runs **in place on Corpus** — never back to the funnel or the Configuration gate:

1. Phase-1 cheap check against the stored `scan_signature` for new/changed files.
2. Nothing new → silent "up to date".
3. New files → **additive** phase-2 parse of the delta only (existing samples are never re-parsed or re-grouped, so curation is preserved); the sheet grows in place, and new flags raise the review banner.
4. Errors → an inline notice/toast, not the full-screen Scan-failed page.

The **only gated rescan** is a **Configuration / pattern change** (via ⚙) — editing patterns re-runs phase 2, because pattern changes can re-include/exclude files.

### 6.4 Routing map

- `/` → `/experiments` (gallery, the nav-default).
- `/experiments/new` — directory picker (draft; no row until Approve).
- `/experiments/:id/corpus` — home base (state-driven, §6.2).
- `/experiments/:id/config` — Configuration (the first-run gate *and* the later ⚙ edit reuse this route).
- `/experiments/:id/grouping` — grouping-review surface (from the banner; returns to corpus when clean).
- Focus + Loupe are experiment-scoped routes under the experiment.
- `/series` (folio) + scoping/builder — cross-experiment; reached from the top `Series` section and via "send to series" from a corpus selection. A member opens via `Enter` → Focus, with the return thread reading `‹ Series`.
- **Switch** experiment via the gallery, series via the folio (no in-context switcher).
- **Retired**: the flat `/samples` root, the `?beamtime` chip, the `Samples` stage tab, the second wordmark/home, the `ExperimentTopNav` vs `CorpusTopbar` duplication (collapse to one `TopNav`).

### 6.5 New surfaces (per-page behaviors)

The surfaces this design introduces (drafted as mockups, §9). Refined surfaces (contact sheet, Focus, Loupe, Series, the Configuration *body*) are reused unchanged.

- **Experiments gallery** (`/experiments`) — cards grouped under **year** headings (newest first) + a slim left **timeline rail**; per-card state chip (up-to-date / N-to-review / scanning); `+ New experiment`. The gallery *is* the experiment switcher.
- **New experiment** (`/experiments/new`) — directory picker: path autocomplete + the two checks (exists, not-already-an-experiment); `Review →` commits the path; no indexing here.
- **Configuration** (`/experiments/:id/config`) — first-run gate *and* later ⚙ edit. Runs phase-1 indexing in place (progress strip + latest-files stream). Sources patterns + auto-derived geometry read as **one aligned key/value table** with **quiet gray** provenance tags (prp/setup/computed); **undo/redo** on the editable fields (reuses ⌘Z/⌘⇧Z). Approve = creation boundary (disabled "Indexing…" until phase-1 done).
- **Scan & review groups** (combined; first ingest + the `/grouping` route) — see §6.2.
- **Scan failed** — Open Configuration (primary); a **scrollable** list of all unmatched files grouped by miss type, each paired with the **nearest existing file**; an **adaptive pattern test** (one field per affected type — image/metadata/integration — clearing independently) + "Apply all in Configuration"; "Ingest N that parsed" is a real button with a two-stage in-place confirm.
- **Corpus assembly** — the daily-loop home: top bar + experiment header (name, dir, rescan status + Rescan, stat masthead) + the review banner + the reused contact sheet (2D cursor) + the dock.

## 7. Scope

**In scope:** the unified two-tier shell; one `TopNav`; the contextual dock component + its per-surface population; the one unified keyboard set (§4) wired through the central registry + the roving-grid handler; routing collapse (one shell, retire `/samples` + redirects); wiring the per-experiment corpus sheet (`ExperimentCorpusPage` → the existing `SheetTable`, scoped).

**Out of scope (untouched):** the contact sheet table, the Focus plate + assignment sidebar, the Loupe frame/side panel, the Series folio/builder/scoping internals, the Configuration body. Their *layouts*, columns, and in-surface command controls (scale toggle, auto-fit, +Peak, export, reorder, undo/redo, confirm/cancel) stay exactly as refined.

**Scope nuance — keyboard:** delivering the unified keyboard set necessarily adjusts some in-surface *key handling* (re-scoping Focus candidate navigation; the Corpus `Enter`/`L` policy; removing the grid's tab-in requirement). This is distinct from surface *visual/layout*, which is untouched — the changes are expressed through `shortcuts.ts` + the roving-grid handler, not by re-laying-out any surface.

## 8. Open / deferred

- Whether the dock is hidden vs. present-but-quiet on a surface with nothing contextual to show (resolved for now: always present, because the up-link + cursor always have something to say).
- Final placement of the grouping-review banner relative to the corpus dock.
- **Backend dependency (not yet built):** the directory-picker + phase-1 index endpoints — `/api/fs/suggest`, `/api/fs/validate`, and the ephemeral phase-1 directory index — do not exist server-side; only the frontend `suggestPaths`/`validatePath` fetchers (E1 stubs) do. The funnel's picker + Configuration-indexing depend on building these.
- Migration sequencing (this is a follow-on plan, not part of this spec).

## 9. Reference

Interactive mockups (rendered at the exact Print tokens) developed during this session, served locally at `/tmp/shell-mockups/`:
- `bars.html` — the converged dock set per surface (Corpus / Focus / Loupe / Series).
- `glyph.html` — the key-chip treatment (↵, scoped frost).
- `flow.html` — the screen-flow & states diagram (funnel, creation boundary, daily loop, rescan, series lane).
- `pages1.html` — the funnel surfaces (gallery, picker, Configuration, Scan & review groups, Scan failed).
- `pages2.html` — the Corpus assembly.

These are throwaway exploration artifacts, not committed source.
