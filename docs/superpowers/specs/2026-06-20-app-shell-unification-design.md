# App Shell Unification — Design

> Status: **DRAFT (review incorporated)** · 2026-06-20 · branch `ingestion-redesign`
> Scope: the app **shell, the navigation between surfaces, and one unified keyboard set** — NOT the refined surfaces' layouts.
> Rev 2 folds in a 5-dimension spec review (24 confirmed findings): ⚙ moved off the top bar (§3.1), dock/keymap cull verbs reconciled per-surface (§3.3/§4), the roving grid is **dropped** for the sheet in favour of nested 1-D cursors (§4.3), the keyboard set is noted as superseding the 2026-06-13 lock (§4), and §6 is corrected to the **live backend** — the progress rail (`broadcast_progress!`, `ingest_status`, the async create route) already exists; only the phase-1 manifest endpoints are net-new.

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

A single lean bar: **wordmark** (`Himalaya · SAXS`) · **section tabs** (`Experiments` / `Series`). No experiment switcher (removed — switching goes through the gallery).

**⚙ Configuration is NOT a top-bar item** (review A). Configuration is experiment-scoped (`/experiments/:id/config`, §6.4) and meaningless on the gallery / Series where no experiment is current. It lives **with the experiment's identity in the page header** (the same place the bar deliberately keeps experiment identity off the global bar). So the top bar is genuinely route-stable — it never gains or loses items — and the ⚙ edit affordance rides the experiment header, present exactly when (and only when) an `:id` is in context.

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
| **Loupe** | `‹ Corpus` | Sample ↑↓ · Frame ‹› | Drop · Keep · Set representative · Restore | **Focus** |
| **Series** | `‹ All series` | Sample ↑↓ | — | **Focus** |

**Cull verbs are per-surface** (review B), and the dock table is authoritative — the §4 keymap matches it exactly: **Corpus** has no representative concept (the sheet operates on samples/frames, no per-sample representative frame), so its verbs are `{Drop, Keep, Restore}` and `R` is **not** bound on Corpus. **Loupe** operates within one sample's frames, so it gets the full `{Drop, Keep, Set representative, Restore}` and binds `R`. `Restore` (`Backspace`) exists on both.

Notes: Corpus keeps both destination buttons (Loupe and Focus) — they are the dock homes for the `L` and `Enter` keys, so a verbose "Open &lt;sample&gt; in Focus" label is unnecessary (clicking a row still opens it too). Series names the current sample (swatch + label + source experiment) so "which sample Focus opens" is unambiguous, and highlights that sample's trace on the overlay.

## 4. Keyboard model (one unified set)

The whole app shares a single keymap. It was designed holistically; existing bindings were **not** held sacred. Four rules:

1. **Bare arrows navigate** — `↑/↓` = which **sample**; `←/→` = which **detail within it** (a *frame* on Corpus/Loupe, a *candidate* on Focus). First-class — active without tabbing into a grid.
2. **Bare letters = verbs** on the current sample/frame (mnemonic).
3. **`Enter` = drill into Focus** (the primary "go"); **`Esc` = back up one level**, a *ladder* that dismisses the innermost transient state first.
4. **Modifiers are consistent**: `Alt` = move/reorder · `⌘/Ctrl` = meta (undo, find, select-all, confirm) · `Shift` = extend/range · `Space` = toggle-select. All bare keys suppress while typing in a text field; **`?`** shows the keymap overlay.

> **This keymap supersedes the 2026-06-13 keyboard-shortcut-library lock** (review C). That lock standardized the *previous* greenfield gestures (`[`/`]` sample-step, `↑/↓` candidate-nav, `Alt+arrow` reorder); this design deliberately reassigns the navigation axis app-wide and is the new source of truth. The reorder, undo/redo, and find bindings from that lock are preserved unchanged (§4.2).
>
> **Focus exposure axis (review C):** Focus shows the sample's **representative exposure**; `←/→` on Focus is the *candidate* cursor, not an exposure stepper. Stepping between a sample's exposures stays where it already lives — the **rail-scoped** exposure control (the existing `FO-EXPSKIP` affordance) — and is intentionally *not* on bare arrows. The bare-arrow axis is one consistent thing everywhere (sample ↑↓ / detail ←→); exposure-within-sample is a rail action, not a global gesture.

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
| `X` / `K` | Drop / Keep (the selection, else the current frame) — Corpus, Loupe |
| `R` | Set the current frame as representative — **Loupe only** (no representative on Corpus) |
| `Backspace` | Restore — clear the verdict (selection / current) — Corpus, Loupe |
| `Esc` | Clear selection (first rung of Back) — Corpus, Loupe |

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

- **One registry, one semantic ID per key (review G9).** Express the whole set through `src/print/shell/shortcuts.ts`. `matchShortcut` is a **flat global registry** (no surface scope) and resolves first-wins — so the same physical key MUST map to a single registry ID, and the *page* interprets it. Bind `←/→` to one `prevDetail`/`nextDetail` pair (each surface reads "detail" as frame or candidate); bind `Enter` to one `drillIn` ID; bind `Space` to one `toggleSelect`. Do **not** register two IDs on `ArrowLeft`/`Enter` and disambiguate by surface — first-wins would silently break one surface.
- **The roving grid is removed for the sheet (review D — central simplification).** With the open-button column gone, the contact sheet's only interactive targets are the **sample row** (`↑/↓`) and its **exposures** (`←/→`) — nested 1-D cursors, not a 2-D grid. So `src/lib/grid/useRovingGrid.ts` + `rovingGrid.ts` are dropped, and `SheetTable.tsx` (`roving` prop) + `SampleTableRow.tsx` (`useRovingGridContext`) revert to a plain table. Blast radius is exactly those 3 files (sole consumers; inert context default → nothing else can break). The cursor becomes a page-level `{sampleIndex, frameIndex}` and there is a **single** window-level arrow/Enter handler — eliminating the dual-handler double-dispatch problem entirely (no grid handler to arbitrate against).
- **The one handler DECLINES** (returns false, the existing convention) when (a) `document.activeElement` is inside a text input / open modal (NavModal, CustomIndexModal, Series value fields — the existing structural `tagName`/`contenteditable` gate, which also covers the new bare `Space`/`Enter`), and (b) for `Enter`, the target is a native interactive element (`button`/`a`/`[role=button]`/sort header). State this as a testable rule, not "prevent double-handling."
- **Pointer parity.** Clicking a row or an exposure sets the same page-level `{sampleIndex, frameIndex}` the keyboard drives — one cursor source, so mouse and keyboard cannot diverge. The dock steppers call the same cursor setter.
- **Per-row Restore (review D).** A dropped row gets a real inline `<button>Restore</button>` (today Restore lives only in the floating multi-select bar; a single dropped row has no affordance). This gives `Backspace` a visible target and a pointer path.
- **`?` overlay combo (review G8).** Special-case `e.key === '?'` in `eventCombo` to emit a stable `?` token regardless of the Shift bit (layout-robust, matches the bare-`?` notation) rather than registering `Shift+/` (which collides with `/` find under normalization). Add a normalization unit test for `?` alongside the existing CapsLock-`X` case. The overlay component itself is net-new (today `KbdLegend` is per-group inline, not a global modal) — build it registry-driven so it lists the live keymap.
- **`Alt` guard.** Because bare `↑/↓` is now sample-step, the handler must ignore events when `Alt` is held so `Alt+arrow` reorder is never shadowed.
- **Esc ladders preserved** per surface (Corpus: clear-selection → up; Focus: clear-preview → disarm +Peak → up). "Clear selection" is one rung. On the gallery/home (top of the ladder) `Esc` is a no-op.

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

- **Picker — commit a path only.** The directory picker just commits a path (autocomplete via a new `GET /api/fs/suggest` + two cheap checks: the directory exists AND isn't already an experiment — the create route already 409s on a duplicate `data_dir`; the picker also guards inline). No indexing here.
- **Phase 1 — file manifest (on Configuration; ephemeral, read-only).** Arriving at Configuration runs a cheap server-side **file** index: `readdir` + glob-match against the image/metadata/integration patterns → `{total, matched-by-type, unmatched grouped by miss-type + the nearest existing file per miss}`. It **does not parse PRPs and writes nothing** — PRP/geometry parsing is deliberately phase 2 (this is the cheap-listing/full-parse split you asked for). Because it's pure listing it is fast even over SMB. New endpoint: **`GET /api/fs/manifest?path=…&{image,metadata,integration}_pattern=…`**, keyed by `(path, patterns)` — editing a pattern just re-fetches (request-scoped, naturally abortable). Configuration shows the manifest result + an indexing spinner while the request is in flight, and **Approve sits disabled until it returns.**
  - *Backend status:* the `/api/fs/*` endpoints are **net-new** (today only the `validatePath`/`suggestPaths` frontend fetchers exist as E1 stubs, hitting nothing). They factor the glob-match half out of the existing `scan_directory` (which already separates listing from parsing). This is the **only genuinely-new backend surface** the funnel needs.
- **Approve = the creation boundary, and it is essentially the route that already exists.** `POST /api/experiments` today already: validates `isdir`, rejects a duplicate dir (409), inserts the row with `ingest_status='scanning'`, returns **202**, and `Threads.@spawn`s `scan_and_group!` — i.e. it is *already* create-then-async-scan. So there is **no draft DB row and no separate approve route**: the pre-approval draft (path + edited patterns + the phase-1 manifest) lives entirely **client-side**, and Approve is this POST. The one change: the create handler should **accept the edited patterns in its body** (today patterns are a follow-up PATCH) so Configuration's edits flow straight into the phase-2 scan.
- **Phase 2 = `scan_and_group!`, on the combined Scan & review groups surface** (§6.2). It parses PRP/geometry, computes grouping, commits the structure in one fast transaction, then analyzes **per-exposure outside the transaction** (peaks/phase commit row-by-row). The structure is therefore queryable the instant the fast write lands; analysis streams in. **Confirm groups** commits when the scan finishes → **Corpus** (or Corpus directly if there were no flags).

**Cancel** is available any time before Approve and discards the **client-side** draft (and aborts the in-flight phase-1 `fetch`); nothing is written server-side because no row exists yet. **After Approve the experiment is permanent — no in-app cancel or delete** (a `DELETE /api/experiments/:id` route exists but the funnel does not surface it; it stays for CLI/admin). Existing experiments skip the funnel (gallery → their Corpus).

### 6.2 Scan + grouping (combined) and Corpus states

**First ingest** uses one continuous **Scan & review groups** surface. The "live" behaviour rides the mechanism that already exists, not a transaction rewrite: the grouping **structure** (loads/samples/exposures) commits in one fast transaction and is queryable immediately, then per-exposure **analysis streams in** (peaks/phase commit row-by-row outside the transaction). A settled sample is reviewable (its flags → Keep separate / Merge / Split via the existing structural-edit events) the moment the structure lands, while analysis fills the rest in. Safe because ingest is **additive** (a settled load is never re-touched) and edits go through the **event log** exactly as post-ingest.

The progress signal is **already plumbed**: `broadcast_progress!` emits transient `ingest_{started,progress,complete,failed}` SSE frames (no `user_actions` row), the frontend has the matching `ingest_*` cache arms + `ingestInFlight`, and `ingest_status` is a column the create route already drives. The **only backend change for "live"** is to give `scan_and_group!` an `on_progress` callback fired (throttled) inside its per-exposure analysis loop, which the create handler wires to `broadcast_progress!(kind="ingest_progress", processed, total)`. We deliberately **keep the single structural transaction** (no per-load commit) — it buys nothing here (the structural write is already fast) and would trade away the clean "structure appears whole, then peaks fill in" story; `ingest_complete` stays the authoritative terminal frame (intermediate ticks may drop under the 64-deep channel cap, which the helper already documents).

The final **Confirm groups** stays disabled until `ingest_status === 'complete'` (a later exposure can still raise a flag), then commits → Corpus. The **same component is the standalone `/grouping` route** reached from the Corpus banner after a rescan.

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

1. Phase-1 cheap check for new files — `cheap_change_check` (exists today) compares the on-disk file count against the persisted exposure count. *Note:* the `scan_signature` column exists but is currently **unused** by the cheap check (set to `NULL` on a pattern change as a rescan signal with no reader); strengthening change-detection to hash file metadata into `scan_signature` is a small optional follow-up, not required by this design.
2. Nothing new → silent "up to date".
3. New files → **additive** phase-2 parse of the delta only (existing samples are never re-parsed or re-grouped, so curation is preserved); the sheet grows in place, and new flags raise the review banner.
4. Errors → an inline notice/toast, not the full-screen Scan-failed page.

The **only gated rescan** is a **Configuration / pattern change** (via ⚙) — editing patterns re-runs phase 2, because pattern changes can re-include/exclude files.

### 6.4 Routing map

- `/` → `/experiments` (gallery, the nav-default).
- `/experiments/new` — directory picker (draft; no row until Approve).
- `/experiments/:id/corpus` — home base (state-driven, §6.2).
- `/experiments/:id/config` — Configuration (the first-run gate *and* the later edit, the latter reached via the **⚙ in the experiment header**, §3.1 — same route either way).
- `/experiments/:id/grouping` — grouping-review surface (from the banner; returns to corpus when clean).
- Focus + Loupe are **experiment-scoped** routes nested under the experiment: `/experiments/:id/sample/:sampleId` (Focus) and `/experiments/:id/sample/:sampleId/loupe` (Loupe). The current legacy `/sample/:id` redirects to its experiment-scoped form during the routing collapse.
- `/series` (folio) + scoping/builder — cross-experiment; reached from the top `Series` section and via "send to series" from a corpus selection. A Series member opens Focus at its own experiment-scoped route carrying a `from=series` return marker, so the up-link reads `‹ Series` (rather than `‹ Corpus`) for that visit — resolving the dock-vs-route back-target ambiguity.
- **Switch** experiment via the gallery, series via the folio (no in-context switcher).
- **Retired**: the flat `/samples` root, the `?beamtime` chip, the `Samples` stage tab, the second wordmark/home, the `ExperimentTopNav` vs `CorpusTopbar` duplication (collapse to one `TopNav`).

### 6.5 New surfaces (per-page behaviors)

The surfaces this design introduces (drafted as mockups, §9). Refined surfaces (Focus, Loupe, Series, the Configuration *body*) are reused unchanged; the contact sheet is reused with the interaction changes in §7 nuance (open-button column removed, roving grid dropped, per-row Restore added) — its table layout is otherwise unchanged.

- **Experiments gallery** (`/experiments`) — cards grouped under **year** headings (newest first) + a slim left **timeline rail**; per-card state chip (up-to-date / N-to-review / scanning); `+ New experiment`. The gallery *is* the experiment switcher.
- **New experiment** (`/experiments/new`) — directory picker: path autocomplete + the two checks (exists, not-already-an-experiment); `Review →` commits the path; no indexing here.
- **Configuration** (`/experiments/:id/config`) — first-run gate *and* later ⚙ edit (reached from the experiment header, §3.1), same route in two modes:
  - **First-run (pre-approval):** runs the phase-1 **file** manifest in place (spinner + matched-by-type counts + unmatched list). The editable **patterns** are the focus; **undo/redo** on those fields (reuses ⌘Z/⌘⇧Z). Geometry is **not yet shown** — it derives from PRP in phase 2, which hasn't run (phase 1 lists files, it doesn't parse). Approve = creation boundary (disabled until the manifest returns).
  - **Later ⚙ edit (post-ingest):** patterns **plus** the auto-derived geometry/sources as **one aligned key/value table** with **quiet gray** provenance tags (prp/setup/computed) and per-field override; saving re-runs the gated phase-2 rescan (§6.3).
- **Scan & review groups** (combined; first ingest + the `/grouping` route) — see §6.2.
- **Scan failed** — Open Configuration (primary); a **scrollable** list of all unmatched files grouped by miss type, each paired with the **nearest existing file**; an **adaptive pattern test** (one field per affected type — image/metadata/integration — clearing independently) + "Apply all in Configuration"; "Ingest N that parsed" is a real button with a two-stage in-place confirm.
- **Corpus assembly** — the daily-loop home: top bar + experiment header (name, dir, rescan status + Rescan, stat masthead) + the review banner + the reused contact sheet (nested sample/frame cursor, §4.3) + the dock.

## 7. Scope

**In scope:** the unified two-tier shell; one `TopNav`; the contextual dock component + its per-surface population; the one unified keyboard set (§4) wired through the central registry + a single page-level cursor handler; **removing the roving grid** (`useRovingGrid`/`rovingGrid`, the `SheetTable` `roving` prop, the `SampleTableRow` context) in favour of nested 1-D cursors (§4.3); routing collapse (one shell, retire `/samples` + redirects); wiring the per-experiment corpus sheet (`ExperimentCorpusPage` → the existing `SheetTable`, scoped); the net-new backend (`GET /api/fs/manifest` + `/api/fs/suggest`; `scan_and_group!` `on_progress` callback; patterns in the create body — §6.1/§6.2).

**Out of scope (untouched layouts):** the Focus plate + assignment sidebar, the Loupe frame/side panel, the Series folio/builder/scoping internals, the Configuration body. Their *layouts* and in-surface command controls (scale toggle, auto-fit, +Peak, export, reorder, undo/redo, confirm/cancel) stay exactly as refined.

**Scope nuance — what this design *does* reach into** (so the plan doesn't treat these as untouched):
- **The contact sheet** keeps its table layout/columns, but: the open-button column is removed (rows/`Enter` open Focus), the roving grid is dropped, and a per-row `Restore` button is added to dropped rows (§4.3). These are interaction changes, not a re-layout.
- **FocusPage key handling** is edited — candidate nav re-scopes `↑/↓ → ←/→`, and the bare-arrow sample axis is added. This is distinct from Focus's *visual/layout*, which is untouched; the changes live in `shortcuts.ts` + the page cursor handler.

## 8. Open / deferred

- Whether the dock is hidden vs. present-but-quiet on a surface with nothing contextual to show (resolved for now: always present, because the up-link + cursor always have something to say).
- Final placement of the grouping-review banner relative to the corpus dock.
- **Backend — what already exists vs. net-new** (verified against live source, rev 2):
  - *Already built (reuse):* `POST /api/experiments` (validate-isdir → 409-on-dup → create `ingest_status='scanning'` → 202 → async `scan_and_group!`); `broadcast_progress!` + the four transient `ingest_*` SSE frames + the frontend `ingest_*` cache arms + `ingestInFlight`; the `ingest_status` state column; `scan_and_group!` (idempotent dedup-key ingest, fast structural tx + per-exposure analysis outside it); `cheap_change_check`; the tiered-backoff rescan scheduler; `DELETE /api/experiments/:id` (not surfaced by the funnel).
  - *Net-new (build):* `GET /api/fs/manifest` (phase-1 file listing) + `GET /api/fs/suggest` (picker autocomplete) — the `validatePath`/`suggestPaths` frontend fetchers are E1 stubs hitting nothing today; an `on_progress` callback on `scan_and_group!` (for `ingest_progress` ticks); accepting patterns in the create body; the `?` keymap-overlay component; the per-row Restore button.
- Optional: hashing file metadata into the unused `scan_signature` column for stronger rescan change-detection (§6.3).
- Migration sequencing (this is a follow-on plan, not part of this spec).

## 9. Reference

Interactive mockups (rendered at the exact Print tokens) developed during this session, served locally at `/tmp/shell-mockups/`:
- `bars.html` — the converged dock set per surface (Corpus / Focus / Loupe / Series).
- `glyph.html` — the key-chip treatment (↵, scoped frost).
- `flow.html` — the screen-flow & states diagram (funnel, creation boundary, daily loop, rescan, series lane).
- `pages1.html` — the funnel surfaces (gallery, picker, Configuration, Scan & review groups, Scan failed).
- `pages2.html` — the Corpus assembly.

These are throwaway exploration artifacts, not committed source.
