# App-Shell Unification — Canonical Spec

**Date:** 2026-06-22 · branch `ingestion-redesign` · Status: **THE single source of truth**

This is the **only** app-shell / ingestion spec. It consolidates and replaces three earlier docs,
which were **mined and deleted** on 2026-06-22 (their full text is in git history):

- `2026-06-15-ingestion-redesign-design.md` (data model + ingest pipeline rationale → §11–§12)
- `2026-06-20-app-shell-unification-design.md` (shell/nav/keyboard architecture → §3, §7, §8)
- `2026-06-21-ingestion-appshell-consolidated.md` (milestone/delta tracker → §10)

This document is the **target**: it transcribes the canonical mockups at copy-and-pixel level, folds
in the data-model/backend internals from the deleted specs, and folds in every decision from the
2026-06-22 live walkthrough + a 10-source conversation-mining sweep. It states what is built vs. what
remains. (The phase plans under `docs/superpowers/plans/2026-06-18-*` still reference the deleted
specs by name; those are historical implementation records, left as-is.)

> **Reading rule — this spec IS the design target.** Every surface description is what to *build*.
> Where the spec cites shipped code (`file:line`, "shipped …", "fidelity gap", "§10 Remaining"), that
> documents the **current state and the gap to close** — it is **never** a license to downgrade the
> target to match what happens to be shipped. When shipped and the mockup disagree, **the mockup (and
> the ratified decisions in this doc) win**; the shipped state is the work queue, not the goal. The
> only deviations *from* the mockup that are canonical are the ones explicitly ratified here as a
> decision (e.g. the kept Tags column, §5.6).
>
> **Code citations (`file:line`) are indicative anchors, not contracts.** Line numbers drift as the
> branch moves; when a number and a named symbol disagree, trust the **symbol**. Don't treat a stale
> line number as a spec defect. Frontend (`.ts`/`.tsx`) cites in the backend sections especially are
> indicative — re-verify against live code before relying on an exact line.

**Canonical visual reference (the design target):** `docs/redesign-mockups/app-shell/`. The
optimized renders are authoritative. (The stale Jun-18 top-level `docs/redesign-mockups/*.html` set was
**deleted** 2026-06-22 — only the `app-shell/` set remains; history is in git.)

| Surface | Render | Source |
|---|---|---|
| Experiments gallery | `p1-gallery.png` | `pages1.html` §1 |
| New-experiment picker | `p1-new.png` | `pages1.html` §2 |
| Configuration (first-run) | `p1-config.png` | `pages1.html` §3 |
| Scan + grouping (combined) | `p1-grouping.png` / `p1-scanning.png` | `pages1.html` §4 |
| Scan failed | `p1-failed.png` | `pages1.html` §5 |
| Corpus assembly (daily home) | `p2-corpus.png` | `pages2.html` |
| Dock per surface | `b1`–`b4.png` | `bars.html` |
| Screen flow & states | `flow.png` | `flow.html` |
| Key-chip treatment | `glyph.png` | `glyph.html` |

---

## 1. What this is, in one paragraph

Himalaya is **experiment-first and fully unified**: one shell, one top bar, one contextual bottom
dock, one keyboard set. The experiment is home base. You **point Himalaya at a directory and it
figures the rest out** (Gradescope-style: dump a dir, it discovers the layout, scans the exposures,
groups them into samples, and hands you a review gate). The corpus contact sheet is the daily-loop
home; Configuration, grouping-review, Focus, and Loupe are all experiment-scoped; Series is the one
cross-experiment surface. There is no global flat corpus and no in-context container switcher — you
switch experiments through the gallery and series through the folio.

---

## 2. Design language

Every surface renders at "The Print" `@theme` tokens. The mockups encode the exact values; these
are the load-bearing ones (frontend tokens already match — listed so the spec is self-contained).

**Color (OKLCH, never `#000`/`#fff`, every neutral tinted warm ~hue 80):**
```
--paper        0.978 0.006 85     page background
--paper-sunk   0.951 0.008 84     sunk wells, row fills, key fields
--plate        0.992 0.004 90     cards, sheet, dock
--ink          0.265 0.013 68     primary text
--ink-soft     0.467 0.012 70     secondary text, labels
--ink-faint    0.640 0.010 74     tertiary, hints, "/ total" denominators
--hair         0.882 0.008 80     hairline dividers
--hair-strong  0.806 0.010 78     card borders, controls
--accent       0.555 0.150 38     terracotta — the ONE accent (kickers, primary, Drop, into-Focus)
--success      0.520 0.120 162    sage-green — Keep, "up to date", clean (code names it "sage")
--warning      0.620 0.130 70     amber — flags, "N to review", scan-failed, ambiguity
--frame-edge   0.150 0.010 55     detector-thumb background
```
Phase chips carry their own per-phase hues, authored as a TS `PHASE_PALETTE` record of OKLCH literals
in `src/phases.ts`:32-40 (there is **no** `--phase-*` CSS var). `phases.ts` is the one sanctioned
color-authoring file outside `ui/`.

**Type:** serif `Newsreader` 500 for display/names (experiment titles, sample names, year heads);
sans `Plus Jakarta Sans` for body/UI; mono `ui-monospace` for paths, ids, counts, file patterns,
and the uppercase-tracked kickers/labels. Hierarchy is scale + weight, never decoration.

**Surfaces:** cards are `--plate` + `--hair-strong` border + `--shadow-plate-lift` (a 1px inset
white highlight + a soft drop; the token is defined in `styles.css`). The dock and action bars get a
soft *upward* shadow so they read as furniture lifting off the page from below — note there is **no
`--shadow-dock` token today**; `Dock.tsx`:36 uses an inline
`shadow-[0_-2px_8px_0_rgba(0,0,0,0.06)]` (8px blur). Either define the token or treat that inline
value as canonical. **No side-stripe accents** (machine-enforced by `lint:design`, the `no-side-stripe`
rule in `check-design.mjs`). **No nested cards, no gradient text** (design conventions, NOT
machine-enforced — `check-design.mjs` has no rule for either).

**Copy law:** no em dashes in user-facing prose (use `·`, commas, periods). Enforced by
`test/print-copy-no-emdash.test.ts`. The placeholder literal `—` is sanctioned.

---

## 3. Information architecture

### 3.1 The unified shell — two tiers

**Top tier (lean, route-stable).** A single `TopNav`: wordmark `Himalaya · SAXS` (serif) · two
section tabs `Experiments` / `Series` (uppercase; each tab carries a terracotta dot unconditionally,
`TopNav.tsx`:54 — the **active** tab additionally gets a `--paper-sunk` fill + `text-ink`, `:49`) · a
flex spacer. **The `TopNav` carries NO ⚙ gear** (`TopNav.tsx`, the JSDoc at `:18-22`
says "no gear"); the bar **never gains or loses items** across routes.

**⚙ is NOT a global nav item.** It is the experiment-scoped Configuration affordance and is
meaningless on the gallery / Series. It rides the **experiment header** only
(`ExperimentShell.tsx`:140-147, `data-testid="experiment-config-gear"`), present exactly when an
`:id` is in context. The gallery and Series render no gear at all.

**Bottom tier — the dock (contextual).** A persistent **light** bar (plate + hairline + a soft upward
shadow `shadow-[0_-2px_8px_0_rgba(0,0,0,0.06)]`, the canonical inline value — there is no
`--shadow-dock` token, see §2). It is the generalized evolution of the old floating `CullBar`. It owns
**only connective verbs** — up-link, cursor steppers, cull, destinations — never a refined surface's
own command controls. Full grammar in §7.

There is exactly **one** `TopNav` and **one** `Dock`. The legacy `CorpusShell` / `CorpusTopbar` /
`ExperimentTopNav` triplet and the chrome-less home are all collapsed into this. (Done — §10.)

### 3.2 Routing map

```
/                          → /experiments               (gallery is the default)
/experiments               gallery / home               (the experiment switcher)
/experiments/new           directory picker (draft; NO db row until Approve)
/experiments/new/config    first-run Configuration (draft route — no :id yet)
/experiments/:id/corpus    home base — STATE-DRIVEN (§6.6)
/experiments/:id/config    Configuration (first-run gate AND later ⚙ edit; same route, two modes)
/experiments/:id/grouping  Scan + grouping review (combined surface; takeover, OUTSIDE the shell)
/sample/:sampleId          Focus            (flat leaf route; param key is `sampleId`)
/sample/:sampleId/loupe    Loupe            (flat leaf route)
/series                    Series folio + scoping + builder (cross-experiment)
```

**Retired:** the flat `/samples` global corpus, the `?beamtime` filter chip, the `Samples` stage
tab, the second wordmark/home, `SamplesPage`. `/samples` → `/experiments` redirect kept only as a
courtesy.

**Why Focus/Loupe are flat:** a sample id is globally unique, and Focus/Loupe are reached from the
corpus *and* from cross-experiment Series. A flat address is honest and avoids a fetch-then-redirect
loading flash. The chrome (header, `‹ Corpus` up-link) resolves the experiment from the loaded
sample's `experiment_id` (which the page fetches anyway). A Series→Focus visit carries a
`?from=series` query param (honored by **Focus only**, `FocusPage.tsx`:511,863) so the Focus up-link
reads `‹ Series`; the marker, not the path shape, resolves the back-target. Loupe's up-link is fixed
`‹ Corpus` and its Esc uses `goBack()` (`LoupePage.tsx`:366), so a Series→Loupe return goes through
history, not the marker.

**Why `/grouping` is a takeover (outside `ExperimentShell`):** it is full-screen and renders its own
single header. Nesting it inside the shell double-stacks headers. The corpus "scanning" state
**redirects** to `/grouping`; the funnel's Approve lands on `/grouping`; it returns to
`/corpus` when clean / confirmed.

### 3.3 Settled IA decisions — do not re-litigate

1. **One shell, experiment-first.** Global `/samples` + beamtime chip retired.
2. **No Corpus | Configuration tab bar.** Corpus *is* the experiment home; Configuration lives
   behind the ⚙ gear in the experiment header.
3. **Series stays a top-level peer** (it genuinely spans experiments).
4. **No in-context switcher.** Switch experiments via the gallery, series via the folio. The dock
   never changes which container you are in.
5. **Create-on-approval, permanent.** A new experiment is a client-side draft until Approve. Approve
   creates the row and starts the parse. Cancel (pre-approval only) discards it. **After approval
   there is no in-app cancel or delete** (the `DELETE` route exists for CLI/admin; the funnel does
   not surface it).
6. **Refined surface *layouts* are out of scope** (Focus plate + assignment sidebar, Loupe frame +
   side panel, Series folio/builder/scoping internals, the Configuration *body* widgets). This design
   touches their chrome, navigation, and keyboard wiring — not their internal layout.

> Note on decision 5: "no in-app cancel/delete after approval" reads together with the known
> crash-mid-scan edge (§6.6 / §10 G4) — an experiment stuck at `ingest_status='scanning'` after a crash
> has no in-app recovery today; that gap is tracked, not a contradiction of decision 5.

---

## 4. The two-phase ingest funnel

This is the spine of the redesign. Transcribed from `flow.html` (the optimized flowchart) and
`pages1.html`.

### 4.1 The flow (flowchart, verbatim intent)

```
Experiments ──+ New──▶ New experiment ──create·index files──▶ Configuration ──confirm·parse──▶ Scanning… ──flags──▶ Grouping review ──resolved──▶ Corpus
  (gallery)             (directory picker)   │                  (check files·        (full parse·                  (Gradescope gate)            (home base·
                                             │                   patterns·geometry)   unfold)                                                     sheet·banner·↻)
                                             │                                          │                                                            │
                              ⟵ draft · created here ⟶ (the Approve gate)              ├─error──▶ Scan failed ──retry──▶ (back to Scanning)         │
                                                                                        └─clean (no flags)───────────────────────────────────────▶ ┘
  Corpus ──⚙ edit later → re-parse──▶ Configuration          Corpus ◀──review banner──▶ Grouping review
  Corpus ──Enter──▶ Focus     Corpus ──L──▶ Loupe     Focus ⇄ Loupe (L / Enter)
  Series folio ──+ new series──▶ scoping ──build──▶ builder ──open a member → Focus
  top nav: Experiments ⇄ Series
```

Legend: solid = primary forward · dashed = branch/back/secondary · accent = into Focus/series. Every
screen also has a `‹` up-link + `Esc` (up one level).

### 4.2 The two phases (the load-bearing split)

- **Phase ① — index (cheap, read-only, no writes).** Runs on **Configuration**. A pure filesystem
  listing: `readdir` + glob-match against the image / metadata / integration patterns →
  `{total, matched-by-type, unmatched grouped by miss-type}` + a geometry preview derived from a
  capped sample of PRP/setup files. **It writes nothing.** Because it is pure listing it is fast even
  over SMB. Endpoint: `GET /api/fs/manifest`. Editing a pattern re-fetches (request-scoped,
  abortable). Approve sits disabled until it returns.

- **Phase ② — parse (the real scan).** Runs after Approve, on the combined Scan + grouping surface.
  `scan_and_group!` parses PRP/geometry, computes grouping, commits the **structure** in one fast
  transaction (queryable instantly), then analyzes **per-exposure outside the transaction**
  (peaks/phase commit row-by-row, streaming in). Editing Configuration later re-runs phase ②.

**Approve is the creation boundary and it is the existing `POST /api/experiments`.** That route
already: validates `isdir` → 409s on duplicate `data_dir` → inserts `ingest_status='scanning'` →
returns 202 → `Threads.@spawn scan_and_group!`. So there is **no draft DB row and no separate
approve route** — the pre-approval draft (root + resolved sub-paths + edited patterns + the phase-1
manifest) lives entirely **client-side** in `useDraftExperiment`. The one backend addition was:
accept the resolved `{name, data_dir, analysis_dir, patterns}` in the create body so Configuration's
edits flow straight into phase ②.

### 4.3 The directory-resolution model (the 2026-06-22 refinement)

**You point at ONE directory — the experiment root — and Himalaya structurally resolves everything
inside it.** This replaces the old "type every path by hand / `experiment.toml`" model. `.toml`
config is **deprecated going forward** (the funnel is toml-independent; removing the `config.jl`/toml
*code* is a separate future cleanup).

- **`GET /api/fs/resolve?path=…`** (the wire param is `path`, `routes_fs.jl`:167 / `api.ts`:273; its
  value is the experiment root) returns `{data_dir, analysis_dir, setup_file, setup_ambiguous,
  image_pattern, metadata_pattern, integration_pattern}` and the resolved `name`. *(Verify the handler
  serializes `:name` — one review pass found the response NamedTuple omits it, `routes_fs.jl`:170-177;
  if so, either add `:name => r.name` or have Configuration prefill the Name field from the root
  basename. The Name is editable regardless.)* It is **name-based and fast** (a couple of
  shallow `readdir`s + bounded stats, `routes_fs.jl`:110-136 — *not* a deep walk; the deep-walk version
  was unusably slow over SMB, the name-based rewrite is near-instant. Exact timings are environment-
  dependent; the contract is "no deep walk," not a millisecond figure).
- **Path display is root-relative (ratified 2026-06-22).** The experiment **root is absolute** (typed
  on the picker, shown under the Configuration title). The **data / analysis locations display and
  edit relative to that root** (`data`, `analysis/automatic_analysis/tot_files`). The Edit affordance
  shows a dim fixed `…/<rootBase>/` prefix + an editable suffix. A leading `/` in the suffix still
  accepts an absolute override (e.g. a shared analysis dir outside the root). **The backend always
  receives absolute paths** — the relative form is a display transform only.
- **Setup file is auto-resolved**, surfaced as an editable field **only when ambiguous**
  (`setup_ambiguous` = none found, OR setup files span more than one directory; multiple files in a
  single dir auto-pick the latest by name, `routes_fs.jl`:130-132). Otherwise a quiet caption: shipped
  `Geometry read from <basename>` (`ConfigurationPage.tsx`:489-491); the ambiguous-branch copy is at
  `:473-477`. (Caption copy is editable; the rule is "name the geometry source," not a literal string.)
- **Integration-layout auto-detection (the SSRL convention).** `_detect_integration_layout` probes a
  sample TIF stem against the analysis subtree and recognizes the canonical SSRL layout, returning
  both the directory and the patterns:
  - **TOT convention** — `analysis/automatic_analysis/tot_files/{base}_tot.dat` (per-base totals,
    frame-stripped). Keys patterns off the actual frame suffix:
    `image={name}_0_001.tif` · `metadata={name}_0_001.prp` · `integration={name}_tot.dat`.
  - **Per-frame fallback** — `{stem}.dat` → `{name}.*`.
  - else fall back to `root/analysis`.

  This is **pattern-driven, no stem-equality assumption**: with `image={name}_0_001.tif` the stem
  *is* the base, so `{name}_tot.dat` resolves to `{base}_tot.dat` directly. Per-base falls out
  naturally; extra calibration frames (`_0_002`…`_0_005`) are excluded by the `_0_001` image pattern
  (their data is already in the base total). The resolver seeds these patterns into the draft;
  Configuration shows them pre-filled and editable. **mini and every other SSRL dataset resolve with
  zero hand-wiring.**

---

## 5. Surface specs — transcribed from the mockups

Severity for the *gap to live*: **P0** core daily loop · **P1** spec'd surface wrong/missing ·
**P2** degraded · **P3** nit.

### 5.1 Experiments gallery — `/experiments` · `p1-gallery.png`

Header: kicker `Experiments` (terracotta, +4px below to the title) · serif display `All experiments`
· a terracotta **`+ New experiment`** button right-aligned.

Body = a **timeline rail + year-grouped card grid**:
- **Left rail** (slim, ~120px): heading `TIMELINE`; one row per year (serif year + mono `N sessions`),
  a vertical hairline with dots, the active year's dot filled terracotta. (Kept deliberately slim —
  drop it if it reads as noise.) *Typography is a fidelity target: shipped renders the year sans
  (`text-body font-medium`), the count `text-ink-soft`, all dots `neutral`, and no hairline rail
  (`ExperimentsHomePage.tsx`:132-147). Build toward serif year / mono count / terracotta active-year
  dot / hairline-with-dots.* **Rail count semantics (target):** `N sessions` should be the **macro-
  session sum** for the year (consistent with the masthead's `sessions` meaning, §5.6); shipped uses
  `items.length` = experiment-count-in-year (`:141`), a divergence to fix.
- **Year sections** (newest first): a serif year head with a bottom hairline, then a 2-up card grid.
- **Card** (`xc`): the name (`SSRL April 2026 · 1p7m`, **middot not hyphen**) + a **state chip**
  top-right + a **date-range line** (below) + a counts line (`170 samples · 682 exp · 13 loads`). **No
  `data_dir`** on the card (ratified, Jonathan 2026-06-22 — the path is noise here; it lives in the
  corpus header once you open the experiment). *Typography note: the mockup shows a **serif** name + **bold**
  count numbers; shipped renders the name sans (`text-body font-medium`, `ExperimentsHomePage.tsx`:180)
  and counts plain `text-ink` (`:201-206`) — a fidelity gap to close toward the mockup (serif name, per
  §2's serif-for-names rule).* State chips:
  - `up to date` (success, outline)
  - `N to review` (warning, outline) — the flagged-sample count (`review_count`, a **flat list-only**
    field; absent from the detail endpoint — §9)
  - `scanning…` (accent, outline) — the card body becomes `indexing N files…` (target; shipped is
    `indexing…` with no count, `ExperimentsHomePage.tsx`:198, shown on `isScanning || !stats`)
- The whole card is a **`Card as="button"`** (the door pattern — not a clickable `<div>`). Click →
  that experiment's `/corpus`.

**Empty state** (an `<EmptyState title="No experiments yet">`, `ExperimentsHomePage.tsx`:111-127 — no
`experiments-empty` testid; the only testid is `new-experiment-cta` at `:103`): the header
`+ New experiment` is **disabled** (redundant with the banner); a centered empty banner with a
**50%-larger** primary CTA (`Button size="lg"`, a
design-system size, not a scale hack), measure capped `max-w-[34ch]`, copy: *"Point Himalaya to your
experiment directory. It will discover your setup, scan your exposures, group them into samples.
Easy."*

**Backend prerequisite (§9):** the list endpoint must carry per-card stats + `started_at` for year
grouping + `review_count` (the nested-`stats` + flat-`review_count` shape is detailed in §9).

**Card date range (ratified, Jonathan 2026-06-22).** Add a quiet **date-range** line to the card: the
experiment's acquisition span, **month + day, year OMITTED** (the year is already the section heading).
Format `Mon D–D` same-month (`April 12–14`) or `Mon D – Mon D` cross-month (`April 28 – May 2`). Derive
it client-side with **no new backend field**: start = `stats.started_at` (min exposure ts), end =
`started_at + stats.span_hours` (= max exposure ts). Place it as a `text-caption` line under the name.
*(The earlier "third name clause / lipid mesophase series" was from the **stale**
`experiments-home.html` / `ingestion-prototype.html`, not the canonical p1-gallery — there is no
description line and none is wanted. Final card = **name · state chip · date range · counts** (no dir,
no description).)*

### 5.2 New-experiment picker — `/experiments/new` · `p1-new.png`

Kicker `New experiment` · serif display `Point at a directory` · intro (`max-w-[60ch]`): *"Pick the
experiment folder (the one holding the data and analysis directories). The next step finds the pieces
and lets you review and correct them before anything is created."* (shipped copy,
`NewExperimentPage.tsx`:115-119.)

A **full-width Directory card**:
- a mono path input (the active path, blinking caret), placeholder reads clearly as a hint
  (`Type or paste your experiment directory…`), **wide** (it is the primary input, long paths must
  fit);
- below it, a **live autocomplete dropdown** of matching directories (mono, hairline-separated rows);
- an assistance hint line: *"Start typing and we suggest matches. Tab completes, ↑↓ choose, ↵
  confirms."*;
- the **two pre-flight checks**, right-aligned (green ✓ / red ✗ / faint ◦ pending): `✓ directory
  exists` · `✓ not already an experiment`. The dup check is computed **locally** from
  `useExperiments().data` (`NewExperimentPage.tsx`:60-66), surfaced as the right-aligned ✗ pre-flight
  (`:147-155`) — **not** a line below the field (`DirectoryPickerField` is passed `validation={null}`,
  `:132`). It mirrors the create route's 409-on-dup; if taken, Review is blocked.

A sticky **action bar** (soft upward shadow, the §2 inline value): note *"Nothing is created yet. The
next step indexes and lets you review."* (shipped, `NewExperimentPage.tsx`:166-168 — em-dash-free per
§2's copy law) · spacer · `Cancel` (ghost) · `Review →` (accent, enabled only when both checks pass).

**No indexing here** — Review just commits the **root** to `useDraftExperiment.setRoot` and navigates
to `/experiments/new/config`.

### 5.3 Configuration — `/experiments/:id/config` (+ the `new/config` draft route) · `p1-config.png`

Two modes, one route.

**Header (both modes):** kicker `New experiment · Review before scan` (first-run, capital R,
`ConfigurationPage.tsx`:253) or the experiment kicker (later) · serif title = the **editable Name** (`Input variant="title"`, with a "click to
rename" cue so the affordance is discoverable) · the absolute **root path** in mono below · a
**task-framing paragraph on the right** describing what to do here (review the auto-resolved layout +
geometry, correct anything, Approve) · an **undo/redo** control top-right. *(Wiring status: the
mockup shows undo/redo; `ConfigurationPage` does not yet wire `useShortcuts` for `⌘Z`/`⌘⇧Z` — this
is a §10 Remaining item, not shipped.)*

**Phase-1 indexing line:** `Indexing directory…` + a mono **`N exposures found`** readout (the matched
image count — **no `/total` denominator**, per §6's honesty fix; drop the earlier `· M seen`). Shipped
matches this (`ConfigurationPage.tsx`:283-294, a placeholder `ProgressBar value=0 total=1`).

**Two-card grid:**

- **Sources — file patterns** (editable). Three fields, **capitalized + renamed**: **Exposure
  pattern** · **Metadata pattern** · **Integration pattern** (each shown as the resolved glob, e.g.
  `{name}_0_001.tif`, `{name}_0_001.prp`, `{name}_tot.dat`). The **label typography is distinct from
  hint typography** (`.text-meta` ink labels vs `.text-caption` ink-soft hints). Below: *"Editing a
  pattern re-runs the index."* (shipped, `ConfigurationPage.tsx`:396-398; the "undo/redo applies to
  every field" clause is a §10 target, not shipped). Plus a **per-pattern coverage line** (§6).
  - The **Data directory** is shown **read-only** (it is the root the user already picked — review,
    not type). The **Analysis directory** is shown **read-only with an Edit reveal** exposing a
    root-relative-suffix field. **Target: that field is autocomplete-assisted** — it suggests
    subdirectories under the root (the `DirectoryPickerField` affordance the user asked for in the
    walkthrough: "bring the picker back, don't make analysis dir a bare text field"), so the user is
    guided, not typing blind. *Shipped today is a plain free-text suffix `Input`
    (`ConfigurationPage.tsx`:340-351, no `DirectoryPickerField`) — that is the gap to close, not the
    target.* The **Setup file** field appears **only when ambiguous**.

- **Geometry — auto-derived** (read-only key/value rows). First-run is **flex `justify-between` rows
  in a `divide-y` column** (shipped, `ConfigurationPage.tsx`:405-451); the `104px | value | source-tag`
  **aligned grid** belongs to the **later-⚙ edit mode** (the "one aligned key/value table" below).
  Rows: `beam center 421.4, 836.9 px` `[setup]` · `flight path 1.8095 m` `[setup]` · `pixel size
  172 µm` `[prp]` · `energy 9.0 keV` `[prp]`. Source tags are **quiet gray mono chips** — only `prp`
  and `setup` are emitted today (§11.2; `computed`/`human` are not yet produced). Caption (shipped,
  `ConfigurationPage.tsx`:493-495): *"Derived automatically: prp = file sidecar, setup = beamline log.
  You can override later."* (Per-field override is a later-⚙ affordance; first-run is review-only.)

**Latest files scanned card:** a 2-column mono list of the most-recent matches + a `N matched`
count, so you can sanity-check the right files are arriving. Calibration frames are tagged
(`agbe_… · calibration`).

**Action bar:** note *"Approve creates the experiment and starts the full scan, enabled when indexing
finishes."* · spacer · `Cancel` (ghost) · **`Approve`** (accent; shipped is a **static** "Approve",
disabled while `indexing || !manifest || createMutation.isPending`, `ConfigurationPage.tsx`:543-547 —
the mockup's dynamic `Indexing… N` label, `N = manifest.matched.image`, is a §10 fidelity target, not
shipped). Approve = `POST /api/experiments` with the resolved values → consume `.id` from the 202
partial → navigate to `/experiments/:id/grouping`. **Approve error path (gap to build):** the create
route can 409 on a duplicate `data_dir`; the shipped `createMutation` has **no `onError`**
(`ConfigurationPage.tsx`:153-163), so a 409 currently does nothing. Target: surface the error as a
toast and keep the user on Configuration (do not navigate).

**Later ⚙ edit mode (post-ingest):** same layout, but geometry/sources become **one aligned
key/value table with per-field override**; saving a pattern change re-runs the **gated** phase-②
rescan (§6.6). Geometry is fully populated (phase ② has run).

### 5.4 Scan + grouping review (combined) — `/experiments/:id/grouping` · `p1-grouping.png`

**One surface for first-ingest scanning AND the standalone rescan-review route.** Full-screen
takeover, its own single header (§3.2). *(The grouping components — `GroupingReviewPage`, `LoadFold`,
`SampleFold`, `GroupingBulkBar` — live under `print/components/`, not `print/pages/` where
`ConfigurationPage`/`ExperimentCorpusPage`/`NewExperimentPage`/`ScanFailedPage` live.)*

**Scanning mode** (`inFlight.status==="scanning"`, with an `exp.ingest_status` fallback for
reload-mid-scan):
- kicker `Scanning · review groups as they land` · serif experiment name.
- A **live progress line**: spinner + `Parsing exposures…` + mono `418 / ~682 · 8 loads · 102
  samples` + a progress bar + an **amber** flag count (the mockup is the target; amber is the
  app-wide flag color). *Shipped currently colors it terracotta (`text-print-accent`,
  `GroupingReviewPage.tsx`:295) — a fidelity gap to amber, not an accepted exception.*
- **Loads unfold live**, newest activity visible:
  - **Flagged completed load** (warning border, expanded): a header (`load 06` · `12 samples` · a
    warning pill — mockup `2 need review`, shipped `⚠ N to check`), then one row per flagged slot —
    sample name (`2-2 + LL37 1:1`) · the **flag reason** · inline verbs right-aligned. *Copy: the
    mockup is the design target (`stage jump → maybe same sample as slot 4`, `counter reset mid-load`,
    a `slot N` chip); shipped strings are `position jumps X → Y mm here, likely two samples`
    (split divider, `SampleFold.tsx`:143-159) / `Filename matches <label>. This looks like a reshoot in
    a later load.` (merge prompt, `:116-138`), verbs `Keep separate` / `Merge into that sample` /
    `Split here` / `Split…`, and shipped rows have **no** `slot N` chip.
    Build toward the mockup copy where it reads better; the verbs are already aligned.*
  - **Clean completed load** (collapsed to one line): `load 07` · `10 samples` · a clean marker
    (mockup `✓ clean`, shipped `✓ grouped cleanly`, `LoadFold.tsx`:61) · `▸ expand`.
  - **Still-loading**: a single generic `unfolding…` tail block (`GroupingReviewPage.tsx`:377-384) —
    one placeholder, not one per in-flight load.
- **Action bar:** note *"Review flags as loads land. Confirm unlocks when the scan finishes. Later
  loads can still raise flags."* (period, not em dash, per §2's copy law) · spacer · **`Confirm
  groups`** (accent, **disabled** showing `Confirm groups · scanning…`, `GroupingReviewPage.tsx`:500,
  until `ingest_status==='complete'`).

**Mechanism (no transaction surgery):** "live" rides the existing rail. The structure commits fast +
is queryable immediately; per-exposure analysis streams (already outside the txn). The SSE
`ingest_progress` frames already invalidate the loads cache each tick, so loads unfold with **no new
plumbing**. Safe because ingest is **additive** (a settled load is never re-touched) and edits go
through the **event log** exactly as post-ingest. The structural-edit verbs (Keep separate / Merge /
Split / dismiss, with free undo) are the existing event kinds. `Confirm groups` is the terminal
commit → `/corpus`. A **clean scan (no flags)** can skip straight to `/corpus`.

Post-scan (or the `/grouping` rescan-review route), the same surface shows the **static** review:
the 3-level fold (load → slot → exposure leaf on expand), the split divider, Rename/Split, and the
filter/search (hidden while scanning).

### 5.5 Scan failed — `/experiments/:id/corpus` (failed state) · `p1-failed.png`

> **Spec-ahead-of-code.** This whole surface is the design target from `p1-failed`; the shipped
> `ScanFailedPage.tsx` is a stub (plain `Scan incomplete` in `text-ink`, no warn tone, `:93`; a bare
> stem list, no nearest-file pairing; plain Inputs, no live `✓N/N`; confirm copy
> `Ingest N parsed files?` / `Cancel` / `Confirm`, armed branch `~:160-189`; no `↵` chip). The copy below is the
> **intended** target — it is a §10 Remaining item, not built. Where shipped and mockup disagree the
> mockup wins; the confirm copy to build is the two-stage form below.

Kicker `Scan incomplete` (warning) · serif name · a sentence: *"**658 of 682 parsed.** 24 couldn't
be paired, mostly missing integration traces, a few metadata mismatches. Each row shows the nearest
file that does exist, so you can see what to match."*

> **Correctness requirement (B1):** the failed-state manifest fetch (`ExperimentCorpusPage.tsx`, the
> `fetchManifest` call ~`:90-94`) **must pass the resolved leaf `analysis_dir` and `setup_file`**, exactly
> as Configuration phase-1 does. `fetchManifest` only matches integration `.dat` against `analysis_dir`
> when supplied, else falls back to `data_dir` (`routes_fs.jl`:208-210) where `.dat` never lives — so
> every integration trace reports `unmatched` and the pattern-test can never clear. The shipped call
> passes neither today and must be fixed (same bug class as §6).

**Two-card layout:**
- **Didn't parse — N files** (scrollable, `max-h` with `scroll ↓` hint): unmatched files **grouped by
  miss type**, each row = the unmatched file + the **nearest existing file** with the discrepancy
  highlighted in warning (`nearest HA_5_044_S1991_total.dat`, `nearest RY_BL3_S1780.PRP (vs *.prp)`).
  The list **scrolls** so a fix can be verified against every miss.
  - *Backend miss types (ratified — build it):* the manifest route emits only `metadata` and
    `integration` miss labels today (`routes_fs.jl`:230; `MissType = "metadata" | "integration"`,
    `ScanFailedPage.tsx`:18). The mockup's third group **`Image mismatch`** (e.g. `.tiff` vs `*.tif`,
    `.PRP` vs `*.prp`) is a **required** part of this surface, not optional polish: near-miss
    extensions are a leading cause of a half-failed scan, and a diagnostic screen must see them.
    **Decision (Jonathan, 2026-06-22): do it right** — widen the manifest to detect near-miss image
    extensions and emit an `image` miss label, extend `MissType` (§10 Remaining). The current
    two-group state is a stopgap, not the target.
- **Test patterns** (adaptive): one field **per file type that has misses** (integration / metadata,
  plus image once the backend emits it), each showing a live `✓ N / N` as you edit until it clears. A
  full-width **`Apply all in Configuration →`** pre-fills Configuration.

**Action bar:** note *"A failed scan is usually a configuration issue."* · spacer · an **in-place
two-stage confirm** for partial ingest (`Ingest 658, skip 24?` → `No` / `Yes, ingest 658`) · primary
**`Open configuration ↵`** (accent, frosted ↵).

### 5.6 Corpus assembly — `/experiments/:id/corpus` (clean/has-flags state) · `p2-corpus.png`

The daily-loop home. **No tabs.**

> **Ownership:** the experiment header and the 5-stat masthead are rendered by **`ExperimentShell.tsx`**
> (the layout route wrapping `/experiments/:id/corpus`; header `:85-149`, StatBar masthead `:62-77`).
> **`ExperimentCorpusPage.tsx` renders only** the review banner + sheet + transient bars + dock — it
> must NOT render the header/masthead or they double-stack with the shell.

**Experiment header:** kicker `Experiment` · serif name (`SSRL April 2026 · 1p7m`) · right-aligned
**rescan status** + a `Rescan` button + the **⚙** gear · the absolute `data_dir` in mono below.
*Target: `● scanned 2 min ago · up to date` (success-green dot + relative time). Shipped renders a
plain `text-ink-soft` status word (`Up to date` / `Scanning…` / `Analyzing…` / `Scan failed`,
`ExperimentShell.tsx`:119-131) — no dot, no timestamp. The dot + `N min ago` are §10 fidelity targets;
`N min ago` needs a new backend `last_scanned_at` field (the experiment payload carries none today).*

**5-stat masthead** (vertical-rule separated): `13 loads` · `170 samples` · `682 exposures` ·
`30.8 h span` · `4 sessions`. (`span` = max−min exposure timestamp; `sessions` = macro-session count,
populated at ingest via `_assign_sessions`, >3h gap → new session.)

**Review banner** (amber, when flagged): a warning triangle + *"**3 samples need a grouping check.**
Two stage-position jumps and one counter reset were ambiguous."* + a right-aligned terracotta
**`Review grouping →`** (to `/grouping`). The breakdown sentence is built from the loads-rollup
split/merge flag counts. Uses `bg-warning/10 border-warning text-warning` tokens.

**The sheet** (`SheetTable`, scoped to this experiment): columns
`☐ · SAMPLE · EXPOSURES · KEPT · TAGS · PHASE` (`SampleTableRow.tsx`:76, render `:178-228` — there is
no `(open)` door cell). A row =
- a select checkbox;
- **identity cell**: serif sample name + a mono id (`S2453`) + a **`slot N` chip** (joined from the
  loads rollup `slot_index` by `sample_id`);
- **exposures cell**: a strip of **detector thumbnails** (`--frame-edge` bg) — the current frame
  outlined terracotta, dropped frames at 40% opacity, a `+N` overflow;
- **kept cell** (visible header "Frames kept", `SheetTable.tsx`:192): `3 / 8` (mono, kept bold ink,
  denominator faint);
- **phase cell** (visible header "Phase"): renders the **StatusCell composite**
  (`SampleTableRow.tsx`:224), not a bare chip — a phase chip (`Pn3m` / `Im3m` per-phase hue) OR
  `Form factor` OR faint `not indexed` OR `No exposures`, depending on sample state;
- the current row is terracotta-tinted (`--accent 10%`).

There is **no per-row "Needs review" pill** (that was noise from never passing `screened`; honest
per-sample `screened` drives any indicator now). There is **no per-row Focus door** — Focus is a
dedicated **Dock** action + the `Enter` key. (The `Tags` column is **kept** per Jonathan, beyond the
mockup.)

Below it all: the **Dock** (§7, Corpus grammar).

**Empty / honest states:** an experiment with no grouping flags shows the sheet with no banner; an
experiment mid-scan redirects to `/grouping`; a crashed-mid-scan `ingest_status='scanning'` should
show a scanning surface (a known edge — §10 G4).

---

## 6. Counting & coverage (the 2026-06-22 honesty fix)

The old header read `matched.image / total` = "142/284 indexed" — apples-to-oranges (the denominator
counted sidecars, not exposures) and read like half the files failed. And integration `.dat` was
matched against `data_dir`, where `.dat` never lives, so it was silently always 0 and the 142
"missing integration" entries were computed but never shown.

**The model now:**
- **Headline = "N exposures found"** (the matched image count, the real exposure count). No `/total`
  denominator.
- **Per-pattern coverage line:** `Image 138 · Metadata 138/138 · Integration 138/138` — each type
  matched against the **right directory** (integration against the **analysis subtree**, mirroring
  `grouping.jl scan_directory` → `sidecar(analysis_dir, stem, dat_pattern)`; image/metadata against
  `data_dir`). A type under 100% renders **amber**.
- **Explicit warning** (`role="status"`, amber `text-warning`) when a pattern matches under its
  exposure count: *"Integration matched 0 of 142. Check the integration pattern or the analysis
  directory."* — a truthful, actionable signal, not a hidden zero.

`GET /api/fs/manifest` gained an `analysis_dir` param; integration matches against it (flat join),
falling back to `data_dir` when absent (backward-compatible). **The `analysis_dir` must be the
resolved LEAF integration directory** (what `_detect_integration_layout` returns, `routes_fs.jl`:67-73),
because the match is a flat `readdir` (`:209-210`), not a recursive walk — passing a non-leaf root
yields 0 integration coverage.

**Zero of an essential type, first-run (gap to build):** when **image, metadata, OR integration**
matches **0** everywhere (pattern wrong, per the zero-coverage model above), shipped suppresses the
per-pattern shortfall warning (`exp > 0` guard, `ConfigurationPage.tsx`:93) but still renders `/0`
denominators and leaves Approve enabled. Target: show a single clear line naming the missing type
(*"No exposures matched. Check the exposure pattern."* / *"No metadata matched…"* / *"No integration
traces matched…"*), hide the `/0` denominators, and **disable Approve** — you cannot create an
experiment with zero complete exposures.

**Zero-coverage model (ratified, Jonathan 2026-06-22).** An exposure is the triple (image, metadata,
integration), and **all three are essential "for now"** — there is no experiment without integrations,
and metadata (`.prp`) is required too (it carries the per-exposure timestamp the grouping segmentation
needs, §12.1, plus geometry). *("For now" leaves room to later add a fallback, e.g. file mtime for the
timestamp, which could make metadata optional; until then it is hard-required.)* So:
- **A whole type matches ZERO everywhere** (image, metadata, or integration; pattern wrong or data
  absent) → **hard block**: amber, **disable Approve**, `Check the <type> pattern.` There is **no
  "advisory / informational" tier** — zero of any leg means the exposure triple cannot be formed and
  the experiment cannot be ingested. (This retro-justifies mini before `tot_files`: integration 0 was
  un-ingestable, correctly amber, and Approve should have been blocked.)
- **A type matches PARTIALLY** (some stems complete, some missing one leg) → amber-actionable but **not
  blocking**: the incomplete stems flow to the scan-failed "didn't parse" list and the
  ingest-what-parsed path (§5.5).

### 6.6 Corpus state machine + rescan

The corpus route is **state-driven** (one route, body switches). Discriminator derived at
`ExperimentCorpusPage.tsx`:73-78; the early-returns are at `:339-366` — `scanning` → `<Navigate>` to
`/grouping` (`:339-343`), `rescanning` (`analyzing`) → inline `ProgressBar` (`:345-355`), `failed` →
`ScanFailedPage` (`:357-365`); otherwise the sheet (with/without the flag banner).

| State | Body |
|---|---|
| scanning (first ingest) | redirect to `/grouping` (the combined surface) |
| rescanning (`analyzing`) | progress **inline** on Corpus (a ProgressBar), never the full-screen surface |
| scan failed | the Scan-failed surface (§5.5) |
| has grouping flags | the sheet + the review banner |
| clean | the sheet |

**Rescan runs in place on Corpus** — never back to the funnel:
1. phase-1 cheap check (`cheap_change_check`, on-disk count vs persisted exposure count);
2. nothing new → silent "up to date";
3. new files → **additive** phase-② parse of the delta only (existing samples never re-parsed or
   re-grouped → curation preserved); the sheet grows in place; new flags raise the banner;
4. errors → an inline notice/toast, not the full-screen page.

The **only gated rescan** is a Configuration / **pattern** change (via ⚙), because patterns can
re-include/exclude files. A pattern-edit PATCH chains a **forced** rescan; geometry / non-pattern
edits do not. (The auto-scheduler `_rescan_tick!` drives the same lifecycle on its tiered backoff,
emitting the full `ingest_*` SSE frames with `phase="rescan"` so a rescan shows the inline analyzing
surface, not the takeover.)

---

## 7. The Dock — per-surface grammar

From `bars.html` (`b1`–`b4`). A persistent light bar. Grammar, left → right, **consistent on every
surface**:

```
‹ up-link  │  cursor (steppers w/ N/total readout)  │  cull (where frames exist)  ……  destinations
```

- **up-link** (left anchor, accent): one level up. Present on *every* surface so controls never
  shift. `‹ Experiments` (Corpus) · `‹ Corpus` (Focus, Loupe) · `‹ All series` (Series).
- **cursor**: labeled steppers mirroring the keyboard. **`Sample ↑ N / total ↓`** always; **`Frame ‹
  N / total ›`** where frames exist. The `lbl` ("Sample" / "Frame") recedes via `--ink-faint` at the
  same size; the `N / total` is mono (total faint).
- **cull**: verdict verbs, only where frames exist. Drop (`X`, accent outline) · Keep (`K`, success
  outline) · then per-surface (Restore / Set representative).
- **destinations** (right, primary rightmost): `Loupe` (`L`, neutral pill) · `Focus` (filled accent,
  the unambiguous primary). *The **frosted `↵` key-chip** on the Focus button is a design target — no
  shipped dock renders it yet (`ExperimentCorpusPage.tsx`:541-543, `LoupePage.tsx`:586-593 are plain
  `variant="accent"` "Focus" buttons); adding the `↵` KbKey on Corpus/Loupe/Series is a §10 item.*

**Per-surface table (authoritative — the keymap in §8 matches it exactly):**

| Surface | up-link | cursor | cull | destinations |
|---|---|---|---|---|
| **Corpus** | `‹ Experiments` | Sample ↑↓ · Frame ‹› | Drop · Keep · Restore | Loupe · **Focus** |
| **Focus** | `‹ Corpus` | Sample ↑↓ | — | Loupe |
| **Loupe** | `‹ Corpus` | Sample ↑↓ · Frame ‹› | Drop · Keep · **Set representative** | **Focus** |
| **Series** | `‹ All series` | Sample ↑↓ | — | **Focus** |

Notes: **Corpus has no representative concept** → verbs `{Drop, Keep, Restore}`, `R` unbound. **Loupe**
operates within one sample's frames → full `{Drop, Keep, Set representative}` + binds `R`. When a
selection exists the dock shows a `N frame(s)` count (terracotta square + count) before the verbs.
**Focus** is sample-level only, no cull (frame work lives in Loupe); its one destination `Loupe` is a
quiet neutral button (on Focus the indexing *is* the work). **Series** names the current sample in the
dock (swatch + name + `from <experiment>`) and highlights its trace on the overlay, so "which sample
Focus opens" is unambiguous; no cull (nothing to keep/reject in a series).

> **Series dock — design target, mostly not shipped.** The table above (Series: Sample `↑↓` stepper +
> `Focus` destination + the swatch/name/from-experiment segment) is the target. Shipped, the
> `SeriesBuilderPage` dock has **only the `‹ All series` up-link** (`SeriesBuilderPage.tsx`:536-544),
> and `SeriesScopingPage`/`SeriesFolioPage` render no dock at all. To build: bind the stepper to the
> member-list cursor, `Focus` to the current member's sample, and source the swatch/name/experiment
> from that member. (§10 Remaining.)

**Key-chip treatment** (`glyph.html`): a `currentColor` keycap, centered (flex + `line-height:1`).
**Frost** (translucent `currentColor` fill) is the treatment for a key-chip on a **filled colored
button** — i.e. the `Focus` primary's `↵` **where built** (the chip itself is a design target, not yet
on shipped docks — see destinations above). Neutral and colored-outline buttons get a plain chip whose
key matches the text color. Uniform text size across all dock controls; hierarchy is weight + color, never size. The
unit is **Sample** on every surface, Series included (members are samples).

The dock primitive (`print/ui/Dock.tsx`) is a plain flex row, placement-children only; each page
composes its own segments (no surface enum, no shell-mounted state). Light treatment throughout — the dark,
attention-grabbing look was the *transient* CullBar's, which it earned by appearing only on selection.

**Distinct from the persistent Dock, Corpus has two *transient* floating bars** shown only on
selection (`ExperimentCorpusPage.tsx`:462-478): the exposure-grain **`CullBar`** (Drop/Keep/Restore on
the frame selection) and the sample-grain **`ComposeBar`** (`+ New series` from a sample selection).
These are not folded into the Dock — keep all three.

*File locations (cited once): `Dock` and `ComposeBar` live in **`print/ui/`** (`print/ui/Dock.tsx`,
`print/ui/ComposeBar.tsx` — there is no `print/components/ComposeBar.tsx`); `CullBar` and the grouping
folds live in **`print/components/`**.*

---

## 8. Keyboard model (one unified set)

The whole app shares one keymap, designed holistically (existing bindings were **not** held sacred).
**This supersedes the 2026-06-13 keyboard-shortcut-library lock.** Four rules:

1. **Bare arrows navigate.** `↑/↓` = which **sample**; `←/→` = which **detail within it** (a *frame*
   on Corpus/Loupe, a *candidate* on Focus). First-class — active without tabbing into a grid.
2. **Bare letters = verbs** on the current sample/frame (mnemonic).
3. **`Enter` = drill into Focus** (the primary "go"); **`Esc` = up one level** (a ladder that
   dismisses the innermost transient state first).
4. **Modifiers consistent:** `Alt` = move/reorder · `⌘/Ctrl` = meta (undo, find, select-all,
   confirm) · `Shift` = extend/range · `Space` = toggle-select. Bare keys suppress while typing in a
   text field; **`?`** shows the keymap overlay.

**The map:**

| Key | Action | Where |
|---|---|---|
| `↑` / `↓` | Prev / next **sample** | Corpus, Focus, Loupe, Series |
| `←` / `→` | Prev / next **detail** — frame (Corpus, Loupe) · candidate (Focus) | |
| `Enter` | Open the current sample in **Focus** | Corpus, Loupe, Series |
| `L` | Open the current sample in **Loupe** | Corpus, Focus |
| `Esc` | Back / up one level (ladder) | all |
| `/`, `⌘K` | Find / jump (NavModal) | all |
| `?` | Keyboard-map overlay | all |
| `Space` | Toggle the current frame in the selection | **Corpus** |
| `Shift+←/→` · `⌘A` | Extend / select-all the frame selection | **Corpus** |
| `X` / `K` | Drop / Keep (selection on Corpus, else current frame) | Corpus, Loupe |
| `R` | Set current frame as representative | **Loupe only** |
| `Backspace` | Restore — clear the verdict | Corpus, Loupe |
| `P` | Toggle add-peak mode | Focus |
| `Enter` (Focus) | Apply focused candidate / toggle focused peak | Focus |
| `⌘Z` / `⌘⇧Z` | Undo / redo | Series, Configuration, Scoping (see status) |
| `Alt+↑` / `Alt+↓` | Reorder the sample up / down | Series |
| `A` | Add a sample | Series |
| `⌘Enter` | Confirm / build | Series |

**Ratified (Jonathan, 2026-06-22): the full keymap above is the design target — build it all, do not
trim it to current wiring.** Every binding is straightforwardly buildable (a few new `ShortcutId`s
plus per-page wiring); "not yet wired" below is a scheduling fact, not a design judgment.

**Wiring status (this is the DESIGN; not all of it is shipped — verified against `shortcuts.ts` +
the page `useShortcuts` calls):**
- **Wired:** Corpus `↑/↓` (`prevSample`/`nextSample`), `←/→` (`prevDetail`/`nextDetail`), `Enter`
  (open Focus), `L` (Loupe), `X`/`K` (Drop/Keep), `Backspace` (restore), `Space`/`Shift+arrow`
  selection, `/`+`⌘K` find, `?` overlay, `Esc`; Loupe `←/→`/`↑/↓` cursor + `X`/`K`/`R`/`Backspace`/Esc
  (`LoupePage.tsx`:352-367); Focus `←/→` candidate nav + `L` (`FocusPage.tsx`:487); Series Builder
  `⌘Z`/`⌘⇧Z` (`SeriesBuilderPage.tsx`:155).
- **Loupe scope (by design, not a gap):** Loupe culls the **current** frame only — no batch
  frame-selection model (the mockup `b4` dock shows the Frame stepper + current-frame verbs, no
  selection count). So `Space`/`Shift+arrow`/`⌘A` are **Corpus-only**, not Loupe. `Enter`→Focus on
  Loupe IS in the design (the Loupe dock's Focus destination) but is **not yet wired**
  (`LoupePage.tsx` binds no `openFocus`) — a §10 gap, not a scope narrowing.
- **NOT yet wired (DESIGN TARGET, §10 Remaining):** Series `↑/↓`, `Enter` (open Focus), `A` (add
  sample), `⌘Enter` (confirm/build) — none bound; `A` and `confirm`/`build` also have **no
  `ShortcutId`** in `shortcuts.ts`. Focus `P` (add-peak — `+Peak` is a button click, `FocusPage.tsx`
  ~`:749`; needs a new id) and Focus `Enter` (apply-candidate). Series Scoping has `⌘Z` only, **no
  `⌘⇧Z` redo**. Configuration wires neither undo nor redo.
  - *New registry defs to add (B3):* `addSample: { keys: ['a'], label: 'Add a sample', group: 'Edit' }`
    · `confirm: { keys: ['Mod+Enter'], label: 'Confirm / build', group: 'Edit' }` (`Mod+Enter` is
    currently unbound, no collision) · `addPeak: { keys: ['p'], label: 'Toggle add-peak mode', group:
    'Edit' }`. Keep these **omitted from `OVERLAY_IDS`** (like the other 8 surface verbs) — the overlay
    is a curated subset, not the full registry.
- **The `/grouping` surface is mouse-only (B4):** Keep separate / Merge / Split / dismiss-flag /
  Confirm-groups have **no keyboard bindings** (the bare-key model is for the sheet/Loupe/Focus, not
  the grouping takeover). Note the registry `dismiss` id (`shortcuts.ts`, Escape = "Back / dismiss") is
  **unrelated** to the grouping `grouping_flag_dismissed` verb (§9) — do **not** wire Escape to dismiss
  a flag.

**Changes vs. today:** drop `[`/`]` (→ `↑/↓` everywhere); candidate-nav moves `↑/↓ → ←/→` on Focus;
`Enter` opens Focus on Corpus; new `Space`/`Shift+arrow` selection + `⌘A` + `Backspace` restore + `?`
overlay. Focus's **exposure axis stays rail-scoped** (the existing `FO-EXPSKIP` filmstrip control,
*not* a bare-arrow gesture — bare arrows are one consistent sample/detail axis everywhere).

**Implementation invariants:**
- One registry (`src/print/shell/shortcuts.ts`); `matchShortcut` is a **flat global registry**,
  first-wins → **one semantic ID per physical key**, the *page* interprets it. `←/→` bind one
  `prevDetail`/`nextDetail` pair; `Enter` one `openFocus` (the real id, `shortcuts.ts`:53 — there is no
  `drillIn`); `Space` one `toggleSelect`. Never register two IDs on one key and disambiguate by surface.
- The **roving grid is removed** for the sheet. With the open-button column gone, the sheet's only
  targets are the sample row (`↑/↓`) and its frames (`←/→`) — nested 1-D cursors, not a 2-D grid.
  `lib/grid/useRovingGrid.ts` + `rovingGrid.ts` deleted, `SheetTable` `roving` prop +
  `SampleTableRow` context reverted to a plain table. The cursor is a page-level
  `{sampleIndex, frameIndex}`; pointer clicks set the same cursor (one source → mouse and keyboard
  never diverge); dock steppers call the same setter.
- The page handler **declines** (a) inside a text input / open modal — shipped, via
  `suppressGlobalKeys` (`useShortcuts.ts`:27-42 → `lib/keys.ts`), the structural `tagName`/
  `contenteditable` gate that also covers bare `Space`/`Enter`. Gate **(b)** — decline `Enter` on a
  native interactive target — is a **DESIGN TARGET, not implemented** (`suppressGlobalKeys` has no
  button/`a`/`[role=button]`/sort-header branch; `openFocus` guards only `activeSample == null`,
  `ExperimentCorpusPage.tsx`:249). To build: an early `return false` in each page's `openFocus`
  handler (or a shared helper keyed on `id==='openFocus'`) when `e.target` matches
  `button, a, [role=button], [aria-sort]`.
- **`Alt+arrow` is never shadowed by bare arrows — structurally, with no explicit guard.** `eventCombo`
  prefixes `Alt+` (`shortcuts.ts`:82-93) and `matchShortcut` requires an exact combo (`:96-102`), so a
  bare-arrow id can never match an Alt+arrow event. (There is no `altKey` ref on any page; don't spec
  one.)
- `?` overlay: special-case `e.key === '?'` in `eventCombo` to emit a stable `?` token regardless of
  the Shift bit (avoids the `Shift+/` ↔ `/` collision). The overlay (`KbdOverlay.tsx`) draws its
  glyphs from the registry but lists a **curated subset** (`OVERLAY_IDS`, 13 of the 21 registry ids).
  The 8 omitted are the registry ids `toggleSelect`, `extendPrev`, `extendNext`, `selectAll`, `undo`,
  `redo`, `reorderUp`, `reorderDown` (`KbdOverlay.tsx`:20-34). (Focus `P` / Series `A` / `⌘Enter` are
  not omitted-from-overlay — they are not registry ids *at all* yet; they must be added, see above.)
  *(Note: `KbdOverlay.tsx`:19 cites "spec §4.1"; in this canonical doc the keyboard model is §8.)*
- **Find and help are NOT page-registry dispatched.** The two app-global gestures — find (`/`, `⌘K`)
  and help (`?`) — are bound separately in `hooks/useGlobalShortcuts.ts`:30-51 (hand-rolled key checks,
  hoisted above the shell), not through `matchShortcut`/`useShortcuts`. "The page interprets it"
  applies to page-level registry ids only. (`shortcuts.ts` defines `find` `:64` / `helpOverlay` `:62`,
  but those are not the dispatch path.)
- `Esc` ladders per surface (Corpus: clear-selection → up; Focus: clear-preview → disarm +Peak →
  up). On the gallery/home (top of the ladder) `Esc` is a no-op.

---

## 9. Backend contract

**Already built (reuse):**
- `POST /api/experiments` — validate isdir (**400** if not a dir, `routes_experiments.jl`:171) →
  **409** on duplicate `data_dir` (`:181`) → create `ingest_status='scanning'` → async
  `scan_and_group!`. Accepts `{name, data_dir, analysis_dir, patterns}` in the body, where **`patterns`
  is nested `{ image?, metadata?, integration? }`** — the body keys are `image`/`metadata`/`integration`
  (`routes_experiments.jl`:200-201 reads `ppat(:image)` etc.; `api.ts`:229 types
  `CreateExperimentBody.patterns = { image?; metadata?; integration? }`). The `_pattern` suffix is the
  **DB column** name (`image_pattern`…) and the **separate** `/api/fs/manifest` query-param name — NOT
  the create-body key. An absent pattern stores NULL (`db.jl`:2173); the `{name}.*` default is applied
  **downstream in the scan layer**, not by the create route. `DELETE /api/experiments/{id}` exists
  (`routes_experiments.jl`:502; not surfaced by the funnel). **Returns
  `202` with only a partial `{id, status, name, data_dir}`** (`routes_experiments.jl`:253-255), NOT a
  full `Experiment` — consumers navigate on `created.id` and refetch for full state. (`createExperiment`
  is over-typed as `request<Experiment>` at `api.ts`:232; a `CreateExperimentResponse` partial is more
  honest.)
- `broadcast_progress!` + the four transient `ingest_{started,progress,complete,failed}` SSE frames
  (curation event name, **no** `user_actions` row) + the frontend `ingest_*` cache arms +
  `ingestInFlight`; `ingest_status` column; `phase="rescan"` discriminator.
- `scan_and_group!` — idempotent dedup-key ingest, fast structural txn + per-exposure analysis
  outside it, `on_progress` callback for ticks; `_assign_sessions` (macro-session segmentation at
  ingest); `cheap_change_check`; the tiered-backoff `_rescan_tick!` scheduler.
- Structural-edit events via the event log; `get_loads_rollup` + flags. **Event vocabulary:**
  `sample_renamed`, `exposure_moved`, `sample_created`, `sample_split`, `grouping_flag_dismissed`.
  **There is no first-class `sample_merged` event** (`applyRemoteToCache.ts`:403 confirms) — a merge is
  the frontend `mergeSamplesMutator` fanned out server-side as `exposure_moved` (per exposure) +
  `sample_renamed` (survivor) + `grouping_flag_dismissed`.
- **Own-op cache reconciliation is wired** (not a to-do): each mutator's `onSuccess` in
  `lib/queue/mutators/grouping.ts` calls `invalidateLoads` (`:37-40`; merge `:126`, rename/move `:161`,
  split `:196`, dismiss `:223`, undo-dismiss `:251`), invalidating **loads + samples + corpusSamples** —
  NOT in the `queries.ts` hook wrappers (`:846-866` pass no `onSuccess`). This is necessary because the
  own-op replay path runs `applyPostStateOnly` only (`replayCoordinator.ts`:47) and never calls
  `applyRemoteToCache` (foreign-only, `:87`), so the SSE arm alone would leave an optimistic negative-id
  split sample unreconciled.
- `GET /api/fs/resolve` (structural layout + pattern detection, §4.3); `GET /api/fs/manifest`
  (phase-1 listing + geometry preview + coverage, §4.2/§6); `GET /api/fs/suggest` / validate (picker).
- `DELETE /api/experiments/:id` (not surfaced by the funnel).

**Stat additions for the gallery/masthead** (the only net-new backend the surfaces needed):
- Each list row carries a **nested** `stats: { loads, samples, exposures, sessions, span_hours,
  started_at }` object (folded in via `_experiment_stats`, `routes_experiments.jl`:79-99, set as
  `d[:stats]` at `:151`) PLUS a flat sibling **`review_count`** (added by the LIST handler only,
  `_experiment_review_count` at `:111`, attached at `:264`). The frontend reads
  `e.stats.loads/.samples/.exposures/.sessions/.span_hours/.started_at` (`ExperimentShell.tsx`:62
  `exp.data?.stats`, `ExperimentsHomePage.tsx`:163,201-205). NOT flat `load_count`/etc.
  - Mechanism: `_experiment_row_to_json` folds in `stats` only; `review_count` is added separately by
    the list handler, so the **detail** endpoint `GET /api/experiments/:id` carries `stats` but **not**
    `review_count` (`routes_experiments.jl`:270-279).
- `span_hours` (in `stats`) + `session_id` populated at source (`_assign_sessions` +
  `backfill_load_sessions!` at migrate time so migrated DBs get sessions without a rescan).

---

## 10. Implementation status (what's done vs. what remains)

This branch has implemented the bulk of the design through milestones M0–M5 (all committed,
render-verified live, branch **stays unmerged** awaiting Jonathan).

**Done & committed:**
- **Shell + nav:** unified `TopNav`, routing collapse, `CorpusShell`/`CorpusTopbar`/`ExperimentTopNav`
  retired, `SamplesPage` retired (cross-cutting e2e sweep), flat Focus/Loupe routes.
- **Keyboard + sheet:** re-axised `shortcuts.ts`, roving grid removed, page-level cursor, `?` overlay.
  **Corpus is fully wired; Loupe `R` and Focus `←/→`+`L` wired. NOT yet wired (see §8 status): Series
  `↑/↓`/`Enter`/`A`/`⌘Enter`, Focus `P`/`Enter`, Config undo/redo, Scoping redo.**
- **Net-new components (built & committed, were mis-listed as net-new in an earlier draft):**
  `AcquisitionTimeline` (`b5378f34`), `useUndoStack` (`47a328b5`), `Dropdown`/`Menu` (`8e35b478`),
  `Toast` action slot, `PageFrame` `home`/`experiment` width keys. Only **`CullBar` parameterization**
  remains (§13).
- **Dock:** `ui/Dock` on Corpus / Focus / Loupe with the b1–b4 grammar (labeled steppers + N/total
  readouts + right-anchored destinations, Drop/Keep colored outlines). **Series dock ships only the
  `‹ All series` up-link** (stepper/Focus/sample-segment are §7 design targets, Remaining); the
  **frosted `↵` chip** on Focus buttons is not yet rendered anywhere (§7).
- **Backend stats (M1):** list-endpoint stats, `span_hours`, sessions-at-source.
- **Gallery (M2):** year-grouped card gallery + timeline rail + state chips + counts; empty state.
- **Corpus IA (M3a):** tabs retired → ⚙ gear, 5-stat masthead.
- **Corpus body (M3b/c):** real sheet (thumbnails / kept / phase / slot-chip), cull + Dock, amber
  explanatory banner, no per-row door, Tags kept.
- **Funnel (M4):** picker → p1-new (card + dual checks + footer); first-run Configuration → p1-config
  (geometry streams in before Approve); the **structural resolver** (`/api/fs/resolve`,
  experiment-root → autopopulated sub-paths, editable, root-relative display).
- **Grouping live-unfold (M5):** combined scan+review surface, live unfold via existing SSE, Confirm
  gate, takeover route.
- **2026-06-22 walkthrough refinements (UNCOMMITTED on branch, awaiting review):** notes 1–10 (gallery
  empty-state, picker layout + assistance, Configuration labels/read-only/guidance, path-relative
  model); finding 11 (coverage counting + integration-vs-analysis + warning); finding 12 (tot_files
  copied into mini + the SSRL integration-layout auto-detection in the resolver).

**Remaining:**
- **M6 polish:** the later-⚙ Configuration edit mode (per-field geometry override + aligned
  table + gated rescan-on-save + undo/redo wiring); dock polish sweep; Configuration fidelity vs
  `p1-config`.
- **Scan-failed surface (stub today):** build `p1-failed` fidelity — fetch the real manifest, nearest-
  file pairing, adaptive pattern-test card with live `✓N/N`, warn-toned kicker, the two-stage in-place
  partial-ingest confirm, the `↵` chip (§5.5).
- **Keyboard gaps (§8 status):** Series `↑/↓`/`Enter`/`A`/`⌘Enter` (needs new `addSample`/`confirm`
  registry ids), Focus `P` (needs `addPeak`) + `Enter`-apply, Configuration undo/redo, Series Scoping
  redo.
- **Backend miss types:** widen the manifest to emit an `image` miss label + extend `MissType` so the
  scan-failed "Image mismatch" group can render (§5.5).
- **Dock fidelity (§7):** the Series dock (stepper + Focus + swatch/name/from-experiment; only the
  up-link ships), the frosted `↵` key-chip on Focus destinations (Corpus/Loupe/Series), and the
  `Enter`-on-native-interactive decline gate (§8 gate b).
- **`CullBar` parameterization** (§13) — the one genuinely-outstanding component item. Compose rule:
  when the optional `actions: {label,onClick,variant}[]` is provided, render it **in place of** the
  fixed Drop/Keep/Restore/Clear set; otherwise keep the fixed handlers. `variant` ∈
  {`accent`,`success`,`ghostInverse`}.
- **The heavy end-to-end live walk:** Approve → full phase-② scan over real `tot_files` integration
  data (138 exposures) → grouping unfold → corpus render. The funnel/preview is verified; a real
  ingest has not yet been walked.
- **Deferred / open:** the gallery card third-clause meta (§5.1); the zero-everywhere coverage
  warning behavior (§6); series `series_samples.position`-gap compaction after merge (§11.4 — shipped
  leaves gaps; decide intended-final vs follow-up); the cross-version SSE deploy sequencing (§11.5);
  crash-mid-scan `ingest_status='scanning'` corpus state (G4); the toml/`config.jl` code removal (toml
  deprecated going forward); a single e2e tying an SSE frame to the corpus ProgressBar.

**Standing constraints:** branch stays unmerged; never `git add -A`; verify by rendering live;
refined surface layouts (Focus/Loupe/Series/Configuration body) stay untouched except for keyboard +
chrome wiring.

---

## 11. Data model & migration

This layer was developed in the 2026-06-15 ingestion design and the 2026-06-18 grouping session. It
is load-bearing for the funnel and grouping surfaces and is consolidated here (the source docs are
deleted).

> **Status: implemented & committed on-branch** (anchors): `loads` table (`db.jl`:78);
> `migrate_exposures_experiment_id!`; `migrate_samples_name_collapse!`; merge dedup
> (`routes_grouping.jl`:150-206); `grouping_flag_dismissed` + undo (`routes_grouping.jl`:333-387);
> `_assign_sessions` (3 h, `grouping.jl`:522); frontend op-queue `SCHEMA_VERSION=5` (`persistence.ts`).
> Genuinely-remaining: the cross-version SSE-frame deploy sequencing (§11.5) and the
> `series_samples.position`-gap compaction decision (§11.4).

### 11.1 The mental model

- **An experiment is a soft lens, not a hard partition.** Navigation is experiment-first (you enter
  through an experiment), but only the **corpus** is experiment-scoped. Samples, Series, Focus, and
  Loupe are *not* architecturally walled by experiment — which is exactly why Focus/Loupe are **flat
  global routes** (§3.2) and Series spans experiments. An experiment is "a kind of folder," a user
  concept over a directory.
- **Load ▸ Sample ▸ Exposure is the hierarchy.** A *load* is one acquisition sweep (a rack run); a
  *sample* is a slot within it; an *exposure* is one detector frame. Samples are **pure containers** —
  all analysis (peaks, indices, phase, assignments) is keyed to `exposure_id`. **Regrouping is pure
  bookkeeping: it triggers NO reanalysis** (moving an exposure between samples changes a foreign key,
  nothing recomputes). Structural edits are constrained to within one experiment.
- **No sample `role`/`type` field.** References (agbe calibration, empty capillary) are **ordinary
  samples named accordingly**, not a first-class classification. Ingest everything; the human names
  and curates. `agbe_*` frames are excluded from per-base totals only by the image pattern, not by a
  role flag.

### 11.2 Schema additions (additive, no destructive rewrite)

- **`loads` table (new first-class entity):** `id, experiment_id, load_index, session_id, start_time,
  end_time, frame_count, note`. Makes the Load ▸ Sample ▸ Exposure fold queryable. `session_id`
  populated at ingest by `_assign_sessions` (macro-session segmentation; >3 h gap → new session).
- **`exposures`:** add a **denormalized `experiment_id`** (because `sample_id` is ephemeral during
  grouping and `analyze_exposure!` needs the experiment). The **dedup key moves to
  `(experiment_id, filename)`**.
  - *Backfill = fail-fast:* migration **aborts and lists** any exposure with no derivable
    `experiment_id` — never silently tombstones.
  - *Filename caveat (live-found):* `exposures.filename` is the on-disk stem **minus** the `_0_001`
    frame suffix; relink/dedup must reconcile full vs frame-stripped stems (`dbkey_of_stem`).
- **Per-field geometry source columns:** each geometry field carries a `*_source` provenance value.
  The backend currently emits only **`default` | `prp` | `setup`** (`geometry.jl`:163,171,179-202;
  legacy/unset rows → `default`, `routes_experiments.jl`:141-146). `computed` and `human` (override)
  are **designed but not yet emitted** — do not switch UI on values the backend never sends. Derived
  analysis tables are untouched by regrouping.
- **Grouping flag is NEVER a column on `samples`.** It is **derived/recomputed every grouping pass**
  and returned on the loads roll-up. The only durable flag state is the `grouping_flag_dismissed`
  **event** (§11.4).

### 11.3 Sample identity & naming

- **Cross-rescan identity = the `(load_id, slot_index)` coordinate**, never a text label. This makes
  renames safe and lets a sample be re-derived from its coordinate (rebuild-not-log-derivable), not
  folded from the event log.
- **`display_name → name` collapse.** The two text columns collapse to one `name` (plain TEXT,
  whitespace-trimmed only). **Drop the old `[A-Za-z0-9._-]` charset rule** (samples are addressed by
  numeric id, no slug constraint) and **drop the `UNIQUE` index** on the name.
  - *Migration ordering (load-bearing, THREE steps):* inside `migrate_samples_name_collapse!`, guarded
    by `if "display_name" in existing` (`db.jl`:699-707), the body is **(1) DROP INDEX
    `samples_unique_name` → (2) DROP COLUMN `name` (guarded by `if "name" in existing`) → (3) RENAME
    `display_name` → `name`**. All three matter: step (1) because SQLite re-points a surviving index
    onto the renamed column (a unique index left in place would re-impose label uniqueness and reject
    legitimate duplicates like `HA85 (S01P15)`); step (2) because without dropping the pre-existing
    `name` column the RENAME throws a duplicate-column error. The whole collapse runs after
    `migrate_samples_naming!`, in **one** migration commit (no dual-accept window).
  - *Naming write sites:* `ingest.jl`:~424-443 (the cell-adoption path, `name_source = CASE WHEN
    load_id IS NULL THEN 'user' ELSE name_source END` at `:434`). Note `db.jl::create_sample!` is
    **already converted** to `name` (`name::AbstractString`, no `display_name`) — verify on migration,
    don't re-edit. Plus any test fixture passing `display_name:`.
  - *Fresh-DB caveat (B5):* the fresh `CREATE TABLE IF NOT EXISTS samples` (`db.jl`:92-93) still
    declares **both** `name TEXT` and `display_name TEXT`; the one-column end state is reached on fresh
    DBs too, via the same `migrate_samples_name_collapse!` (`if "display_name" in existing`, `db.jl`:699
    — runs on fresh DBs). Do **not** drop `display_name` from `CREATE TABLE` without also removing the
    collapse migration's `display_name` guard — **remove both together or neither**, or the RENAME never
    fires and `display_name` is left dangling.
  - *Name provenance is internal.* `name_source` (`auto`/`user`) powers never-clobber; it is **never a
    Configuration badge** (only *geometry* fields show a visible source tag).
- **Default auto-name = `<label> (SNNPMM)`** — a filename-derived label + a season/episode coordinate
  parenthetical (`S` = set/load index, `P` = slot position, zero-padded): `HA85 (S01P15)`,
  `JC C04 (S02P05)`.
- **Cross-load reshoot identity comes ONLY from the filename label** / notebook / human — **never from
  slot position** (slot numbers are not stable across loads; a counter reset to 1 signals a fresh
  rack). Default: each load's slots are new samples; sameness is surfaced as a **Merge prompt**. The
  human retains sole authority on merge/split; the UI never commits silently.

### 11.4 Structural edits, flags, undo

- **Grouping flag is backend-produced.** The merge/split flag (and its reason) is computed by the
  grouping core and returned on the loads roll-up. The frontend reads `sample.flag`, **never
  re-derives** it (it reflects the grouper's own clustering confidence, e.g. a single-frame fallback,
  which the frontend can't reconstruct).
- **Dismiss ("Keep separate") is a durable event** (`grouping_flag_dismissed`) — it survives rescans
  (reconciled into the roll-up), not session-local.
- **Flagged count — target: single source; shipped: divergent (B6).** The *flag* itself (`sample.flag`)
  is single-sourced by the backend rollup (`derive_sample_flags` / `get_loads_rollup`). But the
  frontend currently re-derives the *count* per surface with **inconsistent predicates** —
  `LoadFold.tsx` and `GroupingReviewPage.tsx` filter on truthiness (`s.flag`), while
  `ExperimentCorpusPage.tsx` filters on `s.flag !== null` — so the counts CAN diverge. **Target (do it):
  add one exported `flaggedCount(loads)` helper and normalize all three call sites to one predicate**,
  so the banner / gallery card / grouping summary cannot disagree.
- **Merge dedup (before the repoint UPDATE):**
  - `series_samples`: `DELETE` the loser's membership where the survivor already belongs to that
    series, then repoint; a warning toast names the affected series (an error would block merging any
    series member). **Then COMPACT** (ratified, Jonathan 2026-06-22): renumber the surviving rows to
    `0..n-1` ordered by `position` per affected `series_id`, inside the same merge `with_idempotency`
    block. Shipped leaves gaps today (`routes_grouping.jl`:150-205) — that is the gap to close.
    *Rationale: `position` is NOT pure sort order — the frontend always renumbers densely on save, and
    the backend add-to-series path appends via `:position => length(samples)` (`db.jl`:1260) which
    assumes density (a gap makes it recompute a colliding position → `UNIQUE(series_id, position)`
    violation). Compacting keeps the DB consistent with the dense-position invariant the rest of the
    system already enforces.* Add a test asserting positions stay `0..n-1` after a merge.
  - `sample_tags` (`UNIQUE(sample_id, key)`): drop the loser's tag where the survivor holds that key
    (survivor wins).
  A saved Series is preserved on merge by **re-pointing membership, never by deleting a sample**.
- **Undo is session-local** (frontend `useUndoStack`, cleared on reload). `user_actions` rows are an
  **audit trail, not a server-side multi-row undo API** — a single `undoes_event_id` can't reverse a
  multi-row reassignment.
- **Never-clobber.** User edits (sample name, status, curation, geometry override) survive rescans;
  rescan only refreshes derived facts. The scan UPDATE is **gated on `name_source='user'`**.

### 11.5 Migration of pre-rework production DBs

- Migrated experiments collapse `display_name` (e.g. `2-2 + LL37 1:1`, `HEPES Only`) into `name` with
  `name_source='user'` — so a migrated DB reads as if freshly curated. A fresh scan would auto-derive
  `<label> (SNNPMM)` names; the migration **keeps the science** (the human labels).
- **Op-queue `SCHEMA_VERSION` is already `5`** (frontend `persistence.ts`:16, bumped 4→5 in Phase E1
  Task 1b — *done*, not pending; there is no backend integer schema version, the backend uses
  sentinel-gated `MIGRATION_*` in `db.jl` + SQLite `PRAGMA schema_version`). On a version bump the
  persistence framework **drops** stale frontend ops with the existing "N edits couldn't be restored"
  toast; it does **not** rewrite payloads. A pre-deploy `update_sample` op carrying `display_name` is
  simply dropped (the user redoes the rename).
- **Deploy sequencing (cross-version SSE hazard):** `applyRemoteToCache`'s `update_sample` branch does
  a blind untyped spread `{ ...old, ...(payload ?? {}) }` (`applyRemoteToCache.ts`:333) — it does NOT
  flip any key. So a pre-deploy SSE frame still carrying `display_name` would splice a stale field onto
  a `Sample` now keyed off `name`. Mitigation is **deploy sequencing**, not a key flip. (The `name`
  read lives in the separate `sample_renamed` branch, `:374`.)

---

## 12. Backend internals (grouping, geometry, SSE, scheduler)

Expands §9 with the load-bearing internals the surfaces depend on.

### 12.1 Grouping algorithm

- **Frame count is variable — never assume 4.** Infer each acquisition by clustering *consecutive
  same-stage-position* frames, however many (real data has 3- and 6-frame acquisitions).
- **Load segmentation is gap-RELATIVE, not a fixed cutoff.** A gap > `gap_k` (=10.0) × `median(consecutive
  gaps)` starts a new load (`_segment_loads_with_flag`, `grouping.jl`:111-160; `:99` is the
  `_segment_loads` wrapper); if no gap exceeds the threshold it is `:unimodal_fallback` (whole dir =
  one load). The observed 53 s-within / 533 s-between bimodality is the *motivation* for gap-relative
  detection, **not** a hard 300 s constant (that value does not exist in code). The only absolute time
  threshold anywhere is `_assign_sessions`'s 3.0 h macro-session boundary (`grouping.jl`:522).
- **Within-load split suggestion:** when stage position jumps mid-group, the grouping core emits the
  split flag with `jump_from` / `jump_to`. The `position jumps X → Y mm` divider itself is a **frontend
  render** (`SampleFold.tsx`:149) of that backend flag, not emitted by `grouping.jl`.

### 12.2 Integration-trace naming (the truncated-stem fallback) — load-bearing

Integration `.dat` traces are named by the **truncated** stem (frame suffix dropped):
`JC_C01_1_S2453_tot.dat`, not `JC_C01_1_S2453_0_001_tot.dat`. The frame suffix that is *correct* for
the per-frame TIF model is *wrong* for the per-acquisition trace name. The resolution is **one shared
helper, `resolve_trace_path` (`pipeline.jl`:895-906)** — it tries the exact name first, then strips
the frame suffix (regex `r"_\d+_\d+$"`) and retries. The three call sites all route through it (do NOT
re-implement): `analyze_exposure!` (`pipeline.jl`:992), `routes_trace.jl`:27, and `series.jl`:422.
(This is the runtime counterpart to the resolver's `{name}_tot.dat` pattern in §4.3.)

### 12.3 Geometry derivation

- **Flight path:** source the setup file's **AgBe-calibrated `Mean distance`** (source `setup`) when
  present, falling back to PRP `Pipe length` (source `prp`). The two differ ~6.4 % (1809.5 vs
  1700 mm) and propagate into *every* q value; the calibrated number is correct. **This gap is EXPECTED
  and intentionally NOT flagged** — AgBe calibration refines the nominal value (option A, ratified
  2026-06-19, `geometry.jl`:211-217); the calibrated value is used silently.
- **Beam center** comes from the **analysis-dir setup files** (`source='setup'`, `geometry.jl`:178-196),
  derived from the agbe calibration. **Pixel pitch** comes from the **PRP detector-model lookup**
  (`detector_pixel_size_um`, `source='prp'`, `geometry.jl`:160-176) — NOT the setup files.
- **Energy / flight-path / detector** parse from the PRP. All geometry is **experiment-level and
  constant** across samples (it won't drift mid-run).
- **Discrepancy flags fire ONLY for intra-source variance** — PRP `Pipe length` varying across PRPs,
  multiple detector models, or energy varying — emitted by the backend as
  `manifest.discrepancies = [{field, message}]` (`routes_fs.jl`:309). The setup-vs-PRP flight-path gap
  is **not** such a discrepancy (previous point). *The Configuration **"geometry check found N issues"
  banner** that renders these is **unbuilt** (a frontend §10 target); the string is illustrative, not a
  literal — the backend already produces the data.*
- **Per-exposure, store only `timestamp` + `exposure_time`** from the PRP. Motors / slits /
  attenuator / intensity are noise — not stored, not surfaced.

### 12.4 Progress SSE (precise contract)

`broadcast_progress!` calls **`_try_put!` directly** (not `broadcast_event!`): no `user_actions` row,
no FK contract, no `event_id`. It rides the existing **`curation`** SSE event name; the frame carries
a positive non-zero `experiment_id`. The frontend `ingest_*` cases discriminate on **`kind`** and
read **`payload.experiment_id`**, never `remote.entity_id`; the four `ingest_*` switch arms land
**before** the `default:` arm in `applyRemoteToCache.ts`. The channel cap is 64 → a ~680-exposure
scan can drop intermediate frames → **treat `ingest_complete` as the authoritative terminal state.**
`phase="rescan"` rides only the **started + progress** frames; the terminal `ingest_complete` /
`ingest_failed` omit it (`server.jl`:346-356), and the frontend clears `ingestInFlight` on any
terminal frame regardless of phase.

### 12.5 Auto-rescan scheduler (tiered backoff)

A per-experiment poll `Timer` (body `server.jl`:262-291, `_rescan_tick!` `:317`) starts when the first
scan completes and runs the cheap change-check first. The scheduler `start_rescan_scheduler!`
(`server.jl`:261) is armed from the first-scan-complete path (`routes_experiments.jl`:236, right after
the `ingest_complete` broadcast), re-armed on a detected change (`server.jl`:368/389), and on a forced
pattern-edit rescan (`routes_experiments.jl`:469). **Tiered backoff (concrete defaults,
`server.jl`:317-322):** `fast_interval=3600 s` (hourly)
while exposures still arrive → after `ticks_before_daily=6` empty ticks → `daily_interval=86400 s` →
after `ticks_before_stop=3` empty daily ticks → **stopped** (manual-only). A manual rescan or new
files re-arms the fast tier; bounded ~1 week.
- **Persist tier state** in the DB (`experiments.last_scan_tier` + `consecutive_empty_ticks`) so a
  server restart doesn't reset every quiet experiment to fast tier.
- **Reentrancy:** guard with a **per-experiment `ReentrantLock`** (NOT the global `_DB_WRITE_LOCK`);
  run the heavy scan/analyze **off the timer callback** (`@spawn`) so it doesn't stall the libuv timer
  loop. `_rescan_tick!` emits the `ingest_*` lifecycle tagged `phase="rescan"` on the started+progress
  frames (the terminal frame omits phase — §12.4).

### 12.6 Own-op cache reconciliation (structural mutators)

**Already built — see §9** ("Own-op cache reconciliation is wired"). Each mutator's `onSuccess` in
`lib/queue/mutators/grouping.ts` calls `invalidateLoads` (invalidating loads + samples +
corpusSamples), because the own-op replay path runs `applyPostStateOnly` only and never calls
`applyRemoteToCache`. Kept here as a pointer so the internals index is complete.

---

## 13. Component inventory

**Already built & committed** (an earlier draft wrongly listed these as net-new; they exist):
- **`AcquisitionTimeline`** — a distinct chart type (a `TracePlot` cannot render it): an
  exposure-count-per-slot **grouped bar chart, one cluster per session**, labeled with dates + frame
  counts, on d3/SVG `PlotFrame`/`Axis` scaffolding. `src/print/components/AcquisitionTimeline.tsx` +
  `plot/AcquisitionChart.tsx` (`b5378f34`).
- **`useUndoStack`** — `src/hooks/useUndoStack.ts` (`47a328b5`), already consumed by
  `SeriesScopingPage.tsx`:2,314 (`HistoryEntry` at `:79`). A shared hook for grouping-review, geometry
  overrides, and scoping.
- **`Dropdown` / `Menu`** — `src/print/ui/{Dropdown,Menu}.tsx` (`8e35b478`); the reusable dropdown for
  the per-exposure **Move…** picker and other list pickers.
- **`Toast` action slot** — `Toast.tsx` has `action?: ToastAction` `{label,onClick}`, wired in
  `ToastContainer` (the function exported from `Toast.tsx:49` — there is no separate
  `ToastContainer.tsx` file; lines `:9,67-71,121-127`).
- **`PageFrame` width keys** — `home` (`max-w-[1080px]`, `PageFrame.tsx`:12) / `experiment`
  (`max-w-[1280px]`, `:13`) exist.
- **`GroupingBulkBar`** — the grouping-review bulk bar is **already its own component**
  (`print/components/GroupingBulkBar.tsx`, props `{count, noun, primaryLabel, primaryEnabled,
  onPrimary, onClear}`, used at `GroupingReviewPage.tsx`:471-478) — it is NOT a CullBar reuse.

**Genuinely-outstanding component items (two):**
- **`CullBar` parameterization** — `CullBar` today exposes fixed `onReject`/`onKeep`/`onRestore`/
  `onClear` with hardcoded Drop/Keep/Restore/Clear labels (`CullBar.tsx`:42-45,92-101). Add an optional
  `actions: {label, onClick, variant}[]` for **general reuse** (NOT to replace `GroupingBulkBar`, which
  already exists).
- **`KbKey` frost variant** — `KbKey` (`src/print/ui/KbKey.tsx`) hardcodes
  `bg-plate text-ink-soft border-hair-strong` with only `{children, className}` props, so it **cannot**
  render the §7 frosted `↵` chip on a filled colored button. Add a `frost`/`variant` mode (translucent
  `currentColor`, e.g. `bg-current/15 text-current border-current/45`) so the chip inherits the button's
  text color (the Focus accent fill). This is the primitive the §7/§10 frosted-`↵` target depends on.

---

## 14. Refinements absorbed from the mining sweep (2026-06-22)

Smaller surface refinements confirmed across the recent sessions, folded into the relevant surfaces
above; collected here so none is lost:

**Gallery:** drop the phase-distribution bar from ready cards (noise) · selecting an experiment
persists to Zustand so return visits route straight in.

**Picker:** the directory input placeholder must read as a prompt, not pre-filled · the card is
full-width with the two checks on the right · the assistance hint line is *"Start typing and we
suggest matches. Tab completes, ↑↓ choose, ↵ confirms."* · the duplicate-dir check is computed locally
and shown as the right-aligned ✗ pre-flight (see §5.2), not a line below the field.

**Configuration:** the analysis-directory Edit is an **autocomplete-assisted** root-relative-suffix
field (target; suggests subdirs under the root, §5.3 — shipped is a plain `Input`, the gap) and must
not look blank when empty · the Approve action is the `Button` primitive
(`variant='accent'`), Cancel is `ghost` (no naked `<button>`) · provenance chips are **neutral gray**,
never phase-hued (a color is a label; accent is reserved for human override) · geometry overrides get
a card-level "↺ Undo last change" **and** a per-field "↺ revert to PRP/setup value."

**Grouping:** a **fuzzy/glob filter** (`HA8*`, `JC C0?`, substring) with **checkboxes that persist
across filtering** (so a cross-load merge is check → filter → check → Merge, no scrolling at
170-sample scale) · the three-level fold shows **every filename** (mono) + H-position + time, **no
`+N` truncation** (showing every filename is what makes a split *verifiable*) · a per-exposure
**Move…** affordance · **no phase/analysis info on this surface** (its only job is "did everything
load and group correctly"; loading reads as "frames loading," not "analyzing") · optional per-load
H×V position sparkline (marked "very extra").

**Corpus:** a manual **Rescan** button (`triggerScan(id, …, force=false)` = additive), disabled while
scanning, no shortcut · a per-row **`Restore`** `<button>` on dropped rows (gives `Backspace` a
visible target + a11y) · rescans render the **inline `ProgressBar`** (analyzing), not the
GroupingReviewPage takeover.

**Configuration → scan plumbing:** editing a pattern and saving must **chain a forced rescan**
(`POST /{id}/scan?force=true`) — `cheap_change_check` is file-count-only and ignores pattern changes,
so without `force` a pattern edit does nothing · `POST /api/experiments` accepts nested
`patterns: { image?, metadata?, integration? }` in the body (see §9 — body keys are unsuffixed, the
`_pattern` suffix is the DB column) so edits flow into phase 2 without a separate PATCH.

**Conformance defects found live (must not regress):** `fetchManifest` must **single-encode** its
query — build it with **one** `new URLSearchParams({path, …})` (`api.ts`:304-312, which is correct
today); do NOT *additionally* `encodeURIComponent` the values before adding them (the historical bug
was double-wrapping → double-encode → 400 → whole funnel broken, NOT `URLSearchParams` itself) · the
data-dir field must show the **root-relative suffix**, never default to the project directory (the
shipped impl regressed this and confused the user).

---

## 15. Build-time decisions (all RATIFIED, Jonathan 2026-06-22)

Small implementation choices, resolved so they don't block. All confirmed 2026-06-22 (surfaced by the
spec-review sweep).

| # | Decision | Ruling |
|---|---|---|
| D1 | `data_dir === root` display (§4.3/§5.3) — `stripRoot` returns `"."` | Render **`(root)`** (friendlier than a bare `.`). |
| D2 | Corpus dock `Restore` key-chip (§7) | **Add the `Backspace` `KbKey` chip** for consistency with Drop `X` / Keep `K`. |
| D3 | `createExperiment` return type (§9) | **Retype to `CreateExperimentResponse { id; status; name; data_dir }`** (the real 202 partial); update `api.ts`. |
| D4 | Zero of an essential type (§6) — image, metadata, OR integration `=== 0` everywhere | **RATIFIED:** all three legs essential ("for now"). Headline **`No <type> matched`** + **disable Approve** + hide `/0` denominators. No advisory tier. |
| D5 | Single `flaggedCount` (§11.4) | **Build one exported `flaggedCount(loads)` helper** + normalize the predicate across the 3 call sites. |
| D6 | `series_samples.position` gaps after merge (§11.4) | **RATIFIED: compact** — renumber survivors to `0..n-1` per affected `series_id` in the merge block + test. (`position` is not pure sort order: the `length(samples)` append idiom, `db.jl:1260`, assumes density.) |
| D7 | Indexing readout `M seen` (§5.3 vs §6) | **Dropped** — `N exposures found`, no denominator (already folded into §5.3/§6). |

The canonical key for an experiment is its **`data_dir`** (the backend 409/dup compares `data_dir`,
`routes_experiments.jl`:181-188; the picker's check of both `root` and `root/data` is a heuristic for
it). The Approve 409 error path (§5.3) reads the body `{error, path, experiment_id}` and toasts
"This directory already backs an experiment" without navigating.
