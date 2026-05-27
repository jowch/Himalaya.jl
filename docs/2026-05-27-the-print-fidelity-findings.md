# "The Print" — consolidated fidelity findings

**Date:** 2026-05-27 · validates `main` @ `c5591c1`
**Companions:** the code-level gap report (`docs/2026-05-27-redesign-gap-report.md`), the remediation plan (`docs/superpowers/plans/2026-05-27-the-print-remediation-plan.md`), the canonical design reference (`DESIGN.md`, now "The Print"), and the five mockups (`docs/redesign-mockups/*.html`).
**Method:** a live impeccable visual pass (app run against the dev-copy DB, Playwright over all five surfaces + focus workspace, compared to the mockups) plus three parallel shipped-vs-mockup fidelity audits, one per workflow stage. This is the **evidence base** the remediation issues (#221–#233) reference by finding ID.

**How to read this.** Finding IDs are **scoped to their surface section** (each audit numbered independently, so e.g. a sample-table `M-1` and the focus `L-6` are unrelated to a builder finding of a similar letter). Always resolve an ID *within the section named by its owning issue*. The matrix below maps every finding to its owning issue; start there.

**Severity:** P1 workflow-breaking / promised-and-missing · P2 degraded but functional · P3 polish.
**Status:** net-new (unflagged before this pass) · known-deferral (cross-ref #207/#208/#209/#162) · intentional (deliberate simplification, noted).

---

## 0. Cross-cutting themes (read first)

1. **The dark↔Print control seam (dominant).** The surfaces were migrated to "The Print" (paper, terracotta kickers, the contact-sheet/trace-plate metaphors) but are built largely from **legacy dark-"Darkroom" components rendered unchanged on a light shell**. Those components consume the neutral ramp (`bg-bg`, `bg-bg-elevated`, `text-fg`, `text-fg-muted`, `border-border`) and the ice-blue `--color-accent` (hue 220). Because the app default is still `theme:"dark"` (`state.ts:284`), they resolve dark on warm paper. **Root-cause fix = R0a (#221)**, which resolves the majority of theme findings *by inheritance* (redefine the token values, retire the dark theme). Theme findings below are mostly "verify after R0a."
2. **Newsreader serif is not loaded at all** (`styles.css:1-4` loads only Plus Jakarta Sans). Every mockup display heading uses Newsreader; the serif title is a defining "The Print" trait. App-wide miss, cheapest high-identity fix. → R0a (load) + per-surface application.
3. **The phase palette is dark-tuned** (`phases.ts`, L≈0.78–0.80, low chroma) not paper-tuned (mockup L≈0.50–0.58, AA on `--plate`). → R0b (#222).

---

## 1. Finding → issue matrix

| Surface section | Finding IDs | Owner issue |
|---|---|---|
| Token system (cross-cutting theme) | ST T-1..T-8, FW A-1/A-3/A-4, SB B-A/B-B/B-C, SC S-B, toggles L-5/L-7 | **R0a #221** |
| Phase palette | FW A-2 | **R0b #222** |
| Residual semantic tokens | ST T-4/T-7/T-8, SB B-C | **R0c #223** |
| Sample table — contact sheet | ST M-1/M-2/M-3/M-6/M-7/M-10, L-3/L-4/L-5/L-8/L-11 | **R1 #224** |
| Sample table — loupe | ST M-8/M-9, L-9/L-10 | **R2 #225** |
| Focus — header | FW L-6/L-7/L-8 | **R3 #226** |
| Focus — rail | FW L-9/L-10/L-11/L-12 | **R4 #227** |
| Focus — affordances | FW F-11/F-12/F-13/F-14 | **R5 #228** |
| Series — folio | SF F-A/F-B/F-C/F-D/F-G/F-H, X-1 | **R6 #229** |
| Series — scoping | SC S-A/S-C/S-D/S-E/S-F/S-G/S-H | **R7 #230** |
| Series — builder | SB B-F/B-G/B-H/B-I/B-J | **R8 #231** |
| Docs — DESIGN.md drift | DOC-1 | **R9 #232** |
| Backend — image route | BE-1 | **R10 #233** |
| (deferred) reflections table + q-link row | FW F-8/F-10 | **#209** |
| (deferred) heatmap + peak-tracking | SB B-D/B-E | **#208** |
| (deferred) tag on-ramp | ST M-4/M-5 | **#207** |

---

## 2. Sample table (Contact sheet + Loupe) — `sample-table.html`

### Theme (→ R0a / R0c)
| ID | Sev | Location (shipped → mockup) | Shipped vs mockup | Status |
|---|---|---|---|---|
| T-1 | P1 | `LoupeFrame.tsx:35,49` → mockup `.big-frame{background:var(--frame-edge)}` | Big frame backed by dark `--color-bg` (cool, hue 250); "Dropped" label `text-bg` invisible on accent. Mockup = warm `--frame-edge` (hue 55). | net-new |
| T-2 | P1 | `LoupeFrame.tsx:35`, `LoupeSidebar.tsx:86,107,117,133,152` `border-border` → `--hair`/`--hair-strong` | All loupe borders are dark `--color-border` (oklch 0.27) → near-black hairline on paper. | net-new |
| T-3 | P1 | `LoupeSidebar.tsx:86,108,134,153` `bg-bg-subtle` → `--paper-sunk`/`--plate` | Verdict box + tag pills filled near-black on light paper. | net-new |
| T-4 | P1 | `LoupeSidebar.tsx:91` `bg-success`; `:120-123` `accent`; `LoupePage.tsx:157,170` → `--im3m`/terracotta | Kept-dot ice-green not sage; all accents ice-blue not terracotta. **Semantic remap → R0c.** | net-new |
| T-5 | P1 | `LoupeFrame.tsx:48`, `LoupeSidebar.tsx:91,121,123` `bg-accent`/`text-accent` | "Dropped" badge + "Representative" mark ice-blue, not terracotta (the surface's one brand accent). | net-new |
| T-6 | P1 | `LoupePage.tsx:150,157` | Not-found card: dark border + ice-blue link. | net-new |
| T-7 | P2 | `ThumbnailGallery.tsx:62-63,76` `ring-accent`, `bg-accent/85`, `text-bg` → `.thumb`/`.thumb.sel` | Shared strip uses ice-blue selection + dark hover rings (component shared with dark Inspect, never re-themed). **→ R0c.** | net-new |
| T-8 | P2 | `LoupeFrame.tsx:54` `text-ink-faint` caption over dark frame | Dark caption on dark frame → low contrast; mockup uses a light frame-tag. **→ R0c.** | net-new |

### Missing affordances (→ R1; tagging → #207)
| ID | Sev | Location → mockup | Shipped vs mockup | Status |
|---|---|---|---|---|
| M-1 | P1 | `SamplesPage.tsx:60-68` → `.head .progress` (sample-table.html:480-484) | No "N/M samples screened" count + terracotta progress bar. | gap report F-3 |
| M-2 | P1 | `ContactSheetRow.tsx:197-201` → `.screened-mark` (:180-188) | No per-sample hollow/filled screened dot, no unscreened row tint. Comment defers to #162. | F-3 / #162 |
| M-3 | P1 | `ContactSheetRow.tsx:226-249` → floating `.cullbar` (:570-575) | Cull/batch-reject is a per-row inline bar, not the floating bottom-center bar; selection never crosses rows (intentional per-sample query fan-out). | F-3 (partial-intentional) |
| M-4 | P1 | `ContactSheetRow.tsx:282-291` `+ tag` disabled "coming soon" | No tagging on-ramp. | **#207** / F-1 |
| M-5 | P2 | `ContactSheetRow.tsx:282-291` always-on `+ tag` vs mockup hover-reveal on tagged rows | (Moot while disabled.) | **#207** |
| M-6 | P2 | `ContactSheetRow.tsx:296-298` hardcoded "Not indexed" → `.pchip` | Static gray string for every row; no phase chips, not even the leading hollow dot. | net-new (seam intentional; styling gap) |
| M-7 | P3 | `ContactSheetRow.tsx` thumb → `.thumb .fno` (:214-218) | No zero-padded frame-number badge on thumbnails. | net-new |
| M-8 | P3 | `LoupeSidebar.tsx:66-80` → `.meta-list` (:887-895) | Meta is Filename/Kind/Frame/Status; mockup is frame/integration-time/collected/**signal-strength bars**. Signal meter absent. | net-new |
| M-9 | P3 | `LoupePage/SamplesPage` topbar → `.seg #view-seg` (:460-463) | No "Contact sheet \| Loupe" segmented switch (split into two routes + a back link). | net-new (route-based acceptable) |
| M-10 | P3 | `LoupeFrame.tsx:45-57`, `ContactSheetRow.tsx:84-88` → `.big-x`/`.thumb ✕` (:361-362) | Rejected frames dim only; no hand-skewed two-stroke grease-pencil ✕ SVG. | net-new |

### Layout / typography / copy (→ R1 / R2)
| ID | Sev | Location → mockup | Shipped vs mockup | Status |
|---|---|---|---|---|
| L-1 | P1 | no Newsreader import → mockup serif h1/h2/progress | Serif display face (the type signature) not loaded; headings fall to sans. | net-new (→ R0a) |
| L-2 | P2 | `LoupePage.tsx:175` `text-2xl` → 26px Newsreader 500 | `text-2xl` not in the 5-step scale → Tailwind default 24px sans (off-scale). | net-new |
| L-3 | P2 | `SamplesPage.tsx:60-68` → `.head h1`+`.sub` (:470-485) | Missing the beamtime `<h1>` + descriptive framing sentence (only kicker + one-line scope). | net-new |
| L-4 | P2 | `ContactSheetRow.tsx:51` `aspect-[3/4]` → square `.thumb` 62px | Portrait thumbs vs square detector frames (distort/letterbox). | gap report F-4 |
| L-5 | P2 | `ContactSheetRow.tsx:22` grid → `.COLS` (:159) | Column widths diverge; no per-row min-height 92px. | net-new |
| L-8 | P3 | `SamplesPage.tsx:117-122` legend → `.kb-legend` (:500-506) | 4 plain-text hints vs 5 with styled keycap chips (adds ⇧-click range). | net-new |
| L-9 | P3 | `LoupePage.tsx:178-181` subtitle → `lh-sub` (:858-859) | Shows experiment·name; mockup shows sample-id · "exposure N of M". | net-new |
| L-10 | P3 | `LoupeSidebar.tsx:153` `key: value` → bare value (:915) | Tag pills show `key: value`; mockup tags are free-form single tokens. | net-new |
| L-11 | P3 | `ContactSheetRow.tsx:208` strip `h-16` → 92px / 62px thumbs | Row rhythm off from the square-frame contact sheet. | net-new |
| (L-6/L-7 here are sample-table layout polish: page max-width/centered shell, contained "sheet" card surface — fold into R1.) |

---

## 3. Focus workspace — `focus-workspace.html`

**Systemic:** built from legacy dark components on a Print shell (`CorpusShell` paints `bg-paper`, inner `PlotCard`/`TraceViewer`/`PhasePanel`/`IndicesCard`/`MillerPlot` consume dark tokens). Most rows below share the R0a fix.

### Theme (→ R0a / R0b)
| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| A-1 | P1 | `styles.css:6-68` `@theme` + `state.ts:284` default `"dark"` | Dark tokens are the live default; near-white `text-fg`/`text-fg-muted` on cream paper. | net-new (→ R0a) |
| A-2 | P1 | `phases.ts:23-30` → mockup `--pn3m/--lam/--im3m/--hex` | Phase colors dark-tuned (L 0.78-0.80, low chroma); mockup paper-tuned (L 0.50-0.58, AA on plate). | net-new (→ R0b) |
| A-3 | P1 | `TraceViewer.tsx:235,408,739`; `PlotCard.tsx:627`; `DetectorRingOverlay.tsx:110` | q-link/accent = ice-blue hue 220; mockup terracotta hue 38. The headline q-link lights the wrong hue. | net-new |
| A-4 | P2 | `styles.css:92` `color-scheme:dark` + `:104-115` grain | Native controls dark-scheme; fractal grain (tuned for dark) muddies the paper. | net-new |

### Header (→ R3)
| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| L-6 | P1 | `FocusWorkspaceLayout.tsx:60` (`<PlotCard/>` prop-less) → `PlotCard.tsx:60-104,394-415`; mockup `.plate-head` (:468-473) | Title shows legacy picker "pick an experiment / click to change": route seeds `activeSampleId` but not `activeExperimentId`, so `experimentName` is undefined → italic picker branch. Mockup wants kicker "Integration" + Newsreader serif **sample name** + mono `smp · beamtime · representative exposure`. **Fix:** focus-variant header / `headerSlot` on PlotCard; drop the `onTitleClick→openNavModal` picker affordance here. | live L-6 |
| L-7 | P2 | `PlotCard.tsx:394` `text-title` (sans 15px) → `.plate-head h1` (Newsreader 27px/500) | Sans title vs serif display hero. | net-new |
| L-8 | P2 | `PlotCard.tsx:469-475` → `.kicker` (:147-150) | No terracotta uppercase "Integration" kicker above title. | net-new |

### Rail as output (→ R4; reflections table → #209)
| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| L-9 | P3 | `PhasePanel.tsx:254-257,280,303,327` → rail-h "Phase call"/"Candidate indexings" (:523,528) | No distinct "Phase call" output block (one card/phase: serif name, score bar, lattice, series ratio); no coexistence framing. | net-new |
| L-10 | P2 | `PhasePanel.tsx:46-155,134-140` → `.cand` multi-select `.c-mark` (:916-927) + `.phasecall` (:904) | Candidates are +/− buttons, not multi-select checkboxes; no series ratio √2:√3:√4 surfaced. | net-new |
| L-11 | P3 | `PhasePanel.tsx:293-294` (`hoveredIndex`→accent) → `setPreview(name)` (:662-688) | Candidate hover-preview uses the accent channel; mockup previews in the **phase's own colour** + dims losing peaks on a swap. | net-new |
| L-12 | P2 | `IndicesCard.tsx:36-50` (Miller inset, dark tokens) → mockup rail has no Miller plot | Legacy √N·q Miller inset not in the Print rail (extra dark-token surface). | net-new |
| F-8 | P1 | absent in `FocusWorkspaceLayout.tsx:58-63` → `.lower` two-panel (:488-515) | Reflections-table panel (phase·hkl·q·d) structurally absent. | **#209** |
| F-10 | P2 | q-link wires 2 surfaces → `setHot()` peak+ring+row (:651-655) | q-link lights trace-peak↔ring only; 3rd target (reflection row) absent. | **#209** |

### Affordances (→ R5)
| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| F-11 | P1 | `FocusDetectorPanel.tsx` (comment :17 drops rep controls) → `.expo` strip (:230-244,840-851) | Representative-exposure switcher (mini-detector thumbnails, set-rep, click-to-switch) entirely missing; exposure auto-picked with no user switch. | net-new |
| F-12 | P1 | `FocusWorkspaceLayout.tsx:72-79` Notes `hidden < xl`; no topbar control → `#notes-btn` (:458) + `@media(max-width:1320px)` (:434-438) | Below `xl` Notes margin is simply hidden with **no fallback** → Notes unreachable on narrow viewports (mockup collapses it to a topbar toggle/drawer). | net-new |
| F-13 | P2 | `CorpusTopbar.tsx:128` → `.stepper` (:453-457) | No per-sample stepper ("‹ sample 4 of 9 ›"); the focus surface's primary inter-sample nav. | net-new |
| F-14 | P2 | `CorpusTopbar.tsx:15-19,97-108` | Index stage-tab permanently disabled even on `/sample/:id` → no active stage tab. | net-new (Phase-4 #181 may own) |

*(`FocusNotesMargin.tsx` is already on Print tokens — the one correct inner component. `FocusWorkspaceLayout` grid `1fr_348px_250px` matches the mockup.)*

---

## 4. Series — Folio — `series-folio.html` (→ R6)

Folio uses Print tokens throughout (theme-clean); gaps are missing affordances.

| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| F-A | P1 | `SeriesFolioCard.tsx:43-63` → `figSVG()` (:440) | Grey placeholder box (height ∝ member count), not a live per-series mini-waterfall SVG. Feasible: reuse `MultiTracePlot`/`MemberTraceLayer` geometry as a read-only variant. | N-1 (intentional placeholder, comment :14) |
| F-B | P2 | `SeriesFolioCard.tsx:49-62` → `.ps` (:243) | 3 fixed swatches + "+N more" vs one cell **per sample** (full transition strip), coexistence = 2-phase gradient, + caption ("Pn3m → Lamellar"). | net-new |
| F-C | P2 | `SeriesFolioPage.tsx` (no pill) → `.pill-new`/`.pill-draft` (:221) | No "+N new match" recipe pill; draft is only lowercase text, no pill. | net-new |
| F-D | P2 | `SeriesFolioPage.tsx:79-108` → controls row | Search + 2-way sort (Recent/Title) only; mockup has filter chips (All/Has transition/Cross-experiment) + 3-way sort (Recent/Variable/Largest). | net-new (chips deferred per comment :35) |
| F-G | P2 | `SeriesFolioPage.tsx` header → `+ New series` btn-solid (:288) | No primary "+ New series" action; only the empty-state text leads in. | net-new |
| F-H | P3 | `SeriesFolioCard.tsx` footer → `card-foot` | No footer rule / beamtime provenance / cross-experiment note / edited-timestamp. | net-new |
| X-1 | P2 | no Newsreader; `SeriesFolioCard.tsx:67` etc. | All series headings sans, not serif display. | net-new (→ R0a) |
| (X-2/X-3: display-sized serif titles + kicker→title→sub stack — fold into R6.) |

---

## 5. Series — Scoping — `series-scoping.html` (→ R7)

The biggest single gap. The `(key,value)` batch-write plumbing + build-gate exist; the proposal **surface** does not.

| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| S-A | P1 | `SeriesScopingPage.tsx` (whole page) | Machine-proposes/human-confirms layer absent: no autogroup summary, no parsed "Ordered by" field (only "Ordering variable: —"), no per-row **trace sparklines**, no amber "check the read" flags, no "Himalaya also found" loose-match section, no preview phase strip, no narrative gate footer. `proposeOrdering.ts` exists but on a cold corpus `orderingKey` is undefined → degrades to a bare manual-entry list (and stays there even once tags exist). | live L-8 |
| S-B | P1 | `SeriesScopingPage.tsx:199`, `ScopingConfirmModal.tsx:53` `bg-accent` | Both "Confirm & build" buttons ice-blue; mockup `.build-btn{background:var(--ink);color:var(--paper)}`. | net-new (→ R0a/R7) |
| S-C | P2 | `SeriesScopingPage.tsx:199` `disabled:opacity-40` | Disabled build = opacity fade of wrong-hue button vs mockup's distinct greyed ink style. | net-new |
| S-D | P2 | `SeriesScopingPage.tsx:114` | Full-width flat list vs centered 760px **worksheet plate** (white plate, hair border, shadow). | net-new |
| S-E | P2 | `SeriesScopingPage.tsx:172-182` | Permanent text `<input>` per row vs "confident ink value, re-openable on click; flagged = amber underline + row tint" (confirm-not-fill-out). | net-new |
| S-F | P2 | `SeriesScopingPage.tsx` | No reorder grips, no "low to high" indication, no "N samples · low to high" line, no in-session Undo/⌘Z. | net-new |
| S-G | P3 | `SeriesScopingPage.tsx:116-124` | Header is kicker + "Ordering variable: —"; mockup leads with editable serif series name + autogroup sentence + metadata-as-byproduct footnote (shipped only in the confirm modal). | net-new |
| S-H | P3 | `SeriesScopingPage.tsx:301` | No "Discard" affordance (relies on browser back). | net-new |

---

## 6. Series — Builder — `series-builder.html` (→ R8; heatmap/tracking → #208)

Structurally closest to its mockup, but carries the most theme debt.

| ID | Sev | Location | Shipped vs mockup | Status |
|---|---|---|---|---|
| B-A | P2 | `GroupingModeToggle.tsx:60-62`, `AnnotationToggles.tsx:51-53` `bg-bg-subtle text-fg` | Active state dark-theme tokens; sibling `RepresentationToggle.tsx:29` correctly `bg-ink text-paper`. Mockup `.btn.on{background:var(--ink);color:var(--paper)}`. | live L-5 (→ R0a) |
| B-B | P1 | `MemberMetaRow.tsx` (gutter, :287,312,330,350-369,375,403-412,450) | Entire read-mode metadata gutter on dark tokens (`text-fg*`, `hover:bg-bg-elevated`, `ring-accent`, edit inputs `bg-bg border-border`). Larger than the two toggles. | net-new (→ R0a) |
| B-C | P2 | `MultiTracePlot.tsx:593-594,643-645`; `MemberMetaGutter.tsx:301,355` | Shared render-core tooltip/selection/insertion-line carry dark tokens + ice-blue accent. | net-new (→ R0c) |
| B-D | P1 | `RepresentationToggle.tsx:33-41` | Heatmap button permanently disabled "coming soon (#208)". | **#208** / gap F-6 |
| B-E | P1 | `SeriesBuilderPage.tsx:261-264` (no track toggle) | No "Track reflections" cross-trace peak-tracking layer; only Peak ticks/labels. | **#208** / gap F-7 |
| B-F | P2 | `SeriesBuilderPage`/`SeriesBuilderRail` | No **trace offset slider** and no **log/linear scale toggle** — core "compose the figure" controls. | net-new |
| B-G | P2 | `SeriesBuilderRail.tsx:32-43` | Full-bleed rail-collapse ships, but the floating **offset dock** (keeps offset reachable while reading) is absent (tied to B-F). | net-new |
| B-H | P2 | `SeriesBuilderPage.tsx:97-122` | No `fig-tags` kicker row, no `fig-sub`, no auto figure caption block. | net-new |
| B-I | P2 | `SeriesBuilderRail.tsx` | No autogroup card in the read rail ("Auto-grouped — Himalaya read N samples…" + Confirm/Adjust). | net-new |
| B-J | P2 | `SeriesBuilderPage.tsx:243-289` | No figure-as-**plate** container (white plate, hair border, shadow, max-width 1180→1336px full-bleed). The figure-as-printed-plate metaphor is the surface's whole point. | net-new |
| (B-K/B-L: edit-rail inputs `bg-paper` vs `bg-plate`; rail recessed-margin treatment — P3 polish, fold into R8.) |

*(`SeriesRecipeEditor.tsx` edit-mode rail is largely on-palette; Commit button correctly `bg-print-accent text-paper`.)*

---

## 7. Docs & backend

| ID | Sev | Location | Finding | Owner |
|---|---|---|---|---|
| DOC-1 | P2 | `DESIGN.md` | Documents the retired dark "Darkroom" identity; impeccable loads the wrong design system. (Addressed up front by authoring `DESIGN.md` to The Print; R9 verifies shipped matches.) | **R9 #232** |
| BE-1 | P2 | image route vs `routes_trace.jl:28` | `GET /api/exposures/:id/image` throws unhandled `ArgumentError`→HTTP 500 on missing source file; should 404/placeholder like the trace route. | **R10 #233** |
