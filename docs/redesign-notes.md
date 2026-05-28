# HimalayaUI redesign: working notes

**Status:** in progress, started 2026-05-17. Living document, updated as the page-by-page discussion converges.
**Companion files:** `PRODUCT.md` (strategic context), `DESIGN.md` (current visual baseline).

This captures a redesign exploration of HimalayaUI's UX. It is a record of a workflow-driven discussion, not a final spec. Sections marked *open* are not yet settled.

## 1. Why this redesign

The surface prompt: the current UI feels "simultaneously too much and too little." A critique grounded in the production action log (1,421 logged actions from a working copy of the DB) confirmed the complaint.

**Root cause: equal space for unequal work.** The Index workspace gives three co-equal cards (chat | trace | indices) to three wildly unequal activities. When everything is weighted the same, everything competes (too much) and nothing leads (too little).

Evidence:

- **Analyze plus peak curation is 86% of all activity** (analyze 64%, peak edits 22%). Index curation (2%), speculative building (1%), chat (1.5%), comparison (under 1%) are rare by comparison.
- **The ChatCard occupies about 22% of the screen for 1.5% of activity** (20 messages ever, across 11 of 139 samples).
- **620 of 678 exposures (91%) have zero indices.** The IndicesCard's normal state is three empty sections of hint text.

**The deeper motivator.** The real reason for a dramatic redesign is the **Compare page**. It is where the actual science happens, and even after a recent rework it does not hit the mark, because the problem is conceptual, not cosmetic. See section 2.

## 2. How the work actually flows

Domain knowledge from the user, not derivable from code or the action log.

**Unit hierarchy:**

> Experiment -> Sample -> multiple Exposures (the beam passes through different parts of the sample, so signal strength varies) -> each exposure has a 1D integration -> one representative exposure's integration is indexed -> indices -> phase call.

A sample's pattern may contain peaks from **more than one phase** — coexistence — so the phase call is a *set* of phases, not a single answer, and the indexing surface partitions peaks among them.

**The workflow has three stages, and the third is the point:**

1. **Triage.** New data arrives. The scientist reviews raw detector *images* one sample at a time, screening for flares and artifacts. Artifacts are fast to spot by eye. Then one representative exposure per sample is chosen for indexing (accepted exposures are usually redundant).
2. **Index.** The representative integration is indexed: curate peaks, read the phase call.
3. **Cross-sample analysis.** Indexed samples are compared against each other: how the diffraction evolves with concentration, time, dose, control vs. treated. **This is where the real analysis happens.** Stages 1 and 2 are preparation for it.

**The cross-check loop.** While reading an integration (stage 2), the scientist flips back to the raw image: measured ~90% to *verify* (is this trace feature real signal or an artifact?), ~10% to *re-choose* a better exposure. So image and trace are two views of one measurement, not two pages.

**The payoff artifact.** The goal of stage 3 is a standard figure: a **stacked-offset waterfall** of integration traces, stacked vertically with an offset, ordered along a meaningful variable (concentration, time, dose). The redesign exists to make building and reading that figure easy.

## 3. The converged model: three surfaces

The app should be three surfaces, one per workflow stage, each weighted to its own job. This resolves the root cause: nothing is co-equal that shouldn't be, and "empty" becomes signal rather than wasted space.

### The sample table (triage + metadata)

A tabular surface over the whole corpus of samples, spreadsheet-like — experiment is a filter facet, not a container (see section 5, "architecture revisited"). Two jobs on one surface:

**Exposure triage (culling).** The first pass is photo culling:
- **Grid and loupe co-equal** (the user works ~50/50 from each); flipping must be instant.
- **Reject-only.** No "accept." Everything is kept by default; the one action is reject. Silence is the default answer, so attention goes only to bad exposures.
- **Multi-select plus batch reject.** Select several, reject in one action.
- **Keyboard-first**, for speed.
- **No auto-flagging.** Artifacts are fast to spot by eye; the tool should not pre-judge. The surface's job is to make images big, sharp, and fast to flip.

**Sample metadata editing.** The table also carries editable metadata columns (see section 3, "Sample metadata"). Triage and "sample management" are the same surface; the metadata-entry "form" the user asked for becomes columns in a table they already use.

### The focus workspace (indexing one sample)

One sample's image, trace, and indices coexist on one surface, never navigated away from:

- The **trace** (1D integration) is the center and the main object. Keep it large.
- The sample's **detector image** is co-resident, not behind a tab, because the verifying cross-check is the 90% reason for leaving the trace.
- **Image and trace are linked by q.** A trace peak is a ring in the image at the same q. Hover or select a peak and its ring lights up: a real reflection shows a crisp ring, an artifact a blob or streak. Turns a page round-trip into a glance; fits the "show the reasoning" principle.
- A **representative-exposure switcher** handles the 10% re-choose case in place.
- **Indices / the phase call** are the output.
- **Notes / chat** are reachable, not resident: one reach away, never paying standing screen rent.

**Validated, do not redesign:** the interactive trace editing (peak add / remove / exclude) and the active-set / candidate index semantics work well and carry over unchanged. The focus workspace redesign is purely layout and visual hierarchy: reduce clutter, fix "the eye fighting the screen." The peak-finder is high precision but lower recall on weak peaks; surfacing sub-threshold candidates as "ghost" marks was considered and rejected as mostly noise.

Of the three surfaces, the focus workspace works best today. It still benefits from an upgrade, but the sample table and the series builder are higher-priority.

### The series builder (the reconceived Compare)

**The destination.** This is where stage-3 analysis happens. A "series" is a first-class object — a machine-maintained *recipe* (which samples, which ordering variable, the order) plus a frozen *plate* (the exported snapshot); see section 5, "architecture revisited". The Series stage is itself three views: a **folio** (browse every saved series, corpus-wide), a **scoping step** (create one — where the ordering metadata is collected), and the **builder** (edit and read one, as the stacked-offset waterfall below).

- Select a set of samples; the tool **proposes the grouping and order** from name patterns plus stored metadata (machine proposes, human confirms, the confirmation becomes data; the model is Narrative's scene autogrouping).
- Render as a **stacked-offset waterfall**, ordered along the chosen variable. Offset is tunable. This replaces the current overlay model, which is messy past a few traces.
- The per-sample indexing surfaces here as annotation, so the scientist reads how the phase evolves across the series.
- The output includes a publication-quality figure (PNG / SVG export).

The current Compare page's conceptualization is set aside; the rethink is "a series is the unit of analysis," not "Compare is an ad hoc overlay action."

### What dissolves: the Inspect page

"Inspect" as a standalone destination is absorbed. Its cull function joins the sample table; its image-viewing becomes the q-linked cross-check in the focus workspace.

### Sample metadata: the prerequisite, captured progressively

The smart series builder needs to know the ordering variable (concentration, time, dose). Today that is unstructured: it lives in free-text `display_name` values ("HEPES Only", "1-1 + LL37 1:1"); the `sample_tags` key/value table exists but is empty. Parsing names alone is too unreliable to trust.

The mechanism, three principles:

1. **Progressive, not up-front.** Capture metadata when it is needed (at series-building time), not at ingestion. Ingestion stays frictionless.
2. **Byproduct, not chore.** The act of assembling and ordering a series is itself the metadata. The tool guesses groupings from name patterns, the user corrects, the correction is stored.
3. **Two moments, two scopes.** A Himalaya "experiment" is one beamtime / beam configuration — it holds many unrelated samples testing a variety of things, so experiment-wide metadata columns do not generally cohere. The coherent unit for an ordering variable is the *series*. So metadata is collected at two scopes: light, free-form per-sample tags on the sample table ("what is this?", small scope, an invitation to start), and the structured ordering metadata collected at **series-scoping time** — an intermediate step (a scoping view or modal) between the sample table and the comparison screen, where a coherent sub-group has been selected and the relevant variable is finally well-defined. Reasoning about metadata for a deliberately-scoped collection is easier than for a table that shows every beamtime sample equally.

The app earns the metadata by making the payoff (the series analysis) the incentive: users have no reason to follow a standard for its own sake, but they will if it unlocks the analysis. The app defines the standard.

The lab-notebook manifest stays a supported bulk-import path for users who keep structured notes, but it is not required. The manifest itself warrants its own redesign later, deferred until the app's metadata needs are known, so the app drives the schema rather than inheriting it.

## 4. Reference interfaces

- **Adobe Lightroom (Library module):** grid plus loupe, reject-centric flagging (`X` to reject, keep by default). The "no accept, only reject" model already shipping at scale.
- **Photo Mechanic:** the professional culling-speed benchmark; keyboard-driven, grid as home.
- **Narrative Select / Aftershoot:** modern, well-crafted cullers. Also the reference for *autogrouping*: machine proposes groupings, human confirms.
- **cryo-EM micrograph curation (cryoSPARC "Curate Exposures", RELION):** the domain cousin. Same task at scientific scale. Borrow the interaction model, not the craft (these UIs are the "legacy scientific software" anti-reference).

## 5. Decisions, open questions, next

**Settled:**
- Three surfaces: a sample table (triage + metadata), a per-sample focus workspace, a series builder.
- The Inspect page is absorbed.
- The series builder is the analytical destination; the waterfall plot is the payoff artifact.
- Culling: reject-only, grid+loupe co-equal, multi-select batch reject, keyboard-first, no auto-flag.
- Focus workspace interactions are validated; its redesign is layout and visual hierarchy only.
- Metadata is captured progressively in-app (>80% of use), as a byproduct of analysis, never an up-front form.
- Priority of pain: the sample table and series builder first, the focus workspace after.

**"The Print" — the series builder:**
- Direction **"The Print"**: the series builder departs from the dark Darkroom to a light, paper, figure-as-a-plate treatment, on the logic that this surface produces a publication figure and you should compose it on the medium it will be published on. Phase colour carries the surface (the dark app rations one accent; here colour is the science). The user finds this direction beautiful and an improvement; since HimalayaUI has no entrenched identity, "The Print" is a candidate to become the app's identity, and the dark "evolve" comparison was dropped as a hedge against a risk that does not exist.
- The series builder was explored across three files — `series-builder-bold.html` (the first full mockup), `series-builder-variations.html` (four representation variants), and `series-builder-directions.html` (three alternative identities). All three were consolidated into one canonical `series-builder.html` on 2026-05-17 and retired (see "Done — series builder consolidated" below).
- `series-builder-variations.html` carried four variants in one switcher: **Full-bleed** (chrome steps away, figure owns the screen), **Heatmap** (the series as a q–ratio intensity map), **Tracked** (peak-tracking lines + d-spacing callouts), **Drenched** (phase-colour bands, the transition as a colour field). Skinnier plot throughout. The variants composed rather than competed (full-bleed is a mode, drenched a colour treatment, tracked an annotation layer, heatmap an alternate representation) — which is why they could be folded back into one trunk.
- All built on a real cubic-to-lamellar story from the DB (an LL37 titration; the lowest and 1:1 ratios measured, intermediate doses illustrative).

**Variant verdicts (user review of `series-builder-variations.html`):**
- **Heatmap** — liked, "never seen it done this way," potentially useful for highlighting weaker peaks. Keep; test on real data.
- **Full-bleed** — liked, but as a collapsible/toggleable rail, not a separate variant. Confirmed: it is a mode, not a design.
- **Tracked** — possibly useful; works best when peaks genuinely align across the series.

**Alternative directions (`series-builder-directions.html`):**
Since the app has no entrenched identity, three alternatives to "The Print" were built in one switcher, holding the layout skeleton constant so the identity is the only variable:
- **The Console** — dark warm-charcoal; the data is the only light source (traces glow). A beamline instrument readout. Departs from The Print on the light/dark axis.
- **The Drafting Table** — pale cool ground, a full measured grid, monospace, an engineering title block. A crystallographer's drawing. Departs on the warm/editorial vs cool/technical axis.
- **The Specimen Sheet** — warm bone ground, each trace framed and numbered (I–VI) as a catalogued specimen, serif. A monograph plate. Departs on the organizing principle (a catalogue of records, not one hero figure).
The file also implements the rail as a collapsible panel (the full-bleed feedback): collapse the rail to read the plate full-bleed.

**Decided — the identity:** the user reviewed "The Print" against three alternatives (The Console, The Drafting Table, The Specimen Sheet) and chose **"The Print"** as the one that stood out most. It is the series builder's identity and the candidate identity for the whole app: light, warm paper, the figure composed as a publication plate, phase colour carrying the surface.

**The sample table — "The Contact Sheet" (`docs/redesign-mockups/sample-table.html`):**
The sample table in "The Print" language, under a second metaphor that nests inside the first: if the series builder is the finished *print*, the sample table is the **contact sheet** that precedes it. Dark detector frames laid as windows in light paper, screened by eye, marked with a grease pencil (the terracotta accent).

The metaphor resolves the surface's structural tension: a table is row-oriented (good for metadata columns) while culling is image-grid-oriented, and forcing one into the other loses the other. On a contact sheet they were never two jobs. Each row is one sample: a **strip of exposure thumbnails** sitting as a normal table cell, alongside metadata columns (the ordering variable, kept count, status). No mode tax.

- **Two views, one switcher.** *Contact sheet* (all samples, every frame inline; multi-select across rows then batch-drop; double-click opens the loupe) and *Loupe* (one frame at ~500px for spotting flares, with the sample's exposure strip beneath, arrow-key flipping, the metadata sidecar).
- Reject-only, no "accept": every frame is kept until dropped; the grease-pencil ✕ is the only mark. The tool never pre-judges an image — a screened sample's flares are shown dropped because a human dropped them.
- The **representative exposure** (the frame that carries to the Index stage) is a first-class pick — an accent dot on the sheet, a "Set as representative" action in the loupe — making the triage-to-index bridge explicit.
- A **Tags column** carries light, free-form per-sample notes ("what is this?" — the lipid line, the peptide, whether it is a control), shown as quiet pill chips. It is deliberately *not* a structured ordering-variable field: the ordering variable (the 1:0.25 / 1:1 ratios) is collected at series-scoping time on a coherent group, never as an experiment-wide table column. An untagged cell shows a `+ tag` invite — empty is an invitation, not a deficiency. Single-sample tagging is also reachable in the loupe (a "Sample tags" section), the bulk path on the sheet. *(Revised 2026-05-17 — see "Done — sample-table metadata lightened" below.)*
- A row's `screened` state is a hollow vs. filled mark — progress tracking, never a verdict on the images. Built on the same nine-sample LL37 titration; detector frames are SVG-rendered (beamstop, Debye-Scherrer rings, flares).

**The focus workspace (`docs/redesign-mockups/focus-workspace.html`):**
The indexing surface in "The Print" language. Layout and hierarchy only — the trace interactions (peak edit, active-set semantics) are validated and unchanged.

- **The trace is the hero**, a full-width figure plate (Newsreader title = the sample name): the 1D integration with its Bragg peaks marked and Miller-labelled.
- **The detector image is co-resident** in a panel below, never behind a tab — it is the trace's partner, not another page.
- **The q-link is the headline interaction.** A trace peak and a detector ring are the same reflection in two projections; hovering a peak, a ring, or a reflection-table row lights all three (accent). This makes "real signal or artifact?" — the 90% reason the brief says scientists leave the trace — a glance instead of a page round-trip. Built on the insight the user explicitly liked: visualising the rings on the image.
- A **reflections table** (hkl · q · d · ratio) sits beside the detector, the third node of the q-link.
- A sample can hold **more than one phase** (coexistence). Peaks are partitioned among the coexisting phases. The **rail is the output**: the phase call lists *every* phase present, and the candidate indexings are a **multi-select active set**, not radio buttons. Independent phases coexist (Pn3m + Lamellar, disjoint peaks); candidates that explain the same peaks swap (Lamellar ⇄ Hexagonal) — `auto_group` / `remove_subsets` made visible. Toggling a candidate re-labels and recolours peaks and rings live; a peak no active phase explains degrades to "?" / "—". Scoring stays coverage × consistency, per phase.
- A compact **representative-exposure switcher** (mini detector thumbnails) handles the 10% re-choose case in place. **Notes** are a topbar button — reachable, not resident, no standing screen rent.
- Shares components with the sample table (detector renderer, exposure thumbnails, phase chips) — the planned shared-component reconciliation.

**Decided in discussion (sample-table follow-ups):**
- *Process:* one solid mockup per surface first, then refine all three together — shared components (detector renderer, exposure strip, phase chips) and cross-surface flows mean isolated refinement risks divergence.
- *Many exposures:* the sheet row's strip is a teaser (one line, cap with `+N`); the loupe is the unbounded surface (one-at-a-time, arrow-key). Sheet scans samples, loupe scans frames.
- *Two selection scopes:* frame-select (thumbnails → batch reject) and sample-select (a row checkbox → build series / batch metadata) are visually distinct, never merged.
- *Metadata model:* key/value tags underneath; the structured ordering metadata is collected at series-scoping time, not as experiment-wide table columns (see section 3, principle 3). The table carries only light per-sample tags. Bulk entry (fill-down across selected rows, column paste) belongs to the table; single-sample entry to the loupe — shared data entered once, never twice.

**Decided in discussion (focus-workspace follow-ups)** — all three applied to `focus-workspace.html` on 2026-05-17, the first piece of the refinement pass:
- *Candidate hover = preview.* Hovering a candidate indexing previews the peaks it would claim (trace + rings, in the phase's own colour), ephemeral; clicking commits it to the active set. Preserved from the current app, where this behaviour is valued. Reuses the q-link highlight channel. The two hovers are colour-distinct: q-link hover = accent (a cross-reference between views), candidate hover = the phase colour (a preview of an assignment). Hover-preview of a *swap* shows what would be claimed and what would be lost before committing.
- *Detector / reflections equal height.* The lower row is a fixed height keyed to the detector square; the reflections panel matches it and scrolls internally when peaks are many (real samples vary in peak count, so no fragile pixel match). Spare room under a short table carries the foot / per-phase coverage summary.
- *Notes as a width-keyed margin, not a card.* On wide viewports Notes is a narrow, quiet, ruled **margin column** — the marginalia of the plate (fits The Print) — `[ work | rail | notes margin ]`. Below a breakpoint it collapses to the existing topbar toggle / drawer. No middle "reflow" state. It is the first thing to yield as width shrinks, never the trace or the output. Deliberately a margin, not a card, so it does not re-import the "co-equal cards for unequal work" failure. Rationale: "no standing screen rent" was a small-screen concern; on a large monitor the rent is free.

**Done — cross-surface component reconciliation (2026-05-17):**
The three standalone mockups were reconciled so the surfaces read as one app. Because they are deliberately standalone HTML files (no shared import), reconciliation means *consistent output*, not shared code:
- *Phase palette.* One byte-identical `:root` block (`--pn3m` / `--lam` / `--im3m` / `--hex`) carried into all three, marked as the shared block; `--hex` was missing from sample-table and series-builder. The topbar / stage-tabs CSS was already byte-identical and needed no change.
- *Detector renderer.* The instrument constants (centre-brightness coefficient, beamstop dimensions, glow blur) had drifted slightly between sample-table's `detectorSVG` and focus-workspace's `renderDetector`; unified to one canonical set. The colour mode stays per-surface on purpose: triage surfaces (sample table, exposure strips) render rings in phase-neutral gold — the question is image quality — while the focus workspace's big detector renders them in phase colour — the question is identity. Focus-workspace's exposure thumbnails were rebuilt to draw the real ring radii in gold (they had been two generic placeholder rings), matching the sample table's thumbnail language.
- *One coherent dataset.* The three surfaces had used different sample sets; they now share the nine-sample LL37 titration. The sample table lists `smp_01`…`smp_24`, the focus workspace indexes `smp_09` ("sample 4 of 9", was a stray "12 of 39"), and the series builder plots `smp_04`…`smp_18` (its old `C4`/`C5a` codes were the same titration). `smp_09` is the Pn3m+Lamellar coexistence sample in both the focus workspace and the series builder.

All three surfaces now have a mockup in "The Print": series builder, sample table, focus workspace. The broad-strokes pass and the cross-surface reconciliation are complete; what remained then was consolidating the series builder (done, below) and revising the connective-tissue surfaces.

**Decided — architecture revisited (2026-05-17):**
Three bootstrap-era decisions were revisited. All three reframe the same way: a rigid wall — a container, a frozen artifact, an editing permission — gives way to something flexible, and the rigid form is kept only for the job it is genuinely good at. No new hierarchy is introduced — the fix removes a wall, never adds a level.

1. *Experiment: container → provenance facet.* "Everything is scoped to an experiment" was a route / UI convention, not a storage boundary — one SQLite DB already holds many experiments (an `experiments` table; `samples.experiment_id`; the `comparisons` table already carries no `experiment_id`). The top-level scope becomes the whole **corpus** of samples; experiment becomes a chip on each sample and a filter facet. The sample table, series builder, and folio all run corpus-wide; the route layer gains corpus-level endpoints (`/api/samples`, `/api/series`) beside the existing `/api/experiments/{id}/...` — additive, not a rewrite. Experiment stays load-bearing for ingestion (`init` / `reingest` per beamtime directory), for calibration provenance (beam energy / flight path / q-units, stored per experiment), and as the natural default filter right after a beamtime lands. Physics constraint: cross-experiment series are valid because q is absolute (Å⁻¹), but the `q_units` label must be normalized to a canonical unit, and a series spanning beam configurations should say so in its caption — a visible signal, not a wall.

2. *Saved comparison: frozen snapshot → recipe + plate.* Today a saved comparison is a frozen enumeration of pinned exposures (`comparison_members`, each with a non-null `snapshot` JSON; git-style `forked_from_id`). It has no living / criteria layer — which the autogrouping principle (machine proposes the grouping, human confirms) needs. A series gains two layers: a **recipe** (the criteria — which samples, which ordering variable, the order rule — that the machine maintains and proposes, with manual pins / excludes overriding it) and a **plate** (today's frozen, hashable, forkable snapshot — the artifact you export and publish). You compose on the recipe; you snapshot to publish. The `comparisons` table reshapes into `series` and gains the recipe.

3. *Series editing: author gating → ungated, the recipe/plate split does the work.* A bootstrap instinct gated series editing to an author. Dropped — for the same reason as 1 and 2: it is a wall guarding a room already locked from the inside. Everything gating protects is already covered. *Clobbering:* every series mutation flows through the event log (`apply_event!` / `user_actions`), reversible like any other — nothing is destroyed. *Attribution:* `user_actions` already records who did what (the deferred per-user audit view surfaces it). *"Whose figure":* the **plate** is frozen by construction — no permission is needed to keep it from changing — and carries a snapshot stamp ("snapshot by X · date"), which is what makes ungated recipe-editing safe to publish from. *"I want my own version":* `forked_from_id` already exists. The recipe is a hypothesis under test, and hypotheses get better when the whole lab can poke at them; a gated series would also be the one author-locked object in an app whose indexing, peak curation, and phase calls are all already ungated and multiplayer (Plan 7) — an inconsistency with no justification. The ownership-flavoured cues that remain are soft, not walls: a Draft card defaults into its creator's folio view (a sort / visibility default, not a lock); a plate carries its snapshot author; `archived` is a state anyone can set, not a permission. If HimalayaUI ever goes multi-tenant, the wall that earns its place is tenancy (whose corpus / whose DB), not authorship of a series — and not now, given the per-experiment `serve` deployment and the single-lab unit.

*Noted, not revamped:* an exposure is hard-coupled to one on-disk file today; the deferred derived-exposure feature will make an exposure "have a source" rather than "be a file" — so the triage redesign must not assume every frame is a raw detector image. The manifest stays an optional bulk-import path, never the schema authority (already section 3, principle 3 — the same spirit).

**The series folio (`docs/redesign-mockups/series-folio.html`):**
The home of the Series stage — a corpus-wide wall of saved series, the first surface built under the architecture decisions above.
- *A masonry of plate-cards.* Each saved series is a card whose top is a live miniature of its own waterfall. A 6-sample series is a taller plate than a 3-sample one, so the wall is a masonry of distinct figures — the "identical card grid" anti-pattern is escaped *because* the figure is the artifact.
- *The phase strip.* Below each mini-figure, one cell per sample coloured by its phase call (coexistence cells are a two-phase gradient): the series' phase story made scannable across the whole wall, with a caption ("Pn3m → Lamellar", "Pn3m throughout", "No clear phase").
- *The recipe surfaces.* A series carries a recipe; when the machine finds samples that match but are not yet in the series, the card shows a "+N new match" pill. A series still being scoped shows as a dashed **Draft** card ("Recipe" kicker, figure dimmed) — a recipe with no committed plate yet.
- *Corpus-wide.* No experiment in the topbar crumb (a "Himalaya · SAXS" wordmark instead); each card footer carries its beamtime as provenance; a cross-experiment series shows "⇄ April + July · q normalized" — the physics honesty note from the architecture decision made visible.
- *Exploration.* Live search, filter chips (all / has transition / cross-experiment), and sort (recent / variable / largest); the header count reflects the filtered subset; an empty state when nothing matches.
- `+ New series` is the one primary action — it, and sample-table multi-select, lead into the scoping step.

**The series-scoping step (`docs/redesign-mockups/series-scoping.html`):**
The surface between the contact sheet and the series builder — the Series stage's other front door, where a new series' *recipe* is born and the ordering metadata becomes structured data. A review surface, not a form: you arrive with samples already selected, and Himalaya has grouped them, read the ordering variable from their names, and parsed each value.
- *Machine proposes, human confirms.* A confident one-line autogroup summary; the proposed **ordering variable** as the single editable decision; the samples already ordered low-to-high, each with a trace sparkline (the cubic → lamellar transition is visible running down the rows).
- *Attention only where needed.* Confidently-parsed values sit calm in ink; a value Himalaya is unsure of is flagged (amber, dashed, "check the read") with a faint wash on its row. Clicking it accepts the read. The footer tracks the state ("1 value to check" → "All 6 values confirmed — ready to build"), the same attention-routing as triage.
- *Metadata as byproduct.* The footer says it plainly — confirming records the ordering variable on every sample, so the next series that needs it already knows. The redesign's progressive-metadata principle made concrete: confirmation, never a form.
- *The machine reasons about coherence.* A "Himalaya also found" section surfaces samples that match loosely but were not assumed in — e.g. a sample on a different lipid line — with the machine's own counsel ("its own series?"). It proposes; it does not presume.
- *Confirm & build* creates the recipe and opens the builder — **gated**: the button stays disabled while any value is still flagged. A deliberate bit of friction (warranted here) so an uncertain parse cannot be rubber-stamped into the series; the footer says why ("1 value to check before you can build"). Clearing every flag is always possible, so the gate is satisfiable.

With the folio and the scoping step both built, the Series stage is complete in mockup: folio → scoping → builder.

**Decided in discussion (scoping-step follow-ups):**
- *Parsing robustness is a graceful-degradation property, not a parser quality.* Name parsing alone is unreliable (section 3 says so) — so the scoping surface is designed to be correct at *any* parser quality: a good parser means a one-click confirm, a poor one means every row flagged and the surface degrades to a guided entry form, never broken. The parsing *technique* (regex / fuzzy / LLM) is therefore an implementation detail. Three design choices carry even a mediocre parser: parsing happens at scoping time on a coherent *group* (the delta across sibling names is far more legible than any one name in isolation); confirmations become structured `sample_tags`, so the parser's job shrinks as the corpus fills (self-improving); and signals beyond the name (filename patterns, manifest columns, prior tags) feed it. The failure mode to guard against is *confident-wrong* — an unflagged bad parse the human rubber-stamps — mitigated by flagging generously, sanity-checking the resulting ordering (duplicates, broken monotonicity), and reversibility.
- *Confirmations are reversible.* In-session, nothing is durable until "Confirm & build": resolving a flag only dismisses the machine's question, every value stays re-openable (click to re-open, or Undo / ⌘Z steps back the last change), and "Discard" throws the whole session away. After build, a confirmation is a `sample_tag` written through the event log (`apply_event!`, `user_actions`) — reversible the same way every other mutation is, and non-destructive (a value set, always visible wherever the sample appears, always re-settable). The deferred per-user audit view is the surface for browsing and reverting; the log already supports it. Applied to `series-scoping.html` (re-openable values, in-session Undo).

**Done — series builder consolidated (`docs/redesign-mockups/series-builder.html`, 2026-05-17):**
The three exploration files collapsed into one canonical mockup; the variant verdicts were the merge rule.
- *Offset waterfall is the default* — the publication figure the redesign exists to produce: phase-coloured traces, peak ticks, Miller labels on the lamellar trace, the transition bracket up the left margin.
- *Heatmap is a representation toggle* — a Waterfall / Heatmap segment in the rail. The heatmap renders the same series as a q–ratio intensity map (one row per sample in its phase hue, brightness is scattered intensity, peak migration reads as a streak); the offset and tracking controls fold away because they do not apply. Kept on the verdict that it surfaces weak peaks well.
- *Full-bleed is the rail's collapsed state, not a variant* — a collapse control folds the rail away (animated grid column), the plate widens to own the screen, a floating "Compose" tab restores it, and a floating offset dock keeps the one live waterfall control reachable while reading. The verdict said full-bleed is a mode; it is now a mode.
- *Peak tracking is an optional annotation layer* — a toggle, off by default (the verdict: useful, but reads best when peaks genuinely align). On, it links each Miller order across the traces that carry its phase and calls out one d-spacing per phase.
- *Drenched was dropped* — it carried no verdict and was redundant: phase colour already carries the plate without a colour-field treatment behind it.
- Built on the reconciled nine-sample LL37 titration (`smp_04`…`smp_18`), corpus wordmark, the shared phase palette — consistent with the other four mockups. With this, all five surfaces have one canonical mockup each.

**Done — sample-table metadata lightened (`docs/redesign-mockups/sample-table.html`, 2026-05-17):**
The contact sheet's metadata column was the last surface still carrying the bootstrap-era inferred-column model — a structured "LL37 : lipid" ordering-variable column with confirmed (ink, solid underline) vs. inferred-from-name (faint, dashed, "read from name · confirm") values. That model was superseded by section 3 principle 3 (structured ordering metadata is collected at series-scoping time, not as experiment-wide columns) and by the scoping step itself, which *is* the machine-proposes / human-confirms surface — the table was doing that job a second time, less well, and on a now corpus-wide table where one named variable column does not cohere.
- The structured column is replaced by a **Tags column**: light, free-form per-sample notes ("what is this?" — lipid line, peptide, control), quiet pill chips. Tags are pure identity; they deliberately do not encode the ordering variable.
- No confidence states, no "read from name · confirm", no machine inference shown in the table. An untagged cell shows a `+ tag` invite; with tags present the `+` is a quiet appendix revealed on row hover (the same hover-reveal language as the sample name).
- The loupe gains a "Sample tags" section — single-sample tagging belongs there, the bulk/fill-down path belongs to the sheet (the two-scopes rule from the sample-table follow-ups).
- The header copy now says tags are a light, optional note and that the ordering variable is named later, at scoping. This is the same un-rigidifying move as the three architecture decisions: a structured field becomes a light facet, and the structured form is kept only where it genuinely belongs (the scoping step).

**Open / next:**
- ~~The refinement pass across all three surfaces together — reconcile the shared components.~~ *Done 2026-05-17* (focus-workspace follow-ups + cross-surface reconciliation, both above).
- ~~Surface and explore saved comparisons — the series folio.~~ *Done 2026-05-17* (corpus-wide, living-recipe model; entry above).
- ~~Design the intermediate series-scoping step (where the ordering metadata is collected).~~ *Done 2026-05-17* (`series-scoping.html` — machine-proposes / human-confirms, metadata as byproduct; entry above).
- ~~Propagate experiment-as-facet to the older three mockups.~~ *Done 2026-05-17* — all five mockups now carry the corpus wordmark; experiment shows as a facet where it belongs on each surface: a beamtime filter chip in the sample table's topbar, a provenance segment on the sample in the focus workspace (`smp_09 · SSRL Apr 2026 · …`), and a provenance tag on the series in the builder. Stale "in this experiment" copy was retired.
- ~~Consolidate the series builder on "The Print": the offset waterfall as default, heatmap as a toggle, full-bleed as a collapsible-rail mode.~~ *Done 2026-05-17* (`series-builder.html` — offset waterfall default, heatmap toggle, full-bleed rail-collapse, peak tracking as an annotation layer; entry above).
- ~~Revise the sample table's metadata treatment toward series-scoped collection (lighten the inferred-column model that the current mockup shows).~~ *Done 2026-05-17* (`sample-table.html` — structured inferred-variable column replaced by a light free-form Tags column; entry above).
- ~~Translate the redesign into an implementation plan against the real codebase.~~ *Done 2026-05-17* — master plan at [`docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`](superpowers/plans/2026-05-17-himalaya-ui-redesign.md), grounded in a file-anchored survey of the current backend and frontend, then **hardened to v5 across five four-reviewer passes** (backend/SQLite, mutation queue, frontend, architecture) — the fifth came back clean from all four; the plan is review-converged. Six execution phases: Phase 0 corpus listing routes → 1 sample table → 2 series backend → 3 series stage UI → 4 focus workspace → 5 final cutover, with per-surface cutover (each phase retires what it replaces). Headline reviewer finding: the `comparisons`→`series` change is a *copy into new tables*, not a rename — a rename would break event-log replay, so the `comparison*` tables and their dispatcher branches are permanent replay machinery; the copy synthesizes events and folds them through pure-replace dispatcher branches. Phases 0–2 are specified to actionable detail; 3–5 at blueprint level, each to get its own detailed plan when picked up. §11 carries the consolidated reviewer constraints.
- The manifest redesign: deferred, app-informed.
- Superseded: an earlier framing offered two directions, "The Bench" and "The Light Table." Both were re-weightings of one screen and are superseded by the three-surface model, which is derived from the workflow.

**Decided — focus-workspace mockup-vs-shipped divergences retained (#209, 2026-05-27):**
Three items the round-2 fidelity audit (`docs/2026-05-27-the-print-round2-findings.md` R4-N2 / R4-N3 / R5-N1) flagged as "verify intentional" — shipped features the mockup does not show. Reviewed during the milestone-3 focus finish; each is **retained** as a deliberate divergence between the mockup (layout study only) and the shipped surface (real features with real data shapes). The mockup is the visual north star, not a feature-completeness contract.

- *Speculative section retained (R4-N2).* `PhasePanel.tsx:349-375` ships a "Speculative" rail section above-the-fold with an "+ Add speculative" button and a delete affordance on each speculative row; the mockup's rail has only Phase call + Candidate indexings + the rail-note. Kept because the speculative builder is a shipped, validated feature — users build sub-minpeaks indices when the auto-pass misses a phase. Collapsing it to an inline "+ Add speculative" at the foot of Candidate indexings would hide existing speculative rows, breaking the "rail is the output" rule. The empty-state pattern (button visible, no rows) is the same shape the Compose menus use elsewhere.
- *TitleStrip q-range / reset / Copy / Save retained (R4-N3).* `PlotCard.tsx:457-480`'s TitleStrip keeps numeric q-range inputs (min/max + reset) and the figure-export Copy / Save cluster; the mockup `.tools` strip is only `[log | linear] · Auto-fit · + Peak`. Kept because these controls back shipped features — figure export is the analytical handoff that the focus workspace exists for (the same PNG/SVG controls live on the Compare page), and the q-range numeric inputs are the precision affordance for setting bounds when the auto-fit window is wrong. Folding them into a popover was considered and rejected: a popover would hide the export entry-point, the highest-signal action on the surface. The mockup is a layout study, not a feature inventory; the strip earns its width.
- *FocusNotesMargin single-textarea simplification retained (R5-N1).* `FocusNotesMargin.tsx:33-49` ships one focus-gated textarea writing to `sample.notes` (a single string); the mockup `.notes-margin` shows kicker + count badge + each note as a `.note` row with author/timestamp + dashed "Add a note…". Kept because the data model is a single string per sample — there is no per-note author / timestamp / reference field on the backend. Upgrading would mean either (a) inventing a notes thread schema (out of scope for #209) or (b) faking the markup against a single-string source (a visual lie). The textarea is the honest projection of the shipped data shape; the mockup is the projection of a richer model the redesign has not committed to. If the deferred per-user audit view ever lands a notes-thread schema, this becomes the natural rendering target — until then, single textarea matches the data.

All three are the same shape of decision: when the mockup and the shipped data shape disagree, the shipped data shape wins. The mockup informs visual hierarchy; the data shape informs surfaces. Recorded here so the next fidelity pass does not re-litigate them.
