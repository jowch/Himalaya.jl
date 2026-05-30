# Plotting redesign — design spec

> **Status:** brainstormed 2026-05-29 (first-principles concept); **high-fidelity visual pass complete + user-approved 2026-05-30**. The concept below is now realised as two interactive high-fidelity mockups; §8 records what the visual pass built, the decisions it settled, and the handful of points where the *built* design deviates from the original concept (encoding atoms especially — read §8 before §5.1). The implementation plan is the next step, downstream of this doc.
>
> **Mockups (the visual contract for implementation):** `docs/redesign-mockups/2026-05-29-focus-plot.html`, `docs/redesign-mockups/2026-05-29-series-plot.html` — real Print tokens/fonts, analytic SAXS traces, real lattice math; every interaction in §8 is live in them.
>
> **Companions:** `PRODUCT.md` (strategic), `DESIGN.md` ("The Print" visual system), `docs/2026-05-29-frontend-audit-impeccable.md` (the audit that scoped this as program item #2), `docs/redesign-mockups/2026-05-30-plotting-mockup-audit.md` (the impeccable audit of the mockups → the implementation a11y/responsive checklist), `docs/redesign-notes.md` (the redesign history this extends).

## 1. Why this pass

The plotting surfaces are two independently-grown plot engines with five disjoint definitions of "how a peak looks" (already drifted), a colour encoding that breaks the product's hard colour-blind rule, and an interaction grammar that's implicit and overloaded. The earlier `impeccable` brief proposed fixing those *as the mockups framed them*. This pass re-asked the question from first principles and changed the **concept**, not just the execution.

**The governing realisation.** The current design makes the **intensity curve** the hero, but the indexing *judgment* is about **peak positions and their ratios**, and — more importantly — **all quantitatively-plausible candidate indexings fit well**. The hard question is not "does this index fit the peaks" (they all do) but **"which of the plausible indexings actually makes sense?"** — a domain-judgment call that lived entirely in the expert's head in the old UI.

The design response: **do not automate the judgment** (the data is too ambiguous; a tool that decides will be confidently wrong). **Automate the exact, tedious physics** (the Gauss–Bonnet lattice ratio) and make every *way of looking* fast, so the human's judgment is well-served, not replaced. This is the "thinking tool, not oracle" line from `PRODUCT.md` principle 2.

## 2. Two surfaces, one vocabulary

Focus and Series are **different jobs** and get **different compositions**, but share the same **atoms** (peak marks, phase colour, q-axis language) so the app reads as one instrument.

| | **Focus** (one pattern) | **Series** (many patterns) |
|---|---|---|
| Question | "Are these peaks real, and which phase(s) make sense?" | "What changes across the variable?" |
| Subject | **The peaks** (curation) | **The traces** (comparison) |
| Peaks rendered | densely, every one (you curate them) | lightly, as interaction anchors + interpretive accents |
| Output | an **assignment** (phase call) per exposure | a shareable **figure** (decoupled export) |

## 3. Focus surface

**Composition stays close to today's focus workspace** (validated layout); the substance is conceptual, in four moves.

### 3.1 The "assignment" model (replaces "active set")
Each exposure has an **assignment** — a cart that starts **empty** (form factor / no phase) and holds **0..N phases**. You **evaluate candidates and pick** them in, like a shopping cart. The swap/coexist machinery (`auto_group` / `remove_subsets`) still informs ranking under the hood but is **not** surfaced as a taxonomy the user must navigate. Vocabulary: "assignment", not "active set".

Two distinct add affordances:
- **`add` on a candidate row** → drop that ranked candidate into the assignment.
- **`+ custom index` next to the assignment** → build a custom / speculative index yourself (the existing speculative builder), for when no candidate is right.

### 3.2 Ranked, scannable candidates + the one automation
The auto-indexer still produces candidates; they're presented as a **ranked, scannable list** (plausible ones on top), each fast to evaluate (phase dot, score bar, lattice `a`). Human scans and decides — this already worked in the old UI; the tool's job is good ranking + fast evaluation, not auto-deciding.

**The single smart automation: the Gauss–Bonnet ratio.** Once a bicontinuous cubic (Pn3m / Im3m / Ia3d) is in the assignment, the lattice parameter of a coexisting cubic is *predicted* by the Bonnet ratio. Candidates whose measured `a` matches get a **score bump and a `⭙ Bonnet` flag**. This is the principled "in context with the phase I already chose" help — physics, not peak-coverage heuristics. **Bonnet is the only such automation in scope**; the design leaves room for more exact constraints later but commits to none.

> **Explicitly rejected as over-reach** (recorded so the next pass doesn't re-litigate): a "leftover-peaks-drive-the-next-phase" flow and a "combination check" engine. Coexistence is messy — one phase dominant and obvious, the second subtle, lower-signal, with peaks partly buried under the dominant phase's. Mechanising the second-phase call pretends a decidability that isn't there.

### 3.3 Reflection combs replace the reflections table
The reflections table becomes a **comb view**: per assigned phase, predicted reflections are comb teeth aligned to a q position ruler; teeth landing on observed peaks = explained, gaps and extras tell the story. This is the faster "does this assignment explain the peaks" read the table couldn't give. The comb panel toggles to the **indexing-space residual** view (transform the axis so a correct phase falls on a straight line; a deviating peak bends off it — good for spotting one bad peak). The comb tooth is a **third q-link node** (peak ↔ ring ↔ comb tooth).

### 3.4 The detector image stays first-class
The detector is **the real source**, not a supporting strip. Sometimes a feature is clearer on the 2D image than the integration; other times you go to the image to rule out an artifact. Co-resident, q-linked (hover a peak/ring/comb tooth → all light), rings rendered **in phase colour** (restoring the mockup intent the shipped code regressed to uniform gray); neutral when unassigned, terracotta on q-link hot.

### 3.5 Instruments are co-present, the human routes between them
Trace ("is the signal real and sharp?"), detector ("or an artifact?"), combs ("does this assignment explain the peaks?"), indexing-space ("does any peak deviate?"), rail ("which do I commit?"). The tool never collapses these into one verdict; it makes each fast. Not an ultra-minimal surface — rich and co-present.

## 4. Series surface

**Waterfall is the hero.** The stacked-offset waterfall, ordered by the series variable, is the primary view. The heatmap (q–ratio intensity map) **demotes to a secondary toggle** — good for weak-peak/migration reading, but not the common idiom, so it doesn't lead.

**Traces are the subject; peaks are light.** Unlike Focus, the Series job is *comparing traces to see what changes* across the variable. Peaks render **lightly** — small anchors, not dense curation triangles — and serve as **interaction handles** (hover → the phase/indexing for that trace lights up) and interpretive accent. Phase is **interpretive detail**, carried by trace colour and a slim **phase-strip companion** (one cell per sample along the variable, coexistence as a two-phase gradient, with a transition caption).

### 4.1 The on-screen view is a thinking instrument, not a figure builder
Decisive call (and the user's instinct): **the Series surface is not a WYSIWYG figure builder.** It is a thinking instrument in the Print aesthetic, optimised solely for reading what changes. The moment it tries to *be* the publication figure, it inherits the audit's "figure competes with its chrome" failure and invites pixel-fiddling instead of thinking.

### 4.2 Export is a decoupled, zero-config, clean scientific render
"Publication-ready" here means **share with a colleague who doesn't use Himalaya, or drop into a slide** — not a journal submission. The export is therefore a **separate renderer**, deliberately **not** a screen-grab and **not** the Print brand:
- **Style:** clean scientific idiom (GraphPad / Origin) — white background, bold lines, labelled axes (`q (Å⁻¹)`, `Intensity (a.u.) + offset`), sample labels, title, provenance footnote.
- **Config:** **none.** One sensible default that handles most plots. One-click **Copy** (and a Save), no tuning UI. Predictability comes from the export rendering the **same plot elements** as the on-screen view (same samples, order, offset, phase annotation) — so you know what you'll get without composing anything. Tuning can be added later only if a real need appears.

This turns the audit's "two source-of-truth for a mark" liability into a **feature**: the on-screen renderer and the export renderer are *meant* to differ in style (Print vs scientific), so the duplication is correct-by-design — same data, two intended looks.

## 5. Shared atoms & architecture

These underlie both surfaces and resolve the audit's Themes B (encoding) and C (grammar) + the architecture finding.

### 5.1 Encoding: colour = phase, everything else gets a non-colour channel
The product's hard colour-blind rule, applied where colour is densest. **The visual pass refined three of these atoms — the list below is the settled, as-built form; see §8.1 for the rationale on each change:**
- **Colour = phase identity, only.** Never provenance, never state.
- **Silhouette = provenance:** auto peak = filled **downward** triangle (a marker pointing *at* the peak, not a "mountain" sitting on it — settled in the visual pass); manual peak = a diamond.
- **Fill vs outline = real vs in-flight** (optimistic peak = outline-only; already the precedent).
- **Excluded = ghosted** (hollow outline + reduced opacity, no label). **The "terracotta strike" from the original concept was cut** — it read as redundant and was hard to see; a hollow faint glyph distinguishes excluded from unindexed-but-real (filled gray) without the strike.
- **Separate mark = predicted-but-unobserved** reflection: a **hollow caret** in the Focus combs, a **hollow ghost ring at the predicted q** in the Series migration track (§4 / §8.2). Same atom, both surfaces; never a triangle.
- Transient interaction state (q-link hot, candidate preview, dim) rides **halo / size-grow + opacity**, never silhouette. *Note:* the q-link "hot" state **grows the mark and adds a terracotta ring rather than recolouring it terracotta** — terracotta (hue 38) sits too close to Pn3m amber (58) for a hue-only highlight to read; size + ring is hue-independent.
- **`peak-manual` magenta is retired (settled 2026-05-30).** An unindexed manual peak is a **neutral (unindexed-gray) diamond**, not magenta — gray = "no phase yet", the diamond silhouette = "manually added"; once indexed it takes the phase colour. This keeps colour strictly = phase (magenta was colour-for-provenance, a rule break) while silhouette alone carries auto-vs-manual, so it survives colour-blindness. The `--color-peak-manual` token is dropped from the plotting surfaces.
- Series peaks use the same atoms but **lightly** (anchors), per §4.

### 5.2 One mark-builder, one plot core
- **`peakMark()`** — a single builder consumed by Focus on-screen, Series on-screen, **and both export paths**. Collapses today's five drifted definitions (triangle direction, offset px) into one. The *only* allowed divergence is the deliberate Print-vs-scientific style split at the export boundary (§4.2).
- **`<PlotSurface>` primitive** — owns the Observable Plot instance, shared log/linear x-scale, gestures (wheel-zoom / brush / dblclick), invert/apply, and a **rAF-batched resize** (kills the per-frame setState replot in TraceViewer / MillerPlot / MultiTracePlot). Exposes `marks` / `overlay(scales)` / `hitTest()`. Both plot cores consume it — retires ~400 duplicated LOC.

### 5.3 Interaction grammar — settled: select-then-act
**Settled in the visual pass (2026-05-30).** The mockups use **select-then-act** throughout and it was approved: clicking a candidate adds it to the assignment; clicking a peak/ring/comb-tooth selects that reflection (q-link); hovering previews without committing (candidate → plot preview; anchor → migration track); the explicit `remove` ×, `+ Peak` arm-toggle, and `+ custom index` are distinct, labelled affordances rather than overloaded clicks. "fix-the-lies" is not pursued. The audit's specific grammar defects (overloaded single-click, theatre "+ Peak", two conflicting resets) are resolved by giving each action its own affordance; the implementation must preserve that one-action-one-affordance discipline rather than re-collapsing them.

### 5.4 Preserve (audit's preserve-list)
`phases.ts` / `coloring.ts` palette engineering (AA-pinned, single source), the serif/sans/mono voice discipline, `formatAxis`, the PlotCard auto-fit heuristic, the **q-link triple** (the best-composed thing in the app — *enrich* it with phase-colour rings, don't regress it), and the modal focus-trap.

## 6. What this is NOT

- **Not** an automated phase-caller. The judgment stays human; the tool ranks, shows, and does the Bonnet arithmetic.
- **Not** a figure builder. The Series surface thinks; export is a fixed clean render.
- **Not** an ultra-minimal surface. The supporting instruments (detector, combs, residual) stay co-present.
- **Not** a layout revolution on Focus — the composition lands close to today; the substance is the four moves in §3.

## 7. Open items / next

1. ~~**High-fidelity visual pass**~~ — **done + approved 2026-05-30.** Both surfaces rendered in real Print tokens/fonts/trace quality; craft judged and tuned over several review rounds. Outcomes in §8.
2. ~~**Interaction grammar**~~ — **settled: select-then-act** (§5.3).
3. ~~**Export default**~~ — **pinned** (§8.2): clean scientific idiom (white ground, 2px bold traces, Arial axis labels `q (Å⁻¹)` / `Intensity (a.u.) + offset`, title, provenance footnote), same samples/order/offset as the view; Copy PNG / Save SVG, no tuning UI.
4. **Implementation plan** — **the immediate next step**, via the writing-plans flow. Build order: `PlotSurface` + `peakMark()` first (the spine), then Focus (assignment cart, ranked candidates + Bonnet, combs + indexing-space residual, hypothetical preview, custom-index builder, phase-colour detector rings), then Series (waterfall thinking-view + phase-strip + migration tracking + decoupled export). The live functional traversal is the acceptance gate. Fold in the mockup audit (`docs/redesign-mockups/2026-05-30-plotting-mockup-audit.md`) as the a11y/responsive checklist: button semantics for the clickable rows, keyboard paths for the q-link / migration track, responsive breakpoints + 44px touch targets, landmark structure, `text-ink-faint` contrast check, em-dash copy sweep.

## 8. Visual-pass outcomes (2026-05-30)

What the high-fidelity pass built and settled. The two mockup files are the visual contract; this section is the prose index to them. Items marked **(new)** were not in the §1–§6 concept and emerged during the pass.

### 8.1 Focus — as built
- **Trace hero** with the §5.1 atoms live: down-triangle auto peaks, magenta diamond manual, ghosted (hollow/faint) excluded peak, `?`-labelled unindexed leftovers. q-line drops from the hovered mark; q-link "hot" grows the mark + adds a terracotta ring (not a recolour).
- **Assignment cart** (right rail): 0..N phase blocks, each with phase name (serif, phase-coloured), score, lattice meta, score bar, and a clean **× remove** in the block header. `+ custom index…` footer reads as a button (hover fill). Contextual note below the cart — appears **only when substantive** (the Bonnet suggestion when Pn3m is in and peaks are unexplained; the Bonnet consistency check `ratio ≈ 1.279` when Pn3m+Im3m coexist); no filler line otherwise.
- **Ranked candidates** with the `⭙ Bonnet` flag + score bump on the coexisting cubic whose lattice matches the Gauss–Bonnet ratio.
- **Combs panel** replacing the table: per-phase teeth on a shared log-q ruler, full-height gridlines tying teeth to observed peaks, hollow-caret for predicted-but-absent, a **leftover row** of unexplained observed peaks. Toggles to the **indexing-space residual** (lollipop) view where deviating peaks bend off the Δq/q = 0 line.
- **Detector** first-class, rings in phase colour (two concentric ring sets under coexistence), q-linked.

### 8.2 Focus — new interactions (new)
- **Hypothetical-assignment preview.** Hovering a candidate previews what adding it would do — rendered **on the plot only** (a dashed ghost comb row + trace highlight), deliberately **not** in the cart (cart mutation on hover caused reflow/jank). The ghost row honestly shows a weak candidate explaining few/zero peaks (all-caret row), which is the disambiguation value.
- **Custom-index builder.** `+ custom index…` opens a modal: symmetry segmented control (Pn3m/Im3m/Ia3d/Lamellar/Hexagonal), a lattice slider + numeric input, and a **live preview comb** with a running **"lands on N of M observed peaks"** fit. Reflections from real physics (`2π√N/a` cubic, `2πn/d` lamellar, `4π√M/(√3·a)` hex). **Snap-to-peaks:** dragging magnetically pulls the lattice to values where a predicted order lands exactly on an observed peak (≈2.6 Å zone, "snapped" indicator), and clicking an observed peak snaps the first reflection onto it; free positioning between snaps. **Add** commits it as a first-class candidate (with swap). *Implementation note:* the snap threshold and density are a one-line tune; revisit if it feels too aggressive on real data.

### 8.3 Series — as built
- **Waterfall hero** ordered by the variable; **heatmap** as a secondary toggle. Per-sample traces from real phase models with drifting lattices, so peaks genuinely migrate.
- **Light anchors as interaction handles (new behaviour).** Anchors are plate-ringed beads (detectable but light). Hovering one **tracks that reflection (phase + order index) across the whole stack**, threading a terracotta connector through every trace it appears in — the migration the waterfall exists to show. A **missing peak** (phase present, that order not observed) renders a **hollow ghost ring at the predicted q** on the baseline and the connector routes through it (same predicted-but-absent atom as the combs); where the phase is *absent entirely*, the track simply ends. Keying on phase+order (not nearest-q) is what makes "absent" drawable.
- **Phase-strip companion** (one cell per sample, coexistence = two-phase gradient, form-factor = hollow dashed cell).
- **"Phases present" reading — derived, not narrated.** Computed from per-sample phase calls: each phase's variable span + lattice trend (`a 205 → 195 Å`), plus `coexistence at …` / `form factor only at …` lines. Generalises to any series; no hand-written "X → Y" story.
- **Member rows** enriched with lattice (`a`/`d`, both under coexistence) + first-peak `q₁`; phase names coloured by phase so coexistence rows self-decode.
- **Form-factor-only members** handled end-to-end: neutral broad-shouldered trace, no anchors, hollow strip cell, "no Bragg peaks · q₁ —" row.
- **Decoupled export** (§4.2): one-click modal rendering the same samples/order/offset in the clean scientific idiom (white, 2px traces, Arial axes/title/footnote); Copy PNG / Save SVG; no tuning.

### 8.4 Deviations from the §1–§6 concept (reconciled)
1. Excluded peak: **ghosted, not struck** (§5.1).
2. Auto peak silhouette: **downward** triangle (§5.1).
3. q-link hot: **grow + ring**, not terracotta recolour (§5.1) — hue-proximity to Pn3m amber.
4. Hypothetical preview: **plot-only**, never the cart (§8.2) — anti-reflow.
5. `peak-manual` magenta **retained** for unindexed-manual peaks (the concept retired it) — open nuance in §5.1, recommend keep.
