# Plotting redesign — design spec

> **Status:** brainstormed 2026-05-29, first-principles pass. Captures the converged *concept* for HimalayaUI's plotting surfaces (Focus single-pattern + Series cross-sample). Visual craft (high-fidelity rendering) and the implementation plan are downstream of this doc, not in it.
>
> **Companions:** `PRODUCT.md` (strategic), `DESIGN.md` ("The Print" visual system), `docs/2026-05-29-frontend-audit-impeccable.md` (the audit that scoped this as program item #2), `docs/redesign-notes.md` (the redesign history this extends).

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
The product's hard colour-blind rule, applied where colour is densest:
- **Colour = phase identity, only.** Never provenance, never state.
- **Silhouette = provenance:** auto peak = filled triangle; manual peak = a distinct glyph (e.g. diamond). Retires magenta-for-manual (frees colour for phase).
- **Fill vs outline = real vs in-flight** (optimistic peak = outline-only; already the precedent).
- **Struck glyph = excluded** (hollow + terracotta strike — the grease-pencil negation), replacing today's invisible "faint gray".
- **Separate mark = predicted-but-unobserved** reflection (a caret/tick, never a triangle).
- Transient interaction state (q-link hot, candidate preview, dim) rides **halo + opacity**, never silhouette.
- Series peaks use the same atoms but **lightly** (anchors), per §4.

### 5.2 One mark-builder, one plot core
- **`peakMark()`** — a single builder consumed by Focus on-screen, Series on-screen, **and both export paths**. Collapses today's five drifted definitions (triangle direction, offset px) into one. The *only* allowed divergence is the deliberate Print-vs-scientific style split at the export boundary (§4.2).
- **`<PlotSurface>` primitive** — owns the Observable Plot instance, shared log/linear x-scale, gestures (wheel-zoom / brush / dblclick), invert/apply, and a **rAF-batched resize** (kills the per-frame setState replot in TraceViewer / MillerPlot / MultiTracePlot). Exposes `marks` / `overlay(scales)` / `hitTest()`. Both plot cores consume it — retires ~400 duplicated LOC.

### 5.3 Interaction grammar (still open — defer to the visual/build pass)
The audit's grammar problems (overloaded single-click, theatre "+ Peak", two conflicting resets) remain to resolve. The earlier brief offered "select-then-act" vs "fix-the-lies"; this was **not** settled in the first-principles pass and should be decided during the high-fidelity build, where the real cursors/affordances can be felt. Recorded as an open decision, not a closed one.

### 5.4 Preserve (audit's preserve-list)
`phases.ts` / `coloring.ts` palette engineering (AA-pinned, single source), the serif/sans/mono voice discipline, `formatAxis`, the PlotCard auto-fit heuristic, the **q-link triple** (the best-composed thing in the app — *enrich* it with phase-colour rings, don't regress it), and the modal focus-trap.

## 6. What this is NOT

- **Not** an automated phase-caller. The judgment stays human; the tool ranks, shows, and does the Bonnet arithmetic.
- **Not** a figure builder. The Series surface thinks; export is a fixed clean render.
- **Not** an ultra-minimal surface. The supporting instruments (detector, combs, residual) stay co-present.
- **Not** a layout revolution on Focus — the composition lands close to today; the substance is the four moves in §3.

## 7. Open items / next

1. **High-fidelity visual pass** — the immediate next step. Render both surfaces in real Print tokens/fonts/trace quality (impeccable's high-fidelity presenter), judge and tune *craft* (the schematic brainstorm medium can't). This spec is the brief for it.
2. **Interaction grammar** (§5.3) — decide select-then-act vs fix-the-lies during that pass.
3. **Export default** — pin the one sensible default (line weight, colour mode, label density) during the visual pass; ship Copy/Save only.
4. **Implementation plan** — after visuals settle, via the writing-plans flow: `PlotSurface` + `peakMark()` first (the spine), then Focus (assignment/cart, combs, Bonnet, phase-colour rings), then Series (waterfall thinking-view + decoupled export), with the live functional traversal as the acceptance gate.
