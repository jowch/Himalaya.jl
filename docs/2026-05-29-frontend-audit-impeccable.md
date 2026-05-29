# Frontend audit — impeccable (2026-05-29)

> Diagnostic pass over the HimalayaUI React/TypeScript frontend, run as the first step of a
> three-part program: **audit → component-library extraction → trace-plotting redesign**. Four
> independent blind assessments (design-system, plotting, cross-cutting technical, holistic
> design-director) plus the deterministic `impeccable detect` scan. This doc is the consolidated,
> prioritized findings set; it feeds the two build specs that follow it.
>
> Scope note: functional/interaction debt (inert controls, unreachable corpus→sample path, disabled
> Index tab) is **not** re-litigated here; it lives in `docs/2026-05-29-functional-interaction-sweep.md`.
> This audit judges the **design, visual, structural, a11y, and consistency quality of what renders.**

## Verdict first: does it look AI-generated?

**No.** Across the holistic and detector passes the result reads as genuinely composed, not
remediated-template. The Print remediation (R0–R10, M1–M4) worked: the contact sheet is a real
instrument (dense table, dark detector windows, floating ink cull-bar, mono frame badges), the
serif/sans/mono voice discipline is *actually followed in the rendered hierarchy*, terracotta accent
is genuinely rationed, and the folio masonry varies card height by member count so it reads as a wall
of distinct figures rather than a card grid. Deterministic scan is nearly clean (4 hits total).
**B+ tier, distinctive, earns the right not to be called AI.**

The seams that remain are *consistency and encoding* problems inside a good system, not a generic
aesthetic. That is exactly the profile that the two follow-on projects are scoped to fix.

## Health scores

| Lens | Dimensions (0–4) | Headline |
|---|---|---|
| Design-system / primitives | token discipline **2** · primitive API **3** · reusability/coupling **2** · consistency **2** | System exists; adoption was never enforced |
| Plotting surfaces | encoding **2** · craft **2.5** · interaction grammar **1.5** · architecture **1.5** | Hue does everything; two duplicated plot cores; incoherent peak grammar |
| Cross-cutting technical | a11y **3** · perf **3** · responsive **2** · theming **3** · anti-patterns **2** | Solid modal a11y; two side-stripe bans; one unresponsive page |
| Holistic (Nielsen) | **29 / 40** | Real-interface band; pulled down by Consistency (#4) and error-help (#9) |

## The three systemic themes (these are the spine of the follow-on work)

### Theme A — "Define a token/role, then ignore it." (→ library extraction)
The design system is well-authored but **unenforced**, so the same concept renders N ways:
- `text-label` role exists, yet **37 hand-rolled kicker labels** reinvent it across 5 sizes (9/10/10.5/11/11.5px) and 4 tracking values.
- `ScoreBar` primitive exists with **zero importers**; the score bar is hand-inlined 3× at divergent dimensions (`PhasePanel.tsx:58`, `:153`, `SamplesPage.tsx:133`). `Input`, `Dot`, `SectionLabel` are also dead-exported.
- `.card` plate recipe exists, yet a plate is **re-derived inline 21×** with drifting radius and border weight; two builder pages inline the full Plate-Lift shadow.
- Radius scale is *stated* (sm 5 / md 7 / full) but **never tokenized** — only `--radius:6px` exists; `.card` hardcodes `12px`; 9 arbitrary `rounded-[Npx]` in use.
- **96 inline `text-[Npx]`** total (≈45 are redundant on-scale restatements; ≈50 are genuinely off-scale).
- **No lint guard** (only `tsc`/vite/vitest/playwright) and **no catalog/Storybook** — nothing prevents the next inline card or stops visual drift.

Promote-to-primitive candidates, with duplication already mapped:
| Candidate | Evidence | Priority |
|---|---|---|
| `SegmentedControl<T>` | `ScaleToggle`/`RepresentationToggle` byte-identical; `GroupingModeToggle` + `AnnotationToggles` use 2 *other* active vocabularies | **highest-value** |
| `PhaseChip` | inlined 4 ways with 10/13/20% tints (`SampleStatusChip:22`, `SpeculativeBuilder:243`, `MemberMetaRow:341`, `PhasePanel:141`) | high |
| `PhaseStrip` | `ScopingPhaseStrip` + `SeriesPhaseStrip` are near-duplicate files | high |
| `ModalShell` | `NavModal`/`ConflictModalShell`/`OnboardingFlow` duplicate backdrop+frame+trap; `ConflictModalShell` already proves the extraction | medium |
| `Kicker` | the 37 kicker labels | medium |
| `IconButton` | dismiss/× inlined 6 ways with drifting padding | low |
| `Card`/`Plate` | the 21 inline plates | medium |

Plus: **adopt-or-delete** the 4 dead primitives before shipping more (the worst failure mode for an
extraction project is shipping abstractions nobody adopts — which already happened once).

**What's already right (preserve as the model):** color/phase tokenization is exemplary — `phases.ts`
is the single source, every consumer goes through `phaseColor()`, AA-on-plate is test-pinned, 0 legacy
tokens remain. The 7 primitives that exist are correctly bounded (prop-driven, no store/query/queue
coupling). Make the type and radius scales behave the way color already does.

**Coupling reality (the extraction's real constraint):** 18/66 components import `useAppState`, 16
import TanStack Query, 2 import `lib/queue` — but that coupling lives in *feature* components. The
*visual* patterns worth extracting are already presentational or trivially separable. The "primitives +
catalog" scope is therefore achievable without touching the queue/SSE architecture.

### Theme B — "Hue does too much; shape does nothing." (→ plotting redesign, + a11y)
The Second-Channel rule (the product's HARD color-blind requirement) is honored on the small surfaces
and **dropped precisely where color is used most**:
- **Peaks**: auto / manual / excluded / optimistic are all the **same downward triangle**, separated only by hue + opacity + stroke (`TraceViewer.tsx:498-527`, `MemberTraceLayer.tsx:264-276`). The central curation distinction in the whole tool is invisible to a color-blind user. **[P1]**
- **Heatmap**: intensity is hue-only with **no colorbar, no legend, and per-row normalization that's never disclosed** (`MemberHeatmapLayer.tsx:138-165`). **[P1]**
- **Phase strips** (folio + scoping): per-segment phase by hue alone, no per-segment label/title/pattern (`SeriesPhaseStrip.tsx:34-43`, `ScopingPhaseStrip.tsx:29-37`). The most-repeated data element on the folio wall, unreadable in grayscale. **[P1]**
- **Detector rings**: uniformly gray (`ink-faint`) even though the trace peaks at the same q are phase-colored — the detector window throws away the phase signal where it reads most viscerally (`DetectorRingOverlay.tsx:108-112`). **[P1]**
- **Toast severity**: color-only `border-l-4`, no icon or word — error feedback, the highest-stakes transient signal, is hue-only (`ui/Toast.tsx:17-23`). **[P1]**

The fix is one principle applied everywhere: **shape/glyph/pattern carries category and state; color is
reserved for phase identity.**

### Theme C — Interaction grammar is implicit and overloaded. (→ plotting redesign)
- **Two conflicting "reset" meanings**: dblclick → full data range (`TraceViewer.tsx:327`), but the zoom indicator + Auto-fit → fit-to-peaks (`PlotCard.tsx:224,499`). Same word, opposite outcome. **[P1]**
- **"+ Peak" mode is theatre**: arming it changes nothing; empty-area click always adds regardless (`PlotCard.tsx:113-116`). **[P2]**
- **Single click is overloaded** into add / exclude / delete by 10px proximity + invisible peak source, with no undo (`TraceViewer.tsx:264-291`). **[P2]**
- **Two different peak grammars across surfaces**: click = curation in TraceViewer, click = display-cycle in MultiTracePlot (`:477-501`). Expert reloads a different mental model between Focus and Series. **[P2]**
- **Zoom/brush/reset undiscoverable**: raw listeners, no cursor change, no hint, indicator only appears after you've already zoomed. **[P2]**
- **Architecture**: TraceViewer (810) and MultiTracePlot (746) **duplicate ~400 LOC** of imperative Observable-Plot lifecycle (scale, wheel-zoom, dblclick, invert, resize, cleanup) with drift — any zoom fix lands twice. **[P2]**

Redesign point of view (from the plotting assessment): one `<PlotSurface>` primitive owning the Plot
instance + shared scale + all gestures, exposing `marks` / `overlay(scales)` / `hitTest()`; encode peak
class as shape and phase as color; an explicit **select-then-act** peak grammar (empty-click adds, peak-click
selects, explicit affordance removes) unified across both surfaces; one honestly-named "home" view; and
give the secondary panels their reasoning back (phase-color rings, label Miller axes, colorbar the heatmap).

## Cross-cutting fixes (independent of the two big projects — can land anytime)

- **[P1] Side-stripe accent border (absolute ban) ×2** — `InfrastructureBanner.tsx:55` and `ui/Toast.tsx:85` both render `rounded-md border-l-4` colored by status. The #1 AI tell, and `Toast.tsx:18` proudly comments it as a feature. Fix together: leading status icon + label + full hairline border; convey kind by icon, not edge hue.
- **[P1] Loupe page doesn't reflow** — `LoupePage.tsx:40,216` use `grid-cols-[1fr_286px]` with no breakpoint; the 286px sidebar squeezes the detector frame on narrow viewports. Copy `FocusWorkspaceLayout`'s `grid-cols-1 lg:grid-cols-[...]` cascade.
- **[P2] Focus-trap omits `<textarea>`** — `hooks/useFocusTrap.ts:3` `FOCUSABLE` selector misses `textarea`, so Tab-wrap is wrong at a textarea boundary (Notes drawer).
- **[P2] ResizeObserver → setState per frame ×3** — `TraceViewer.tsx:194`, `MillerPlot.tsx:50`, `MultiTracePlot.tsx:267` thrash full plot redraws on resize. Centralize a `useResizeKey(rAF)` hook.
- **[P2] Touch targets < 44px** on the compact segmented toggles (`px-1.5 py-0.5`, 10.5px text). (Folds naturally into the `SegmentedControl` primitive.)
- **[P2] `transition-all` / `transition-[max-width]`** — `ThumbnailGallery.tsx:67`, `SeriesBuilderPage.tsx:277` animate layout properties; scope the transition to box-shadow/opacity/transform.
- **[P2] No undo/confirm on the two destructive acts** — batch frame-drop (`CullBar.tsx:34-43`) and phase-call overrule fire optimistically with no toast/undo, while scoping has full undo+⌘Z+confirm. Apply scoping's model to the higher-stakes paths. (The queue already mints reverse-op ids.)
- **[P3] Em-dash claim is mostly a false positive** — of 411 repo occurrences, after stripping comments only ~16 reach output and all but one are the `?? "—"` null placeholder (legitimate). The single real prose em-dash is `InfrastructureBanner.tsx:61` "Couldn't save — try refreshing"; fix that one only. **Do not mass-replace placeholders.**
- **[P3] Duplicate accent token names** — `--color-accent` and `--color-print-accent` are the same value (drift risk); alias one to the other.
- **[P3] Decorative `backdrop-blur-sm`** on the "Indexing" badge (`ThumbnailGallery.tsx:90`) — drop it; the opaque pill suffices. (Modal scrims are fine.)

## Holistic UX notes (inform the redesigns, mostly not blocking)

- **Triage/Samples vocabulary mismatch** — onboarding teaches "Triage" (`OnboardingFlow.tsx:24,29`); the tab says "Samples" (`CorpusTopbar.tsx:20`). Pick one ("Triage" is the better domain word).
- **The header molecule is over-systematized** — terracotta kicker + serif numeral + uppercase micro-label is identical across Samples/Folio/Scoping and is the one place that brushes the dashboard-metric anti-reference. Vary it per surface.
- **The builder figure competes with its chrome** — the publication figure (the journey's emotional peak) shares its plate with a 5-pill kicker + Edit row and a 6-section always-open rail, with Export buried as rail-section-6. Default to full-bleed once a figure has ≥2 members; lift Export to a primary plate affordance.
- **Folio filter/sort chrome outweighs content** — 8 controls gate a small client-held listing; collapse behind a disclosure until ~12+ series.
- **Best-composed surface = Focus workspace.** The q-link triple (trace peak ↔ detector ring ↔ reflection row via one `hoveredQ` channel) is the most "real instrument" thing in the product. Preserve it; the plotting redesign should *enrich* what it carries (phase-color the rings), not regress it.

## How this maps to the program

1. **Library extraction (next)** — Theme A is the brief. Adopt-or-delete dead primitives, then promote
   `SegmentedControl` → `PhaseChip` → `PhaseStrip` → `Kicker`/`Card`/`ModalShell`/`IconButton`, tokenize
   the radius + type scales, add a lint guard and a catalog/gallery page. The `PhaseStrip` primitive is
   where the Theme-B second-channel fix lands for the strips. Folds in the touch-target and undo-affordance
   cross-cutting items.
2. **Plotting redesign (after)** — Themes B + C. One `<PlotSurface>` core, shape-carries-category peak
   encoding, select-then-act grammar unified across Focus + Series, named "home" view, phase-colored rings,
   labeled Miller axes, colorbar'd grayscale-safe heatmap. This is also where PR #278's deferred
   `MultiTracePlot` rename + cross-trace peak grammar resolve.
3. **Standalone cross-cutting fixes** — the side-stripe bans, loupe responsive, focus-trap, ResizeObserver,
   token alias, and the one real em-dash can be a quick independent sweep at any time (they don't depend on
   either big project).

**What to preserve through both rewrites:** the `phases.ts`/`coloring.ts` palette engineering, the
voice discipline, the `formatAxis` cross-surface formatter, the PlotCard auto-fit heuristic, the q-link
cross-highlight channel, and the modal/drawer focus-trap discipline.
