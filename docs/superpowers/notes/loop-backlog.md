# Production-polish loop — backlog & tracker

Durable state for the impeccable production-polish loop. Companion to [2026-06-09-impeccable-design-loop.md](2026-06-09-impeccable-design-loop.md) (the command catalog, law checklist, severity machinery, per-surface map). Each wake re-reads both. One iteration per wake.

## loop-config (resolved with Jonathan 2026-06-09 — binding)

| Key | Value |
|---|---|
| **responsive_scope** | Desktop-first, **honest to ~1024px**. Fix `xl:`→`lg:` rail (Focus+Builder) + stack Loupe grid below a documented min-width. Sub-tablet (<~768px) out of scope; `<600px` scoping foot-row wrap item **dropped**. |
| **duplicate_export** | **Export-figure surfaced now** (P1, presentation/placement). **Duplicate-series deferred** out of this loop (separate `shape`-gated feature; no backend duplicate path today). |
| **finish_bar** | **GOLD / flagship.** Every surface **critique ≥ 36/40 AND audit ≥ 18/20** AND all P0/P1 closed. Genuinely-infeasible flagship items → mark `blocked`, surface, do not force. |
| **ordering_vars** | Scoping "Ordered by" dropdown is **data-driven** from the real manifest/adapter variable set + a "define your own" affordance. Verify against backend before wiring; do not hardcode the mockup's time/dose/temperature set. |
| **branch** | `worktree-greenfield-ui-rebuild` — STAYS UNMERGED. Never merge, never `finishing-a-development-branch`, never `git add -A`, every commit ends with the mandated co-author line. |

## Score trend (per surface — filled by Wave 1 + re-scores)

| Surface | critique /40 | audit /20 | trend | target |
|---|---|---|---|---|
| SamplesPage (contact sheet) | — | — | | ≥36 / ≥18 |
| SeriesFolioPage | — | — | | ≥36 / ≥18 |
| SeriesScopingPage | — | — | | ≥36 / ≥18 |
| SeriesBuilderPage | — | — | | ≥36 / ≥18 |
| FocusPage | — | — | | ≥36 / ≥18 |
| LoupePage | — | — | | ≥36 / ≥18 |

## Wave status

- [x] **Wave 0 — Foundation refresh** — DONE (`c525c93`). DESIGN.md re-rooted to `src/print/ui`, stale two-theme/grain promises killed, status note reframed as living spec, four new sections added (State taxonomy / Spacing & density / Motion / Copy & UX writing) grounded in scanned source; PRODUCT.md a11y de-dark-themed + reduced-motion reframed + named persona seeded. Display token left at 27px (31px = later render-verified task). Doc-only; gate unaffected.
- [ ] **Wave 1 — Scored baseline** (read-only; seeds rows). Per surface: `audit` + `critique`, with `ignore.md` seeded for SAXS vocabulary + a SAXS persona. Record /40, /20, slug; merge P0/P1 findings below.
- [ ] **Wave 2 — P1 fixes top-down**
- [ ] **Wave 3 — P2/P3 enhancement**
- [ ] **Wave 4 — Converge (`polish`) + re-score to GOLD**

## Backlog (seeded from synthesis §4 + decisions; Wave 1 will append audit/critique findings)

Severity: P0 blocking · P1 major/any-WCAG-AA · P2 minor · P3 polish. Status: todo / in-progress / done / blocked.

| id | surface | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|---|
| F-DOC | cross-cut | DESIGN.md/PRODUCT.md doc-debt: re-root to `src/print/`, kill stale two-theme/grain promises, add state-taxonomy/spacing/motion/copy | document+teach | P1 | **done** | `c525c93` Wave 0 |
| F-STATE | cross-cut | No documented state taxonomy (rest/hover/active/focus-visible/disabled/busy/selected/error/read-only; skeleton-vs-spinner-vs-empty-vs-inline-error) | document+harden | P1 | todo | doc-half DONE (`c525c93` State-taxonomy section); harden ENFORCEMENT half remains for Wave 2 |
| F-A11Y | cross-cut | a11y not operationalized: only phase-on-plate contrast pinned; pin ink/accent/status pairs; per-widget keyboard maps/aria | harden | P1 | todo | |
| SC-ORDER | scoping | "Ordered by" is immutable static text; wire data-driven dropdown + custom + restore of-note copy | clarify | P1 | todo | ordering_vars=data-driven+custom; verify manifest/adapter |
| SC-KBD | scoping | drag-reorder has no keyboard alternative (builder has ▲/▼; scoping doesn't) | harden | P1 | todo | |
| BU-EXPORT | builder | surface Export-figure to topbar (mockup series-builder.html:378) | clarify+layout | P1 | todo | Duplicate-series DEFERRED per decision |
| BU-RAIL | builder | work·rail split gated at `xl:` drops rail below work col on 1024–1279px | adapt | P1 | todo | fix WITH FO-RAIL (same change) |
| FO-RAIL | focus | work·rail split gated at `xl:` — same defect | adapt | P1 | todo | fix WITH BU-RAIL |
| LO-GRID | loupe | grid `[minmax(0,1fr)_286px]` has no responsive prefix; crushes detector at narrow widths | adapt | P1 | todo | stack below documented min-width |
| SA-KEYD | contact sheet | window keydown ("X" drops frames) guards only INPUT/TEXTAREA; no isContentEditable/open-modal check, no aria-live, no keyboard cull path | harden | P1 | todo | |
| LO-KEYD | loupe | window keydown (X/R/arrows/Esc) — same guard hole, no aria-live | harden | P1 | todo | shared pattern with SA-KEYD |
| EX-EMPTY | cross-cut | unify bespoke bordered error/not-found divs (Samples, Focus, Loupe) onto EmptyState + first-class action/retry slot | extract | P2 | todo | |
| FOL-ERR | folio | EmptyState used but no action passed on error | extract | P2 | todo | depends EX-EMPTY |
| TY-DISPLAY | cross-cut | restore display-xl step (titles app-wide quieter than designed) | typeset | P3 | todo | reconcile w/ DESIGN.md display token in Wave 0 |
| FOL-SORT | folio | visible "SORT" label missing left of segmented control (mockup series-folio.html:326) | typeset | P2 | todo | |
| LA-COLLIDE | contact sheet | CullBar/ComposeBar can collide with bottom InfrastructureBanner | layout+animate | P2 | todo | |
| EX-SPACING | cross-cut | spacing/density scale underspecified (sm/md/stage only); promote shared-grid/alignment constants to tokens | extract | P2 | todo | |
| AN-MOTION | cross-cut | app essentially static: selection bars snap, drag drop-edge hard hairline, candidate-dim instant, no toast/modal motion | animate | P2 | todo | reduced-motion preserved |
| SA-RETRY | contact sheet | error div has no retry CTA | clarify | P2 | todo | depends EX-EMPTY |
| FO-ERR | focus | bespoke "Sample not found" instead of EmptyState | clarify+extract | P2 | todo | depends EX-EMPTY |
| FO-DIM | focus | candidate-hover losing-peak dim instant; drag drop-edge hard hairline | animate | P2 | todo | |
| LO-FOCUS | loupe | "Back to the sheet"/"Sample not found" rely on UA-default focus, not a spelled-out focus-visible ring | harden | P2 | todo | |
| DI-FOCUSNOTE | focus | candidate rail hint is a dense 2-sentence indexing-theory run-on | distill | P3 | todo | domain-guard: don't dumb down |
| DI-SCOPECOLD | scoping | cold-panel stacked explanatory prose competes for one eye position | distill | P3 | todo | |
| DI-SAHEAD | contact sheet | head body is a long single tagging-philosophy paragraph | distill | P3 | todo | |
| ON-EMPTY | contact sheet | empty corpus dead-ends ("No samples yet") with no add-data/learn path | onboard | P3 | todo | |
| MO-DOC | cross-cut | motion vocabulary undocumented as a system | document+animate | P3 | **done** | `c525c93` Wave 0 Motion section |

## Iteration log

- _2026-06-09 — loop designed, configured (4 decisions), tracker seeded (`docs` commit)._
- _2026-06-09 — **Wave 0 done** (`c525c93`): foundation doc refresh (PRODUCT.md/DESIGN.md). Closed F-DOC, MO-DOC; F-STATE doc-half done. Fresh-implementer cadence + verify-by-diff. Next: Wave 1 scored baseline (6× audit + critique, seed SAXS ignore.md)._
