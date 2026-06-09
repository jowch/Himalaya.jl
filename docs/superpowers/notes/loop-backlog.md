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

## Score trend (per surface — Wave 1 baseline 2026-06-09)

| Surface | critique /40 | audit /20 | trend | target | gap |
|---|---|---|---|---|---|
| SamplesPage (contact sheet) | 31 | 15 | 31 | ≥36 / ≥18 | +5 / +3 |
| SeriesFolioPage | 30 | 17 | 30 | ≥36 / ≥18 | +6 / +1 |
| SeriesScopingPage | 28 | 16 | 28 | ≥36 / ≥18 | +8 / +2 |
| SeriesBuilderPage | 27 | 16 | 27 | ≥36 / ≥18 | +9 / +2 |
| FocusPage | 24 | 16 | 24 | ≥36 / ≥18 | +12 / +2 |
| LoupePage | 30 | 15 | 30 | ≥36 / ≥18 | +6 / +3 |

All six below GOLD on both axes. No P0s anywhere. Baseline reads: theming + anti-patterns are 4/4 almost everywhere (design-guard working); **a11y is the universal drag** (a11y dim 2–3 on every surface). Snapshots persisted under `.impeccable/critique/` (impeccable-native trend; re-score in Wave 4 to show the rise). Slugs: SamplesPage=`imalayaui-…-samplespage-tsx`, Folio=`ayaui-…-seriesfoliopage-tsx`, Scoping=`aui-…-seriesscopingpage-tsx`, Builder=`aui-…-seriesbuilderpage-tsx`, Focus=`himalayaui-…-focuspage-tsx`, Loupe=`himalayaui-…-loupepage-tsx`.

**Note on faithfulness:** Wave-1 critiques ran a single capable evaluation agent per surface applying the real `critique.md`+`audit.md` rubrics over source + rendered screenshots (desktop 1440 + laptop 1120 for the rail/grid surfaces) + the SAXS `ignore.md`, rather than the two-fully-isolated-subagent design/detect split. The deterministic detector (`impeccable detect --json`) was run centrally on every page source and returned `[]` (zero anti-pattern tells) — that is the Assessment-B half, and it corroborates the 4/4 anti-patterns scores. Re-scores in Wave 4 use the same method for trend comparability.

## Wave status

- [x] **Wave 0 — Foundation refresh** — DONE (`c525c93`). DESIGN.md re-rooted to `src/print/ui`, stale two-theme/grain promises killed, four new sections (State taxonomy / Spacing & density / Motion / Copy & UX writing); PRODUCT.md a11y de-dark-themed + named persona seeded.
- [x] **Wave 1 — Scored baseline** — DONE (2026-06-09). SAXS `ignore.md` seeded; 6 surfaces rendered + scored (critique /40 + audit /20); snapshots persisted; new P1s merged below. No P0s.
- [ ] **Wave 2 — P1 fixes top-down** ← IN PROGRESS. ✅ F-LIVE (`c07cdce`). ✅ F-ERRSILENT (closed-by-architecture). ✅ F-CONTRAST (`275c421`). Next cross-cutting P1s: the **rail adapt** (BU-RAIL+FO-RAIL+FO-COMB — layout, render-verify at 1120/1440), **SC-HEAD** (region headings), **F-A11Y** (operationalize contrast pairs as tests — now that ink-soft/ink-faint roles are settled). Then surface-specific P1s.

> **✅ HARNESS STATUS (rebuilt 2026-06-09):** the dev backend died mid-session and its `/tmp/himalaya-uat/` DB copy was reclaimed; **now restored** — sourced a surviving persistent DB at `/Users/me/projects/himalaya-devdata/himalaya.db` (139 samples / 1 series / 3 experiments), `sqlite3 .backup`'d it to a disposable `/tmp/himalaya-uat/himalaya.db`, and restarted `env HIMALAYA_DB_PATH=/tmp/himalaya-uat/himalaya.db julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- serve --port 8091` (from-source). `/api/samples` → HTTP 200. vite :5182 (VITE_API_PORT=8091) proxies to it; prod :8080 untouched. Render-verify is available. The canonical devdata DB stays pristine (we serve the copy). See [[feedback_live_audit_harness]].
- [ ] **Wave 3 — P2/P3 enhancement**
- [ ] **Wave 4 — Converge (`polish`) + re-score to GOLD**

## Backlog

Severity: P0 blocking · P1 major/any-WCAG-AA · P2 minor · P3 polish. Status: todo / in-progress / done / blocked.
`✓W1` = confirmed by the Wave-1 baseline. Prefer items that fix multiple surfaces at once.

### Cross-cutting (the system / multiple surfaces)

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| F-DOC | DESIGN.md/PRODUCT.md doc-debt | document+teach | P1 | **done** | `c525c93` Wave 0 |
| MO-DOC | motion vocabulary undocumented as a system | document+animate | P3 | **done** | `c525c93` Wave 0 Motion section |
| F-LIVE | **No `aria-live`/status-message on mutations & destructive keys** — Samples X-drop, Focus peak/phase/custom-index, Loupe X/R/flip, Builder confirm-progress. One polite+assertive live-region pattern fixes status feedback across 4 surfaces. | harden | P1 | **done** | `c07cdce` Wave 2 · new `lib/announce.ts`+`LiveRegion` SR-only announcer + assertive error toasts; visible-toast/SR-only split; 84 vitest + e2e 35 + build green |
| F-ERRSILENT | **Failed mutations fail silently** — Focus (`queries.ts:888` error read by nobody), Loupe, Builder; legacy of the cancelled conflict-UI decision. Surface plain-language error + retry. | harden | P1 | **closed (resolved-by-architecture)** | **Verify-before-review (Wave 2) found the baseline mis-cited.** Live surfaces mutate via `useQueueMutation`: 4xx→assertive validation toast w/ per-kind copy (`errors.ts:buildValidationMessage`, assertive via `c07cdce`); 5xx/network→auto-retry + `InfrastructureBanner` ("Saving…" >500ms → "Error. Couldn't save." + **Refresh** button >30s). Both asks (surface + retry) satisfied. Cited `queries.ts:888` "error read by nobody" = **retired comparison code** (pin/unpin + saveComparison, Compare retired #177) = dead/out-of-scope. Evidence-based (onError handler + InfrastructureBanner + queue contract tests); forced-error render-verify deferred as contrived. Linked P2s SA-RETRY/LO-ERR/BU-PROGRESS likewise mostly covered — re-triage in Wave 3. |
| F-CONTRAST | **`ink-faint` muted text fails AA (~3.16:1 on paper / 2.92 on paper-sunk)** — used for column headers, sample IDs, "Not indexed", kb legends across surfaces. Darken token to ≥L0.58 or restrict to large/decorative. | colorize | P1 | **done** | **Measured first** (oklch→linear-sRGB→luminance): ink-faint L0.640=**3.16:1** (AA-large only); ink-soft L0.467=6.50:1; ink=14.39. **Owner chose** (AskUserQuestion): keep ink-faint as tuned, reassign small-normal informational TEXT → `text-ink-soft`; leave kickers (AA-large), aria-hidden glyphs, icon controls (non-text 3:1), swatch fills, data-viz `var(--color-ink-faint)`. **50 files, ~72 swaps** + WCAG-context doc comments in styles.css. frontend-reviewer APPROVE (0 misclassifications); render-verified on Focus (captions/legends readable, kickers preserved). Commit `275c421` (gate: 2036 vitest + e2e 35 + build all green). |
| F-A11Y | a11y not operationalized: pin ink/accent/status contrast pairs as tests (like phase-on-plate); per-widget keyboard maps/aria | harden | P1 | todo | corroborated W1; see F-CONTRAST |
| F-STATE | state taxonomy: doc-half done (`c525c93`); enforcement half (rest/hover/active/focus-visible/disabled/busy/selected/error/read-only consistently applied) | harden | P1 | todo | |
| SC-HEAD | **Region labels are non-heading `div`s** (Kicker default) on Scoping (and Folio) — only h1 in the heading tree; pass `as="h2"/"h3"`. | harden | P1 | todo | ✓W1 · WCAG 1.3.1 · spans scoping+folio |
| EX-EMPTY | unify bespoke error/not-found divs (Samples, Focus, Loupe) onto EmptyState + action/retry slot | extract | P2 | todo | |
| EX-SPACING | spacing/density scale underspecified; promote shared-grid/alignment constants to tokens | extract | P2 | todo | |
| AN-MOTION | app essentially static: selection bars snap, drag drop-edge hard hairline, candidate-dim instant, no toast/modal motion | animate | P2 | todo | reduced-motion preserved |
| TY-DISPLAY | restore display-xl (31px) step (titles app-wide quieter than designed) | typeset | P3 | todo | reconcile w/ DESIGN.md display token |

### SamplesPage (contact sheet)

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| SA-SEM | **Non-semantic data table** — grid divs, no `role=table/row/columnheader/cell`, span headers, unlabeled thumbnail `<button>`s | harden | P1 | todo | ✓W1 · WCAG 1.3.1/4.1.2 |
| SA-RESP | **Unbuilt responsive scroll** — row min-width ~1018–1054px overflows near 1024px; sticky-Sample `overflow-x` container the comments describe doesn't exist | adapt | P1 | todo | ✓W1 · WCAG 1.4.10 |
| SA-KEYD | window keydown ("X" drops frames) guards only INPUT/TEXTAREA; no isContentEditable/open-modal check; no keyboard loupe-open (Enter) | harden | P1 | todo | ✓W1 · announce half → F-LIVE |
| SA-RETRY | error div has no retry CTA | clarify | P2 | todo | → EX-EMPTY/F-ERRSILENT |
| DI-SAHEAD | head body is a long single tagging-philosophy paragraph | distill | P3 | todo | |
| ON-EMPTY | empty corpus dead-ends ("No samples yet") with no add-data/learn path | onboard | P3 | todo | |
| LA-COLLIDE | CullBar/ComposeBar can collide with bottom InfrastructureBanner | layout+animate | P2 | todo | |

### SeriesFolioPage

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| FOL-KBD | **Series cards not keyboard-operable** — `onClick` on `<article>`, no role/tabIndex/onKeyDown; primary navigation mouse-only | harden | P1 | todo | ✓W1 · WCAG 2.1.1 |
| FOL-SORT | visible "SORT" label missing left of segmented control (mockup series-folio.html:326) | typeset | P2 | todo | |
| FOL-ERR | EmptyState used but no action passed on error | extract | P2 | todo | → EX-EMPTY |
| FOL-MISC | search input unlabeled (P2, WCAG 3.3.2) + no on-page new-series CTA + indistinguishable empty drafts | harden+craft+clarify | P2 | todo | ✓W1 |

### SeriesScopingPage

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| SC-ORDER | **"Ordered by" renders as an editable input but is inert** (static Field, no options/onClick) — wire data-driven dropdown + custom + restore of-note copy | clarify | P1 | todo | ✓W1 · ordering_vars=data-driven+custom; verify manifest/adapter |
| SC-KBD | **Drag-reorder has no keyboard alternative** (GripHandle aria-hidden; page never adds the path its primitive documents) | harden | P1 | todo | ✓W1 · WCAG 2.1.1 |
| SC-FOLD | **Candidate "also found" list inverts the fold** over the 2-member series + build action | layout | P1 | todo | ✓W1 |
| SC-MISC | candidate rows dead-end (no inline tag-add) + Discard no focus-visible ring + candidate sparkline aria-hidden phase | craft+harden | P2 | todo | ✓W1 |
| DI-SCOPECOLD | cold-panel stacked explanatory prose competes for one eye position | distill | P3 | todo | |

### SeriesBuilderPage

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| BU-DEAD | **Two inert controls in COMPOSE rail** — outline "+ Add sample" Button + collapse "›" chevron render live but get no handler; controls-don't-lie violation; duplicates working native select | harden | P1 | todo | ✓W1 |
| BU-RAIL | **Work·rail split gated at `xl:` (1280px)** drops rail below plate across 1024–1279 laptop band (Confirm/ordering/offset/export off-screen) | adapt | P1 | todo | ✓W1 · fix WITH FO-RAIL |
| BU-EXPORT | surface Export-figure to topbar (mockup series-builder.html:378) — builder exposes only rail "Copy as PNG" | clarify+layout | P1 | todo | loop-config "Export now"; Duplicate DEFERRED |
| BU-PROGRESS | confirm chain no visible progress + generic causeless error copy | clarify | P2 | todo | → F-LIVE/F-ERRSILENT |

### FocusPage

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| FO-COMB | **Comb/reflections panel `hidden lg:flex` drops at <1024px laptop band** (stack under detector, don't hide a primary analysis tool) | adapt | P1 | todo | ✓W1 · distinct from FO-RAIL |
| FO-RAIL | work·rail split gated at `xl:` — confirm whether the assignment rail (not just comb) also drops; fix WITH BU-RAIL | adapt | P1 | todo | verify vs FO-COMB |
| FO-EDIT | trace-edit model undiscoverable (arm "+Peak"/click=add/click-peak=remove/alt-click=exclude) + no keyboard accelerators / plot not keyboard-operable + generic detector canvas label | onboard+harden | P2 | todo | ✓W1 · WCAG 2.1.1/1.1.1 |
| FO-ERR | bespoke "Sample not found" instead of EmptyState | clarify+extract | P2 | todo | → EX-EMPTY |
| FO-DIM | candidate-hover losing-peak dim instant; drag drop-edge hard hairline | animate | P2 | todo | |
| DI-FOCUSNOTE | candidate rail hint is a dense 2-sentence indexing-theory run-on | distill | P3 | todo | domain-guard: don't dumb down |

### LoupePage

| id | opportunity | command | sev | status | notes |
|---|---|---|---|---|---|
| LO-KEYD | **Window keydown (X/R/arrows/Esc) guard incomplete** (no isContentEditable/open-popover) + no aria-live announcement | harden | P1 | todo | ✓W1 · WCAG 2.1.1/4.1.3; announce → F-LIVE |
| LO-FOCUS | "Back to the sheet" bare buttons rely on UA-default focus, not a focus-visible ring | harden | P2 | todo | ✓W1 · WCAG 2.4.7 (was P1; localized) |
| LO-THUMB | filmstrip thumbnails have no accessible frame identity (every one reads "Detector image, button"); add aria-label + aria-current/pressed | harden | P2 | todo | ✓W1 · WCAG 4.1.2 |
| LO-ERR | no failure feedback on drop/representative/tag mutations | harden | P2 | todo | ✓W1 · → F-ERRSILENT |
| LO-GRID | grid `[minmax(0,1fr)_286px]` never stacks — **honest at 1024px+ (in scope), rigid below** | adapt | P2 | todo | ✓W1 DOWNGRADED from P1 (in-scope OK) |

## Iteration log

- _2026-06-09 — loop designed, configured (4 decisions), tracker seeded (`docs` commit)._
- _2026-06-09 — **Wave 0 done** (`c525c93`): foundation doc refresh (PRODUCT.md/DESIGN.md). Closed F-DOC, MO-DOC; F-STATE doc-half done._
- _2026-06-09 — **Wave 2 / F-CONTRAST done** (`275c421`): measured the real contrast first (oklch→linear-sRGB→luminance: ink-faint=3.16:1 AA-large-only, ink-soft=6.50:1). Surfaced the fix-approach as a genuine design decision (AskUserQuestion) — owner chose "keep the token, reassign small-text uses to ink-soft" over darkening. Implementer audited all ~72 small-text usages across 50 files, left kickers/glyphs/icon-controls/data-viz faint, documented the WCAG-context rule in styles.css. frontend-reviewer APPROVE (0 misclassifications, verified by me too); render-verified on Focus; full gate green (2036 vitest, e2e 35, build). The "~6-10 files" estimate was wrong — there were genuinely ~50; the approved principle required all of them. Next: rail-adapt (layout, backend warm)._
- _2026-06-09 — **Wave 2 / harness rebuilt + F-ERRSILENT closed**: dev backend had died (its /tmp DB copy reclaimed); restored it from the surviving persistent `~/projects/himalaya-devdata/himalaya.db` (.backup'd to a disposable /tmp copy, served from-source on :8091, /api 200). Then verify-before-review on the error paths found **F-ERRSILENT already resolved** by `useQueueMutation` (4xx assertive toast + 5xx/network auto-retry + InfrastructureBanner Refresh) — the baseline's `queries.ts:888` cite was retired comparison code. Closed as resolved-by-architecture (no code needed; the honest outcome of checking live source). Next: F-CONTRAST / rail-adapt (backend now warm for render-verify)._
- _2026-06-09 — **Wave 2 / F-LIVE done** (`c07cdce`): cross-cutting aria-live status layer via the subagent cadence (fresh implementer → spec-review SPEC-COMPLIANT + frontend-reviewer APPROVE → 5 follow-ups incl. a real per-region-flip re-speak bug → verified). New `lib/announce.ts`+`LiveRegion` SR-only announcer + assertive error toasts; reused the already-mounted `showToast`/`ToastContainer` (verify-before-review corrected the baseline's "no global toast" claim). Gate: lint:design + tsc + 84 vitest + e2e 35 + build all green. Live render-verify blocked by a dead dev backend (DB copy reclaimed) — waived-with-note for this invisible change; harness must be rebuilt next wake. Partly advanced F-ERRSILENT (queue onError now assertive). Next: F-ERRSILENT remainder / F-CONTRAST / rail-adapt._
- _2026-06-09 — **Wave 1 done** (scored baseline): seeded SAXS `.impeccable/critique/ignore.md`; rendered all 6 surfaces (1440 + 1120 for rail/grid surfaces) via Playwright MCP against the dev-DB server (:5182, prod :8080 untouched); ran `impeccable detect` (all `[]`) + a per-surface critique/audit. Scores: Samples 31/15 · Folio 30/17 · Scoping 28/16 · Builder 27/16 · Focus 24/16 · Loupe 30/15. No P0s. Persisted 6 snapshots. Merged ~10 new P1s (F-LIVE, F-ERRSILENT, F-CONTRAST, SC-HEAD, SA-SEM, SA-RESP, FOL-KBD, SC-FOLD, BU-DEAD, FO-COMB) + confirmed SC-ORDER/SC-KBD/BU-RAIL/LO-KEYD/LO-FOCUS; downgraded LO-GRID→P2. **Next: Wave 2** — open with cross-cutting P1s (F-LIVE → F-ERRSILENT → rail adapt → F-CONTRAST → SC-HEAD)._
