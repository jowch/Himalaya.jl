# Impeccable design loop — synthesis for "The Print" production-polish pass

_Written 2026-06-09. Branch: `worktree-greenfield-ui-rebuild` (STAYS UNMERGED). Scope: drive a multi-iteration production-polish pass over the six greenfield "Print" surfaces (and backend where it affects UX) using the impeccable skillset, held to a tracked, scored bar._

This is a planning/orientation doc. The paste-ready loop prompt lives in §7 and is also returned as `loopPromptDraft`.

---

## 0. Load-bearing environment facts (verified against live source)

These shape every decision below.

- **Context-doc location mismatch.** `PRODUCT.md` (4.6 KB) and `DESIGN.md` (19 KB) live at the **worktree root** (`./PRODUCT.md`, `./DESIGN.md`). The app, the design system, and the design guard live **two levels down** at `packages/HimalayaUI/frontend/`. impeccable's `load-context.mjs` resolves from cwd-root → `.agents/context/` → `docs/`. So **any impeccable command must be run with the context made resolvable** — either `cd` to the worktree root (where the docs are) or set `IMPECCABLE_CONTEXT_DIR=/Users/me/projects/Himalaya.jl/.claude/worktrees/greenfield-ui-rebuild` so the loaded PRODUCT/DESIGN are the real ones, not "no context found." This is the #1 setup gate.
- **No `.impeccable/` yet.** No live-mode config, no critique storage, no design.json sidecar. First evaluate run will create `.impeccable/critique/` for trend persistence. A `.impeccable/live/config.json` must be authored before any `live` session (the app is React/TSX → JSX variant-writing rules).
- **Register is unambiguously PRODUCT.** Authenticated dense scientific instrument (data tables, trace plots, detector windows, assignment rails). `product.md` laws govern, not `brand.md`. PRODUCT.md already declares `Register: product`.
- **The gate that must stay green** (run from `packages/HimalayaUI/frontend/`):
  `npm run lint:design` → `tsc --noEmit -p tsconfig.build.json` → targeted `npm test` → `npm run e2e` → `npm run build`. `lint:design` (check-design.mjs, exit 2 on any violation) is the build's first step and is non-negotiable.
- **Design guard bans** (outside exempt dirs `src/print/{ui,plot,detector,comb,export}/**`): inline `text-[…]`, `rounded-[…]`, raw color literals in color utilities, side-stripes (`border-l/r-N`), any raw color literal on a line. Color-authoring allowlist for rules #3/#5 only: `phases.ts`, `lib/comparison/coloring.ts`, `main.tsx`, `lib/figure-export/**`. **Any new appearance must land in a `src/print/ui` primitive**, consumers pass placement-only className.
- **Import boundary:** `src/print/**` may import NEW (`src/print/**`) + CARRIED (`queries/api/state/lib/hooks`) only. Never `src/components/**` or `src/pages/**` (deleted/relocated).
- **Six surfaces / 11 page files / 93 ui primitives.** Pages: `SamplesPage` (contact sheet), `SeriesFolioPage`, `SeriesScopingPage`, `SeriesBuilderPage`, `FocusPage`, `LoupePage` (+ co-located `*Adapters.ts`/`scopingDerive.ts`). Storybook wired (`npm run storybook`, ~per-primitive `.stories.tsx`). Mockups at `docs/redesign-mockups/*.html`.
- **Verify-by-rendering:** Playwright MCP against a writable DB **copy**, never prod `:8080`. Tests assert `data-*`/semantics, never Tailwind class strings.

---

## 1. Command catalog (organized by category, with relevance to The Print)

Relevance key: **★★★** drives the loop · **★★** situational · **★** low/guarded.

### Setup & Register

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `teach` | Writes strategic `PRODUCT.md` (register/users/purpose/voice/anti-refs/principles/a11y) via interview. | Only to **refresh** PRODUCT.md — reconcile the stale two-theme AA promise (Darkroom was deleted; single light identity), add the SAXS-scientist persona context. | Interview ≥1 round → write template → mandatory `load-context.mjs` re-run. Gate: never overwrite silently. | **★★** PRODUCT.md already strong; only a targeted refresh + persona seed. |
| `document` | Writes/refreshes `DESIGN.md` (six fixed sections + frontmatter) + `.impeccable/design.json` sidecar. Scan mode extracts real tokens. | Re-root DESIGN.md paths to `src/print/`, reconcile dark-theme/shadow drift, add the missing **state taxonomy / spacing-density / motion / copy** sections, frame as extendable not as-built. | Scan tokens → stage frontmatter → one creative-layer ask → write six sections → sidecar → `load-context.mjs` re-run. Gate: anti-refs must reappear as Don'ts; OKLCH accepted with linter warning (project is OKLCH-doctrine). | **★★★** The single biggest doc-debt fix; sharpens the identity lock every later command (and `live`) reads. |

### Build

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `craft` | End-to-end build: shape → (image-gen gated) → build to production bar → **iterate visually** → present. | Mostly **not** as full-build (system is shipped). Its **Step 5 (iterate visually) + production bar** is the model for every polish iteration; its **shape** front-gate scopes each task. | Gates: shape-confirm; on this no-image-gen harness, codex collapses — **the existing mockups ARE the approved north-star contract.** Exit bar: defensible in a studio review; live result has the mockup's major ingredients. | **★★★** as a *discipline* (Step 5 + production bar), **★** as a from-scratch builder. |
| `shape` | Design-planning-only brief via discovery interview; no code. | Front-door for any **net-new** affordance (e.g. surfacing Duplicate/Export, a first-run CTA) where scope/content isn't pinned by docs+mockup. | ≥1 user round → brief → STOP for confirm. Compact 3–5 bullet default. | **★★** for the handful of genuinely-new affordances; most polish is pinned by mockups and skips it. |
| `codex` | Image-gen visual-direction sub-flow (palette → mocks → approval). | **Skipped** — no native image_gen; the canonical `docs/redesign-mockups/*.html` already lock the visual direction. | n/a (its "live result must contain the approved mock's ingredients" rule survives as our fidelity bar). | **★** (mockups replace it). |
| `extract` | Lift recurring patterns into the design system + migrate + delete + document. 3+-same-intent threshold. | Promote the convention-only **shared-grid/alignment constants** and the underspecified **spacing/density scale** into real tokens; unify the **bespoke error/not-found divs** onto one `EmptyState` (+ `action` slot). | Discover system → 3+-same-intent gate → enrich (typed props, a11y) → migrate-then-delete → document. | **★★★** directly matches the closed-look/open-placement contract and the "EmptyState split" seam. |

### Evaluate

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `critique` | Dual-track (LLM design review + 27-pattern detector), synthesized to Nielsen /40 + AI-slop verdict + cognitive-load + persona red-flags + P0–P3, **persisted for trend**. | **Per-surface baseline** at loop start, and re-run after fixes to watch the /40 trend rise. | Two **isolated** sub-agents, each own browser tab → synthesize → persist via `critique-storage.mjs` → trend line. Gate: isolation non-negotiable; `ignore.md` drops domain-jargon false-positives silently. | **★★★** the loop's scored memory. **Caveat:** seed `ignore.md` + a SAXS persona so heuristic-2 (jargon) doesn't penalize Pn3m/Miller/q-value vocabulary. |
| `audit` | Code-level 5-dim scan (a11y/perf/responsive/theming/anti-patterns) → /20 + P0–P3. Documents, doesn't fix. | **Per-surface baseline** alongside critique — maps cleanly onto our hard invariants (theming→check-design, a11y→focus-trap/keyboard, responsive→the `xl:` rail + Loupe grid). | Per-dim 0–4 anchors → report by severity (location + WCAG + fix + suggested command). | **★★★** the technical half of the baseline; its findings seed the backlog directly. |

### Refine

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `polish` | Final pre-ship pass; **first act = align to design system** (classify drift: missing-token / one-off / conceptual). ~22-item checklist. | The **convergence step of every iteration** — each refine/enhance/fix hands off here. Optionally folds the persisted critique snapshot. | Gates: forbidden on non-complete work; MUST align to DS or ASK; detector output = defect evidence only. | **★★★** maps 1:1 onto closed-look/open-placement + "never port legacy presentation." |
| `harden` | Resilience: extreme inputs, long names, large datasets, concurrent ops, error states, a11y resilience, i18n. | The **window-keydown a11y holes** (X-drops, no `isContentEditable`/open-modal guard, no aria-live), **virtualize 1000+ member/trace lists**, double-submit/optimistic-rollback (mutation-queue), per-status error recovery, the undocumented state taxonomy. | Gate: server-side validation non-negotiable; one component erroring must not block the UI; verify against the edge-case checklist. | **★★★** highest-leverage technical pillar for a real-data scientific tool. |
| `distill` | Strip what doesn't earn its place; no nested cards; one spacing scale. **Guard: never oversimplify complex domains.** | The dense run-on rail microcopy (Focus candidate note, scoping cold-panel stacked prose). | Verify: faster task / lower load / still complete / clearer hierarchy. | **★★** useful but apply the domain guard — don't dumb down indexing data. |
| `bolder` | Amplify impact — **for Product = sharper hierarchy/weight-contrast/one-accent, NOT drama.** | Restore the **31px H1** the type scale truncates to 27px (titles app-wide quieter than designed). | Gate: "would they believe AI made this bolder?" → fail. | **★** narrow; theatrics undermine a data tool. |
| `quieter` | Reduce intensity so the tool disappears into the task. | Likely already at the right calm-instrument bar; spot-use if a surface reads loud. | Verify: still functional/distinctive/restrained-not-absent. | **★** the design already targets this. |
| `onboard` | First-value + the 5 empty-state types. | The **empty-corpus dead-end** ("No samples yet" with no path) despite an existing `OnboardingFlow`/`NavModal` in the shell. | Empty states teach; respect dismissals; LocalStorage seen-state. | **★★** real but bounded (P3). |

### Enhance

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `layout` | Space as the primary tool: squint test, rhythm, Flex-1D/Grid-2D, no nested cards, semantic z-index. | **Bottom-fixed collision** (CullBar/ComposeBar vs InfrastructureBanner), spacing-rhythm consistency, placement of newly-surfaced Duplicate/Export. | Squint test; 8–12px siblings / 48–96px sections; semantic z-index. | **★★★** highest enhance leverage for a dense instrument. |
| `animate` | Motion that conveys state; Product = 150–250ms, no page-load choreography, mandatory reduced-motion. | The app is **essentially static** — selection action bars snap, drag drop-edges are hard hairlines, candidate-dim is instant, no toast/modal motion. | Duration table; exponential easing; no bounce; reduced-motion preserves functional anims. | **★★★** clear, bounded wins; fits the 120ms house style. |
| `typeset` | Intentional type scale/hierarchy; tabular-nums for data; fixed rem for app UI. | The **missing 31px display step**, the **folio "SORT" label**, tabular-nums audit on q-values/scores. | Body ≥16px; line 45–75ch; max 2–3 families / 3–4 weights. | **★★★** several concrete, verified fidelity snaps. |
| `colorize` | Strategic palette; Product = Restrained, accent for action/selection/state only. | Mostly **already enforced** (phases live in primitives, no side-stripe, OKLCH tokens). Use as a WCAG/consistency audit, not a redesign. | 60/30/10; WCAG 4.5:1/3:1; no color-alone. | **★** consistency check only. |
| `delight` | Personality at earned moments only. | Narrow: phase-call confidence, reanalyze-complete, error-recovery. AI-slop loading-copy warning applies. | <1s; never blocks; reduced-motion. | **★** one or two earned moments at most. |
| `overdrive` | Push past conventional limits — one extraordinary moment. **Propose 2–3 → ASK before any code.** | Only a **functional-UI** moment if justified (virtual-scroll a huge trace list, View-Transitions dialog morph) — and only behind the propose-then-ask gate + browser verification. | Mandatory banner; propose→ask gate; 5 verify tests; reduced-motion fallback. | **★** highest risk on a scientific tool; default to skip. |

### Fix & Foundations

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `clarify` | Rewrite unclear UX copy; button/error/empty/loading/confirmation conventions. | The **scoping "Ordered by" immutable text** (mockup wants a changeable dropdown + full of-note), **Duplicate/Export labels**, unify error copy + reload CTA. | Verify: comprehension/actionability/brevity/consistency/tone. | **★★★** directly fixes verified fidelity gaps + copy seams. |
| `adapt` | Re-think experience for a new context (responsive). | The **`xl:` (1280px) rail breakpoint** dropping the rail below the work column on 1024–1279px laptops; the **fixed two-column Loupe grid** crushing the detector at narrow widths. | Content-driven breakpoints; input-method detection; test real widths. | **★★★** the verified P1 responsive defects. **Decision needed:** is small-viewport even in scope (see §8)? |
| `optimize` | Performance — measure first, fix the actual bottleneck, measure again. | Only on a **measured** bottleneck (d3 trace-plot/waterfall/heatmap render, long member lists). The doc forbids speculative optimization. | LCP<2.5s, INP<200ms, CLS<0.1; real-device. | **★★** gated on evidence; pair with `harden`'s virtualization. |

### Live & Tooling

| Command | What it does | When WE'd use it | Flow / gates / output | Relevance to The Print |
|---|---|---|---|---|
| `live` | HMR select-element → 3 variants → dial param knobs → accept/carbonize into source. | **Optional accelerator** for fine-tuning dense visual surfaces (spacing/density/type-scale knobs on plates/rails) once a DESIGN.md + `.impeccable/live/config.json` exist. **Default mode** (vary within identity) fits a shipped DS. | Identity-lock → default/departure mode gate → squint test → param budget. Caveat: React-rendered elements often hit wrap's `agent-driven` fallback. | **★★** powerful but setup-heavy; not required for the core loop (Playwright MCP already gives verify-by-rendering). |

---

## 2. Foundational design laws — checklist we hold every surface to

(Product register. These are the bar, not aspirations.)

**Color**
- [ ] OKLCH everywhere; chroma reduced near L 0/100. (Already true — neutrals tinted hue ~55–90.)
- [ ] No `#000`/`#fff` for large areas; every neutral tinted toward brand hue.
- [ ] Restrained strategy: accent (terracotta `0.555 0.150 38`) used for **primary action / current selection / state only**, never decoration. 60/30/10 visual weight.
- [ ] WCAG: body 4.5:1, large/UI 3:1. **Pin the currently-unpinned pairs** (ink/accent/status on paper+plate) as tests like `phases.test.ts` already pins phase-on-plate.
- [ ] No gray-on-color (use darker shade of the bg color or transparency).
- [ ] Color never the sole channel (Second-Channel Rule / Semantic-Dot) — deuteranopia + protanopia verified.

**Theme**
- [ ] Single light "Print" identity is the deliberate scene (publication plate on warm paper). **Reconcile PRODUCT.md's stale two-theme AA promise.**

**Typography**
- [ ] Fixed rem scale (not fluid clamp) — correct for a dense app. Tighter ratio.
- [ ] Body ≥16px equivalent; data/tables may run denser.
- [ ] Serif = titles, Plus Jakarta = chrome, mono = measurements (Serif-Means-Title / Monospace-Means-Measurement). tabular-nums on all aligned numerics (q-values, scores).
- [ ] Hierarchy via scale + weight; titles at intended size (**add the 31px display-xl step**).

**Layout**
- [ ] Predictable grids; consistency IS an affordance. Squint test passes on each surface.
- [ ] No nested cards; cards only for genuinely distinct/actionable content; spacing/dividers within.
- [ ] Responsive behavior is **structural** (collapse rail, stack panels) not fluid type. **Fix the `xl:` rail + fixed Loupe grid.**
- [ ] Semantic z-index scale; bottom-fixed surfaces don't collide.

**Motion**
- [ ] 150–250ms (house style 120ms ease-out); motion conveys **state only**; NO page-load choreography.
- [ ] No bounce/elastic. `prefers-reduced-motion` near-zeroes spatial motion but **preserves functional anims** (spinners/progress/focus).

**Absolute bans** (already enforced by guard + standing mandates)
- [ ] No side-stripe borders (>1px colored border-l/r accent).
- [ ] No gradient text, no glassmorphism-as-default, no hero-metric template, no identical card grids, no modal-as-first-thought.
- [ ] No em dashes (and not "--") in copy.

**AI-slop test (two altitudes)**
- [ ] First-order: theme+palette not guessable from "scientific tool" alone. (Pass — publication-plate identity is specific.)
- [ ] Second-order: aesthetic family not guessable from category+anti-refs. (Pass — sits in the deliberate gap between the four anti-refs.)
- [ ] Product slop test: would a user fluent in Linear/Figma/Notion trust this, or pause at every subtly-off component?

---

## 3. Scoring & severity machinery (how we track the loop)

**Per-surface scored baseline (run once at start, re-run after fixes):**
- `critique` → Nielsen **/40** + rating band (36–40 ship · 28–35 good · 20–27 acceptable · 12–19 poor · 0–11 critical), AI-slop verdict, cognitive-load failure count (0–1 low / 2–3 moderate / 4+ critical), persona red-flags.
- `audit` → **/20** across a11y/perf/responsive/theming/anti-patterns (18–20 excellent · 14–17 good · 10–13 acceptable · 6–9 poor · 0–5 critical).
- Persisted via `critique-storage.mjs` → **trend line per surface slug**. This is the closed-loop signal: re-run shows `24 → 28 → 32`.

**Universal severity ladder (tag every backlog item):**
- **P0 Blocking** — prevents task completion; fix immediately.
- **P1 Major** — significant difficulty/confusion, **or any WCAG AA violation**; fix before "done."
- **P2 Minor** — annoyance with a workaround; next pass.
- **P3 Polish** — nice-to-fix, no real user impact.
- Tie-breaker: "Would a user contact support about this?" → if yes, at least P1.

**Cognitive-load caps (LLM review applies):** ≤4 working-memory items, ≤4 visible options per decision point, ≤5 top-level nav, 1 primary + 1–2 secondary actions. Named violations to watch on dense surfaces: Wall-of-Options, Memory-Bridge, Visual-Noise-Floor, Jargon-Barrier (**domain-calibrated** — Pn3m/Miller/q are correct, not slop).

**Loop target (proposed, confirm with user):** every surface reaches **critique ≥ 32/40 AND audit ≥ 16/20**, and **all P0/P1 resolved**. (See §8 open decisions.)

---

## 4. Surface-by-surface improvement map

Each item: opportunity → impeccable command → severity. (P1s are the loop's first wave.)

### Contact sheet — `SamplesPage.tsx`
- Global `window` keydown ("X" drops selected frames) guards only INPUT/TEXTAREA — no `isContentEditable`, no open-modal/menu check, no aria-live on the bulk-reject result, no keyboard selection path for the cull → **`harden` · P1**.
- Bespoke bordered error div ("Could not load the corpus") instead of shared `EmptyState`; no retry CTA → **`clarify`/`extract` · P2**.
- Empty corpus dead-ends ("No samples in the corpus yet.") with no add-data/learn path → **`onboard` · P3**.
- CullBar snaps in/out at page root, can collide with bottom InfrastructureBanner → **`animate` + `layout` · P2**.
- Head body is a long single tagging-philosophy paragraph → **`distill` · P3**.

### Series folio — `SeriesFolioPage.tsx`
- Visible "SORT" label left of the segmented control is **missing** (only aria-label) — mockup `series-folio.html:326` → **`typeset` · P2**.
- FolioHeader title renders at 27px not the mockup's 31px → **`typeset` · P3** (shared fix).
- Error uses `EmptyState` (good) but folio passes no action → **`extract` (EmptyState action slot) · P2**.

### Series scoping — `SeriesScopingPage.tsx`
- **"Ordered by" renders as immutable static text** (`orderedBy={keyLabel}` with no `orderOptions`/`onOrderSelect`) — the single most important decision on the page; mockup specifies a changeable dropdown + full of-note ("Change it to time, dose, temperature, or define your own") → **`clarify` (wire control + restore copy) · P1**.
- Scoping drag-to-reorder (`useDragReorder`) has **no keyboard alternative at all** (builder has ▲/▼; scoping doesn't) → **`harden` · P1**.
- Cold-panel stacked explanatory prose (panel + 42ch note + button) competes for one eye position → **`distill` · P3**.
- Foot row (`justify-between` Dot+text+Button) wraps awkwardly <600px → **`adapt` · P2** (if narrow-viewport in scope).

### Series builder — `SeriesBuilderPage.tsx`
- Topbar **Duplicate-series (ghost) + Export-figure (solid)** absent; Copy-as-PNG buried in the rail, no Duplicate path at all — mockup `series-builder.html:378` → **`clarify` (then `layout` for placement) · P1**.
- Work·rail split gated at `xl:` (1280px) → rail drops below work column on 1024–1279px → **`adapt` · P1**.

### Focus — `FocusPage.tsx`
- Work·rail split gated at `xl:` (1280px) — same defect → **`adapt` · P1** (fix both pages together).
- Candidate rail hint is a dense two-sentence indexing-theory run-on ("Candidates that explain the same peaks swap; independent phases coexist") on the most technical surface → **`distill` · P3**.
- Bespoke bordered "Sample not found" instead of `EmptyState` → **`clarify`/`extract` · P2**.
- Candidate-hover losing-peak dim is instant; drag drop-edge is a hard hairline → **`animate` · P2**.

### Loupe — `LoupePage.tsx`
- Grid `grid-cols-[minmax(0,1fr)_286px]` has **no responsive prefix** — fixed two-column at every width, crushes the detector frame at narrow widths → **`adapt` · P1**.
- Global `window` keydown ("X" drop, "R" representative, arrows flip, Esc back) — same INPUT/TEXTAREA-only guard, no aria-live, no open-modal check → **`harden` · P1**.
- "Back to the sheet" / "Sample not found" are bare buttons/anchors relying on UA-default focus, not a spelled-out focus-visible ring → **`harden` · P2**.

### Cross-cutting (the system, not a page)
- **No documented state taxonomy** (rest/hover/active/focus-visible/disabled/busy/selected/error/read-only + skeleton-vs-spinner-vs-empty-vs-inline-error rules) → **`document` + `harden` · P1**.
- **DESIGN.md/PRODUCT.md doc-debt**: paths rooted at old `src/components/ui`, stale two-theme AA promise, Plate-Lift shadow drift, written as remediation log not extendable system → **`document` (+ `teach` refresh) · P1**.
- **Spacing/density scale underspecified** (3 named tokens; rhythm ad-hoc) → **`extract` · P2**.
- **A11y not operationalized** (only phase-on-plate contrast pinned; no per-widget keyboard maps/aria) → **`harden` · P1**.
- **31px display-xl step missing** (titles app-wide quieter than designed) → **`typeset` · P3** (one fix, app-wide effect).
- **Motion vocabulary undocumented as a system** → **`document`/`animate` · P3**.

---

## 5. Recommended command sequence (evaluate-first, polish-not-build)

This is an existing, disciplined system. **Do not rebuild.** Establish a scored baseline, fix top-down by severity, converge on `polish`, re-score.

**Wave 0 — Foundation refresh (unblocks everything; do FIRST).**
1. `document` (scan mode) — re-root DESIGN.md to `src/print/`, reconcile dark-theme/shadow drift, **add the state-taxonomy / spacing-density / motion / copy sections**, frame as extendable. Then a light `teach` refresh to fix PRODUCT.md's stale two-theme promise + seed the SAXS persona. Mandatory `load-context.mjs` re-run so the rest of the loop reads the fixed docs. _Rationale: every later command (and `live`) reads these; fixing them first makes all downstream evaluation/refinement correct-by-context._

**Wave 1 — Scored baseline (read-only; seeds the backlog).**
2. `audit` + `critique` per surface (6× each), with `IMPECCABLE_CONTEXT_DIR` set, an `ignore.md` seeded for SAXS domain vocabulary, and a SAXS persona. Persist all twelve snapshots → record the /40 and /20 baseline + trend slugs in the tracker. _Rationale: turns "is it better?" into tracked numbers; the audit/critique P0–P1 findings become the backlog, merged with the verified recon items in §4._

**Wave 2 — P1 fixes, top-down (the bulk of the loop).** For each P1, run the mapped command → `polish` handoff → verify-by-rendering → gate green → commit:
3. `clarify` — scoping "Ordered by" dropdown + copy; builder Duplicate/Export surfacing.
4. `adapt` — `xl:`→`lg:` rail (Focus + Builder together); responsive Loupe grid.
5. `harden` — window-keydown guards + aria-live (Samples, Loupe); scoping keyboard-reorder; a11y operationalization (pin contrast pairs, keyboard maps).
6. `extract` — unify error/not-found onto `EmptyState` (+ `action` slot); promote spacing/density + shared-grid constants into tokens.

**Wave 3 — P2/P3 enhancement.**
7. `typeset` — 31px display-xl step (app-wide), folio "SORT" label, tabular-nums audit.
8. `layout` — bottom-fixed collision (CullBar/ComposeBar vs InfrastructureBanner), Duplicate/Export placement, rhythm consistency.
9. `animate` — selection-bar enter/exit, drag drop-edge settle, candidate-dim, toast/modal (reduced-motion preserved).
10. `distill` — Focus candidate note + scoping cold-panel + Samples head copy.
11. `onboard` — empty-corpus first-run CTA.

**Wave 4 — Converge + re-score.**
12. `polish` final pass per surface (fold the persisted critique snapshot).
13. Re-run `critique` + `audit` per surface → confirm trend rose to target (≥32/40, ≥16/20) and all P0/P1 closed. Surfaces below target re-enter Wave 2/3.

`optimize` and `overdrive` are **on-demand only** — `optimize` solely on a measured bottleneck, `overdrive` only if the user explicitly wants one extraordinary functional moment (propose-then-ask gate).

---

## 6. How `live` fits our verify-by-rendering discipline

We already have a strong verify-by-rendering path (Playwright MCP against a writable DB copy), which the loop uses by default. `live` is an **optional accelerator**, not a requirement:

- **Where it helps:** dialing dense visual params without regeneration — density/spacing/type-scale **param knobs** on plates and rails (`layout`/`typeset`/`colorize` actions in default mode, which is "vary within the existing identity" — exactly right for a shipped DS). The identity-lock latches onto the existing `src/print/ui` tokens even with a thin DESIGN.md (and Wave 0 makes the lock sharper).
- **Setup gates before it works:** (a) the context-doc location must be resolvable (`IMPECCABLE_CONTEXT_DIR` or run from worktree root); (b) author `.impeccable/live/config.json` pointing at the Vite entry, with `commentSyntax: jsx`; (c) variant-writing uses TSX rules (template-literal `<style>`, className/style props).
- **The friction for THIS app:** the six surfaces render elements at runtime from React components/data (trace plot, detector image, gallery rows). `live-wrap.mjs` will frequently return the `agent-driven` fallback (`element_not_found` / `element_not_in_source`) rather than the deterministic path — expect to drive the fallback (find the real component source, scaffold a temporary served-file wrapper, persist on accept). This means `live` is best reserved for **leaf primitives in Storybook** (where wrap is deterministic) rather than live data pages.
- **Carbonize lands in true source** through the real edit path — compatible with the gate, but the accepted CSS must move into a `src/print/ui` primitive (not inline) to satisfy the design guard, and the resulting change still goes through the normal subagent review + gate + commit cadence. **`live` never bypasses the gate or the commit constraints.**

Net: keep Playwright MCP as the loop's primary verification; treat `live` as an opt-in tool for primitive-level tuning runs, gated on the setup above.

---

## 7. The loop prompt

Seeded from this doc's §4 backlog + the Wave-1 scored baseline, it drives one opportunity per iteration through the subagent cadence, verifies by rendering, keeps the gate green, commits per the hard constraints, and self-paces via ScheduleWakeup — stopping for genuine product decisions and at the §3 finish line. Paste-ready:

```
PRODUCTION-POLISH LOOP — "The Print" greenfield UI (branch worktree-greenfield-ui-rebuild). Run one iteration per wake. Read docs/superpowers/notes/2026-06-09-impeccable-design-loop.md (the synthesis doc) before acting — it is your command catalog, law checklist, severity machinery, and per-surface backlog.

ABSOLUTE CONSTRAINTS (violating any = stop and ask):
- NEVER merge this branch. NEVER run finishing-a-development-branch. The branch STAYS UNMERGED.
- NEVER `git add -A` / `git add .`. Stage only specifically-named files by exact path. Never stage *.bones.json, registry.ts, or stray .js.
- Every commit's exact last line MUST be: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- All appearance lives in src/print/ui primitives; consumers pass placement-only className. New appearance utilities outside ui/** fail the guard — put them in a primitive.
- src/print/** imports NEW (src/print/**) + CARRIED (queries/api/state/lib/hooks) ONLY. Never src/components/** or src/pages/**.
- Tests assert data-*/rendered semantics, never Tailwind class strings.

SETUP (every wake, from repo root /Users/me/projects/Himalaya.jl/.claude/worktrees/greenfield-ui-rebuild):
- export IMPECCABLE_CONTEXT_DIR to the worktree root so impeccable loads the real PRODUCT.md/DESIGN.md (they sit at root; the app is at packages/HimalayaUI/frontend).
- All npm/gate commands run from packages/HimalayaUI/frontend.
- Maintain the tracker at docs/superpowers/notes/loop-backlog.md (create on first wake from the §4 backlog in the synthesis doc + the Wave-1 baseline). Each row: id | surface | opportunity | impeccable command | severity (P0–P3) | status (todo/in-progress/done/blocked) | critique/audit slug + last score.

WAKE ALGORITHM:
1. If loop-backlog.md does not exist OR Wave-0/Wave-1 not done: do the next un-done foundation/baseline step, commit, update tracker, ScheduleWakeup, end.
   - Wave 0: run impeccable document (scan mode) to re-root DESIGN.md to src/print/, reconcile the stale two-theme AA promise + Plate-Lift shadow drift, and add State-Taxonomy / Spacing-Density / Motion / Copy sections; then a light impeccable teach refresh of PRODUCT.md (single light identity; seed a SAXS-structural-scientist persona). Commit doc changes (PRODUCT.md, DESIGN.md, .impeccable/design.json by exact path).
   - Wave 1: per surface (SamplesPage, SeriesFolioPage, SeriesScopingPage, SeriesBuilderPage, FocusPage, LoupePage) run impeccable audit + impeccable critique. First seed .impeccable/critique/ignore.md with SAXS domain vocabulary (Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square, Miller, q-value, hkl, reanalyze) so heuristic-2/jargon does NOT penalize correct domain language. Record /40, /20, and slug per surface into the tracker; merge audit P0/P1 findings into the backlog.
2. Otherwise pick the NEXT-HIGHEST-SEVERITY todo row (P0 before P1 before P2 before P3; within a tier, prefer items that fix multiple surfaces at once — e.g. the xl:→lg: rail fix spans Focus+Builder, the 31px display-xl step is app-wide).
3. If the item is a genuine PRODUCT/SCOPE decision, STOP and AskUserQuestion — do not guess. Mark the row blocked until answered.
4. Otherwise execute via the SUBAGENT-DRIVEN cadence, dispatched SEQUENTIALLY:
   a. Run verify-before-review against live source, then dispatch a fresh implementer subagent: give it the exact files, the mapped impeccable command + its reference discipline (from the synthesis doc), the relevant mockup at docs/redesign-mockups/*.html as the fidelity contract, and the constraint list above. The implementer applies the command, hands off to impeccable polish (design-system alignment first), and keeps appearance in ui primitives.
   b. Dispatch a spec-review subagent (does the change match the mockup/brief + resolve the backlog item + respect the laws?).
   c. Dispatch a code-quality-review subagent (closed-look/open-placement honored, import boundary, no class-string assertions, a11y).
   Address review findings (receiving-code-review: verify, don't perform agreement) before proceeding.
5. VERIFY BY RENDERING: serve a WRITABLE COPY of the dev DB via julia-from-source (never the prod DB; leave :8080 untouched), open the surface with Playwright MCP at the in-scope widths, screenshot, and confirm the fix renders and the mockup's major ingredients are present. A screenshot you didn't read back does not count.
6. KEEP THE GATE GREEN (from packages/HimalayaUI/frontend): npm run lint:design → npx tsc --noEmit -p tsconfig.build.json → targeted npm test (the affected spec files) → relevant npm run e2e spec → npm run build. The orchestrator (you) runs the slow suites; do not delegate slow gates to a subagent that will idle. If red, fix before committing.
7. COMMIT: stage ONLY the exact changed files by path (git add <path> <path> …). Write a focused message ending with the mandated co-author line. Do NOT push, do NOT open a PR, do NOT merge.
8. UPDATE tracker: mark the row done, record the command used and any follow-ups discovered.
9. SELF-PACE: ScheduleWakeup to continue the next iteration. End the turn.

RE-SCORE: after every ~4 completed P1/P2 items on a surface, re-run impeccable critique + audit for that surface and append the new score to the tracker trend.

FINISH LINE (stop the loop and report, do NOT merge): the agreed target is met (see the loop-config block at the top of loop-backlog.md). At that point post a summary of scores before/after per surface, the commits made (branch still unmerged), and the remaining P2/P3 backlog for Jonathan to triage. Surfaces still below target re-enter the P1/P2 waves instead of finishing.

GUARDED VERBS: do NOT run impeccable optimize unless a measured render bottleneck exists (measure first). Do NOT run impeccable overdrive unless the user explicitly asks for one extraordinary functional moment (then use its propose-2-3-then-ask gate). impeccable live is optional and only for Storybook leaf primitives — never bypasses the gate or commit rules.
```

The four scope inputs the loop hardcodes (responsive target, Duplicate-series feature, finish-line bar, scoping ordering variables) are settled in the §8 decision record below; the loop-backlog.md header records the resolved answers as a `loop-config` block.

---

## 8. Decision record (resolved with Jonathan, 2026-06-09)

These four were surfaced as genuine product/scope decisions before launch and answered. They are binding loop config.

1. **Responsive scope → desktop-first, honest down to ~1024px.** Fix the work·rail split breakpoint (`xl:`→`lg:`) on Focus + Builder and stack the Loupe grid below a documented min-width. Sub-tablet (<~768px) is explicitly out of scope — the `<600px` scoping foot-row wrap item is dropped, not fixed. Rationale: the verified defects bite in the 1024–1279px laptop-not-maximized range real scientists hit; true mobile is out of character for a dense instrument.

2. **Duplicate / Export → surface Export-figure now; defer Duplicate-series.** Export-figure is a pure presentation/placement P1 (`clarify` label + `layout` placement). Duplicate-series is potentially net-new product behavior (no duplicate path exists; may need backend/state) and leaves this loop as a separate `shape`-gated feature decision — not smuggled in through a fidelity fix.

3. **Finish bar → GOLD / flagship.** Every surface reaches **critique ≥ 36/40 AND audit ≥ 18/20** (the "excellent" bands) AND all P0/P1 backlog rows closed. Where a flagship item is genuinely infeasible without backend/feature work, mark the row `blocked` and surface it rather than forcing it — note blockers as discovered. Surfaces below target re-enter the P1/P2 waves.

4. **Scoping "Ordered by" variables → data-driven + custom.** The dropdown offers whatever ordering variables the backend/manifest actually exposes for the corpus, plus a "define your own" affordance. Verify the real variable set against the manifest/adapter before wiring — do not hardcode the mockup's time/dose/temperature example set (it would lie about capability).
