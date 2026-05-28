# R3-C — Contact-sheet keyboard a11y + chrome cleanup (incl. lock thumbs portrait)

**Issue:** #256 · **Branch:** `r3-c-samples-loupe-portrait` (off `main` @ `6d73d64`)
**Milestone:** HimalayaUI — The Print finish, round 3

## Goal

Close the round-3 residuals on the *entry surface* (contact sheet + loupe): one P1
accessibility blocker, four P2 finish items, and three P3 polish items. All findings
live in four files that ship on the first surface every workflow touches.

Source of truth: the #256 issue body (human-approved spec),
`docs/2026-05-28-the-print-round3-findings.md` (R3-S01/S02/S03/S04/S05/S07/S08 + U-3),
`DESIGN.md` §3 (Fixed-Scale Rule) + §5 (Inputs focus ring), and
`docs/redesign-mockups/sample-table.html`.

## styles.css reuse path (wave-2 disjointness)

**Confirmed: NO new `styles.css` role.** R3-S04 reuses the existing `text-headline-lg`
token (26px, already in the scale, added by R2-T9). The 1px delta from the inlined
`text-[25px]` is invisible at this size. `styles.css` is a shared global that wave-2's
sibling (#257) may also touch — I do not edit it. The Fixed-Scale Rule §3 names "the
progress numeral" explicitly as a serif/title use, and `text-headline-lg` is the serif
26px role, so this is a clean semantic match, not a workaround.

## Scope — files touched

| File | Findings |
|---|---|
| `components/ContactSheetRow.tsx` | R3-S01 (a11y `<button>` + ring + keydown), R3-S05 (accent demote), R3-S07 (tag-input ring), R3-S08 (`+ tag` opacity gate) |
| `components/LoupeSidebar.tsx` | R3-S02 (meta-list rename), R3-S03 (verdict `X` keycap), R3-S07 (tag-input ring) |
| `pages/SamplesPage.tsx` | R3-S04 (`text-[25px]` -> `text-headline-lg`) |
| `components/DetectorImage.tsx` | U-3 (lock thumbs portrait) |

Test files: `test/contact-sheet.test.tsx`, `test/LoupeSidebar.test.tsx`,
`test/DetectorImage.test.tsx` (all exist).

**R3-S06 scope (team-lead amendment 2):** R3-S06 (`DetectorImage.tsx:197`, the legacy
`text-fg-muted` "No image" survivor) is NOT in #256's issue body (verified: 0 matches).
The findings-doc recommended-actions table assigns it to Issue B (detector-window warmth).
I leave it untouched and note it DEFERRED to B in the PR -- even though I am in DetectorImage
for the portrait lock, touching it would cross the wave-disjointness boundary.

## Test-strategy constraints (frontend/test/AGENTS.md)

- **Never assert Tailwind class strings.** Visual findings -> stable `data-*` attributes
  or rendered-text / role / structural assertions.
- The a11y `<button>` promotion (R3-S01) is *behaviorally* testable -- assert role +
  keyboard activation, not classes.
- Regression floors/ceilings, not hard counts where a count is fragile.
- DetectorImage.tsx is inherited by wave-3 (#255) -- keep U-3 surgical: gate inside
  `evaluateOrient`, do NOT touch the shared `decideOrient` helper.

---

## Findings -> implementation

### R3-S01 (P1, a11y) -- `ExposureThumb` root -> `<button type="button">`

`ContactSheetRow.tsx:75-92`. Today the root is `<div onClick onDoubleClick>` with no
`role`/`tabIndex`/focus ring -- mouse-only.

Change:
- Root `<div>` -> `<button type="button">`. The onClick/onDoubleClick handlers move
  unchanged; a button's onClick already fires on Enter/Space (so plain-click *select* is
  keyboard-reachable for free).
- Add `focus-visible:ring-2 focus-visible:ring-print-accent focus-visible:ring-offset-1`
  to the className (the selection ring uses `ring-2 ring-print-accent` already; the
  focus-visible variant is additive and only paints on keyboard focus).
- Keyboard semantics (team-lead amendment 3 -- resolve the double-bind explicitly):
  a native `<button>`'s default activation fires `onClick` (= `onSelect`) on Enter/Space.
  The issue body prescribes Enter/Space -> `onOpenLoupe`. To honor that WITHOUT a
  double-action (select AND open on the same key), the `onKeyDown` handler for Enter/Space
  calls `e.preventDefault()` (suppresses the synthesized activation click, so `onSelect`
  does NOT also fire) then `onOpenLoupe()`. Net coherent model:
    * mouse single-click  -> select (toggle batch selection)
    * mouse double-click  -> open loupe
    * keyboard Enter/Space -> open loupe (the navigation affordance an AT user needs to
      screen a sample; matches the double-click = open analogue, since there is no clean
      single-key analogue of "double-click").
  This is asserted in the test: Enter opens the loupe AND the thumb is NOT left
  batch-selected (proving the auto-click select was suppressed).
- aria-label so AT announces the frame: `Frame ${frameNo}` (+ " (dropped)" when rejected).
- Add `appearance-none` defensively to neutralize UA button chrome (existing box classes
  already specify size/border/bg).

Tests (behavioral, not class):
- `exposure-thumb-N` is a `<button>` element (tagName === "BUTTON").
- keyDown Enter on the thumb navigates to the loupe (assert `loupe-stub`, reuse
  `renderRowRouted`) AND the thumb is NOT batch-selected afterward (preventDefault
  suppressed the auto-click select -> no double-action).
- keyDown Space also opens the loupe (same no-double-select assertion).
- Regression: existing click-select, shift-range, double-click-loupe tests stay green.

### R3-S02 (P2) -- Loupe meta-list: Filename/Kind -> Integration/Collected

`LoupeSidebar.tsx:233-247`. Mockup ships `frame / integration / collected / signal`
(lowercase, instrument-facing). "Kind" leaks SQLite nouns.

Change: the `This exposure` section's MetaRows become (lowercase labels, mockup order):
- `frame` -> frameLabel (existing `loupe-meta-frame`; relabel Frame -> frame)
- `integration` -> "-" stub (testid `loupe-meta-integration`)
- `collected` -> "-" stub (testid `loupe-meta-collected`)
- `signal` -> SignalMeter (existing; relabel Signal -> signal)

Remove the Filename (`loupe-meta-filename`) and Kind (`loupe-meta-kind`) rows entirely.
Stubbing to "-" matches the mockup's mock-data approach until backend integration/collected
fields are plumbed (out of scope here).

Tests (rendered text + testid):
- `loupe-meta-integration` present with text "-"; `loupe-meta-collected` present with "-".
- `loupe-meta-filename` and `loupe-meta-kind` gone (queryByTestId null).
- Labels render lowercase (getByText("integration"), etc.).
- Existing `loupe-meta-frame` / `loupe-meta-signal` tests stay green.

### R3-S03 (P2) -- Verdict toggle regains mono X keycap

`LoupeSidebar.tsx:272-279`. Mockup `v-toggle` reads `Drop X` / `Restore X` with X in mono.

Change: append `<span className="ml-1 font-mono text-[10px] opacity-60">X</span>` to the
Drop/Restore button (after the label).

Note: the `text-[10px]` here is the **issue-body-specified literal** (R3-S03's fix names
that exact span). The Fixed-Scale concern (R3-S04) targets the *progress numeral*, not this
keycap; the keycap idiom is reused verbatim across CullBar / footer-legend and the issue
prescribes it. I follow the issue body exactly. (Flagged so the verify-gate grep does not
mistake it for a survivor.)

Tests: `loupe-drop-toggle` text content includes X; the X is in a distinct mono span
(structural via getByText("X") within the button, not class-string). Kept + rejected both
show the keycap.

### R3-S04 (P2, Fixed-Scale) -- text-[25px] -> text-headline-lg

`SamplesPage.tsx:119` -- `<div className="text-display !text-[25px] text-ink">`.

Change: `text-display !text-[25px]` -> `text-headline-lg` (drops the inline override and
the now-redundant text-display; text-headline-lg is itself a Newsreader serif role at 26px,
so the serif treatment is preserved). Keep `text-ink`.

Tests: existing `screened-progress` count/bar tests assert content (numeral still shows
`screened / total`). The "zero inline text-[Npx]" Done-when is verified by grep at the
verify gate (per AGENTS.md no class-string DOM assertions), not as a runtime test.

### R3-S05 (P3) -- "N dropped" label: accent -> ink-soft + parenthesise

`ContactSheetRow.tsx:378-381`. Bare `text-print-accent` overspends accent on a passive fact.

Change: `text-print-accent` -> `text-ink-soft`; wrap label in parens: `({dropped} dropped)`.

Tests: existing kept-cell test asserts `1 dropped` text -> update to `(1 dropped)`. The
accent->ink-soft demote is visual-only, verified by build (no class assertion).

### R3-S07 (P3) -- Tag-form inputs gain a focus indicator

`LoupeSidebar.tsx:131-167` (loupe tag form) + `ContactSheetRow.tsx:413-457` (row tag form).
DESIGN.md S5 Inputs: 1px accent ring.

Change: on each tag `<input>` (key + value, both forms -- 4 inputs total) add
`focus:ring-1 focus:ring-print-accent/40 rounded-sm`. The inputs currently carry
`outline-none` with no replacement indicator.

Tests: focus-ring is a visual treatment -> no class-string assertion. Inputs stay functional
(existing tag-form open/submit tests cover behavior). Verified visually + by build.

### R3-S08 (P3) -- `+ tag` opacity: reach via row focus-within

`ContactSheetRow.tsx:472`. Today `opacity-0 group-hover:opacity-100 focus:opacity-100` --
the `focus:` cannot fire (button is opacity-0), and keyboard users cannot tab to a `+` they
cannot see arrive.

Change: `opacity-0 group-hover:opacity-100 focus:opacity-100` ->
`opacity-0 group-hover:opacity-100 group-focus-within:opacity-100`. The row wrapper already
carries `group` (`ContactSheetRow.tsx:295`), so any descendant receiving focus reveals `+`.

Tests: behavioral -- the non-empty `+` is a real reachable `<button>` in the DOM (so tab
order includes it). The opacity transition is CSS-only (JSDOM cannot compute it), so assert
presence + that it is a `<button>`. No class-string assertion.

### U-3 (P2, consistency) -- Lock contact-sheet thumbnails to portrait

`DetectorImage.tsx:184-210` rotation logic. The size-driven rotate makes thumbs inconsistent
row-to-row; the contact sheet should read as a uniform grid of portrait windows.

Change (surgical -- does NOT touch shared `decideOrient`): in `evaluateOrient`
(`DetectorImage.tsx:77-94`), short-circuit for `size === "thumb"`: set
`{ orient: "portrait", caps: null }` and return before calling decideOrient. Add `size` to
the useCallback deps. The initial useState default is already `{ orient: "portrait",
caps: null }`, so a thumb is portrait from first paint and never flips. `decideOrient`
(shared with DetectorRingOverlay, #180) is untouched -- the ring overlay only renders on the
loupe/focus full frame, so locking thumbs portrait has no cross-effect.

Tests (uses existing `data-orient` attribute on the wrapper):
- A `size="thumb"` DetectorImage renders `data-orient="portrait"` and never applies a
  rotate transform.
- Regression: `size="full"` orient logic unchanged (existing full tests stay green;
  decideOrient unit tests untouched).

---

## TDD task list (each step: failing test -> minimal impl -> verify -> commit)

1. R3-S01 a11y -- failing tests (thumb is `<button>`; Enter + Space open loupe); promote
   root to `<button>` + focus-visible ring + onKeyDown + aria-label; verify existing
   click/shift/double-click tests stay green; commit.
2. U-3 portrait lock -- failing test (thumb data-orient="portrait", no transform); gate
   `evaluateOrient` on `size === "thumb"`; verify full-frame regression green; commit.
3. R3-S02 meta-list -- failing tests (integration/collected stubs present, lowercase;
   filename/kind gone); rewrite the `This exposure` MetaRows; commit.
4. R3-S03 verdict keycap -- failing test (X mono keycap in Drop/Restore button); append the
   span; commit.
5. R3-S04 scale token -- change text-[25px] -> text-headline-lg; verify existing progress
   test green + grep shows zero text-[Npx] in SamplesPage; commit.
6. R3-S05/S07/S08 polish -- accent demote + parens (update kept-cell text test); tag-input
   focus rings (4 inputs); `+ tag` opacity -> group-focus-within; behavioral test that the
   non-empty `+` is a real reachable `<button>`; commit.

## Verify gate (MUST be green before PR)

```
cd packages/HimalayaUI/frontend
npm test           # Vitest -- all green (capture to file, grep FAIL/passed)
npm run build      # tsc --noEmit + vite build -- green
```

Plus visual Done-whens (below). Then request-pr-review.

## How VISUAL Done-whens are verified

- Keyboardable thumbnails (R3-S01): behavioral test -- element-type + Enter/Space fires
  onOpenLoupe (router nav to loupe-stub). Strongest Done-when, fully automated.
- Meta-list labels (R3-S02): rendered-text assertions (integration/collected lowercase
  present, filename/kind absent) -- text, not class.
- Verdict Drop X keycap (R3-S03): rendered-text + structural (X in a distinct mono span).
- text-[Npx] gate (R3-S04 Done-when, NARROWED per team-lead amendment 1): the three target
  files already carry ~25 PRE-EXISTING pixel literals (text-[10.5px], text-[11.5px], ...)
  that belong to F/#259's project-wide Fixed-Scale sweep and are OUT of #256 scope -- do NOT
  remove them. The gate is: (i) the R3-S04 `text-[25px]` survivor is GONE, and (ii) no NEW
  pixel literal is introduced beyond the issue-prescribed R3-S03 `text-[10px]` keycap. Verified
  by diffing the literal set before/after; I report the diff in the PR.
- Thumbs uniformly portrait (U-3): data-orient="portrait" + no rotate transform -- automated.
- Accent demote / focus rings / opacity gate (R3-S05/S07/S08): CSS-only visual treatments
  JSDOM cannot compute; verified by `npm run build` green + visual check; the behavioral
  parts (parens text, `+` reachable) are unit-tested. Reported as "build-verified + visually
  confirmed", not claimed as runtime-tested for the pure-CSS deltas.

## Risks / notes

- The `<button>` swap (R3-S01) is the only change with regression surface -- a `<button>`
  inside the overflow-x-auto flex strip must keep the same box model. Existing classes
  (relative h-[62px] w-[62px] shrink-0 ... rounded-[3px] border bg-frame-edge) fully specify
  the box; button UA defaults are overridden; `appearance-none` added defensively. Verify via
  build + visual + the surviving click/shift/double-click tests.
- DetectorImage edit is contained to `evaluateOrient` + its dep array -- #255 (wave-3)
  inherits a minimal, well-commented diff.
- No `styles.css` edit. No new scale role. Reuse path confirmed.
