# R5 — Focus affordances (issue #228) implementation plan

Branch: `r5-focus-affordances` off `origin/main`. Scope = F-11/F-12/F-13/F-14 only.
Mockup: `docs/redesign-mockups/focus-workspace.html`. Findings: `docs/2026-05-27-the-print-fidelity-findings.md` §3.

## Design decisions (ambiguity resolved here)

- **Source of truth for the focus surface is the URL** (`/sample/:sampleId`), synced one-way
  into Zustand by `useSyncActiveSampleFromRoute`. So the **stepper navigates the URL**
  (`useNavigate('/sample/:nextId')`), not the store directly — consistent with the existing
  one-way shim and avoids re-introducing a store→URL handshake.
- **Sample ordering for the stepper** = corpus samples filtered to the active sample's
  `experiment_id` (matches the `,`/`.` shortcut's experiment-scoped semantics). No wrap.
- **Representative = the `selected` exposure.** A real persistence mechanism already exists:
  `useSelectExposure(sampleId)` → `PATCH /api/exposures/:id/select`. So F-11 set-rep persists;
  click-to-switch is local view state via `setActiveExposure`. No persistence gap.
- **F-12 Notes drawer**: the topbar hosts a `Notes` toggle button (visible only on `/sample/:id`).
  Toggle state is an ephemeral Zustand flag `notesDrawerOpen`. The drawer body is rendered inside
  `FocusWorkspaceLayout` (which already wires `notesSample` + `updateSample` from the correct
  experiment-scoped cache) so the drawer and the xl margin share one save path. Mirrors the
  mockup: column present `>= xl`, button+drawer `< xl`. The button shows a `•` badge when notes exist.
- **F-14 stage tab**: #181 is CLOSED — the focus workspace cutover shipped; `/index*` redirects to
  `/sample/...`. The Index stage tab in `CorpusTopbar` is still a permanently-disabled button.
  Make it ACTIVE (not a link — there is no canonical `/index` route anymore; it redirects) when the
  route is `/sample/:id`. Minimal, no conflict with #181's retired model.

## Steps (each its own commit, TDD)

1. **state**: add `notesDrawerOpen` + `openNotesDrawer`/`closeNotesDrawer`/`toggleNotesDrawer`
   (ephemeral, not persisted). Closed on sample switch (`setActiveSample`). Test: state.test.
2. **F-13 + F-14 (CorpusTopbar)**: per-sample stepper (prev/next, "sample N of M", display name)
   shown only on `/sample/:id`; Index tab active on that route. Test: CorpusTopbar.test +
   new stepper test. Stepper reads `useCorpusSamples`, navigates URL.
3. **F-12 (Notes toggle + drawer)**: topbar Notes button (sample route only) toggles
   `notesDrawerOpen`; `FocusWorkspaceLayout` renders an `xl:hidden` drawer when open. Tests:
   CorpusTopbar (button visibility/badge), FocusWorkspaceLayout (drawer reachable below xl).
4. **F-11 (rep-exposure switcher)**: thumbnail strip in `FocusDetectorPanel` header — mini
   detector thumbnails per exposure, click switches (`setActiveExposure`), set-rep persists
   (`useSelectExposure`), rep dot on `selected`. Test: FocusDetectorPanel.

## Verify
`npm run build` + `npm test` green. Live screenshots at 127.0.0.1:5206 /sample/2 incl.
narrow-viewport (<1320px) Notes reachability.
