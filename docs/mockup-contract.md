# Mock-up behaviour contract

The workspace is built to match the approved mock-up 1:1, with app-only features replayed on top. This file lists every behaviour the mock-up's code defines: its key handlers, click handlers and render branches, not just its look. Each row says where that behaviour stands in the app. `npm run lint:contracts` (part of `npm run build`) checks the rows, so a behaviour can't drop out of the app without a visible change to this file.

- Mock-up: https://claude.ai/artifact/UoZAK7fkkvCBQ13TQYUXVH
- Mock-up version: `1795147819-6f12`
- Mock-up sha256: `991c30e426078ab8a1bb57cf09b9a259002b201b16ff07af3c58809f98da57b4`

When the mock-up changes, extract its behaviours again from the code and update the version, the hash and the rows in the same PR. The consistency audit compares the pinned version with the live artifact each wave.

## Statuses

| Status | Meaning | Evidence the lint requires |
|---|---|---|
| `built` | The app does it, and a test proves it. | `path › test title`, a path under `packages/HimalayaUI/frontend/` whose file contains that title. |
| `untested` | The app does it, but no test pins it yet. | What shows it works. Turning a row `built` means writing its test. |
| `superseded` | A later product decision replaced it. | The decision: a date, a PR number or a doc path. |
| `planned` | Not built yet, and someone owns it. | The owner: a PR number, an issue (`#n`) or a thread name. |

## Sample list and header

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| list-arrows | ↑ ↓ move to the previous or next sample in the filtered list | built | test/workspace/WorkspacePage.test.tsx › arrow down steps to the next sample in the list |
| list-next-unreviewed | A key jumps to the next unreviewed sample, wrapping around (mock-up J) | superseded | Z, not J (Jonathan 2026-11-21, Gradescope convention); test/workspace/WorkspacePage.test.tsx covers Z |
| list-prev-unreviewed | Shift with that key jumps to the previous unreviewed sample | planned | held small item from the 2026-11-23 mock-up audit (thread "Audit app against the mock-up") |
| list-filter | / focuses the filter; it matches names; Enter opens the first match; Esc leaves it | built | test/workspace/WorkspacePage.test.tsx › / puts the cursor in the filter, where letters type, and ? opens the shortcuts only while the filter is empty |
| list-no-match | An empty filter result says "No samples match." | built | test/workspace/WorkspacePage.test.tsx › says no samples match only when a filter is typed |
| list-status | Each row shows its call, or an unreviewed mark | built | test/workspace/WorkspacePage.test.tsx › updates the list and the count when the corpus status changes |
| list-every-phase | A coexistence call lists every phase in the row | planned | PR 502 |
| list-suspect-mark | A row with an undropped suspect image carries a warning mark | planned | thread "Show combined curve by default" |
| list-follows-call | The row and the progress follow the call as it changes, without a reload | planned | PR 507 (a peak edit that drops the call left them stale) |
| list-cursor | The current row is highlighted and scrolled into view | untested | SampleList scrolls the current row into view on every sample change |
| header-progress | The header shows a progress bar and "k of N reviewed" | built | test/workspace/WorkspacePage.test.tsx › updates the list and the count when the corpus status changes |
| header-help | The header's ? button opens the shortcut list | built | test/workspace/WorkspacePage.test.tsx › ? opens the shortcut list, the header button opens it too, and Escape closes it |
| layout-narrow | Below 980 px the tools pane moves under the plot; below 640 px one column | planned | held small item from the 2026-11-23 mock-up audit (thread "Audit app against the mock-up") |

## Images

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| img-arrows | ← → move between images | built | test/workspace/WorkspacePage.test.tsx › arrow keys step from this sample's last frame into the same sample's frames |
| img-drop | X or the card's Drop button drops or restores the image | built | test/workspace/WorkspacePage.test.tsx › X drops the shown frame and shows the nearest kept frame, so a phase key calls there |
| img-keep-one | Dropping the last kept image is refused ("Keep at least one exposure") | planned | thread "Show combined curve by default" |
| img-kept-hint | The section header says "k of N kept" once one is dropped | built | test/workspace/WorkspacePage.test.tsx › a dropped card offers Restore and the header counts the kept frames |
| img-dropped-card | A dropped card fades and says "dropped" | built | test/workspace/WorkspacePage.test.tsx › a dropped card offers Restore and the header counts the kept frames |
| img-enlarge | Space enlarges the image; Space or Esc closes it | built | test/workspace/WorkspacePage.test.tsx › Space enlarges the shown frame's detector image, and Space or Esc closes it |
| img-enlarge-keys | Inside the enlarged view ← → step and X drops | built | test/workspace/WorkspacePage.test.tsx › flips through frames inside the enlarged view, and dropping one moves to the nearest kept frame without closing it |
| img-dblclick | Double-clicking a card enlarges it | planned | held small item from the 2026-11-23 mock-up audit (thread "Audit app against the mock-up") |
| img-overlay-mean | The plot draws every kept image faintly with their mean in bold | planned | thread "Show combined curve by default" |
| img-show-only | O or "Show only" plots one image, with a banner; O again shows all | planned | thread "Show combined curve by default" |
| img-suspect | A suspect image (low-q excess) is outlined, says why, and is drawn in the warning colour | planned | thread "Show combined curve by default" |
| img-open-suspect | Opening a sample selects its first undropped suspect image | planned | thread "Show combined curve by default" |
| img-rings | Hovering a candidate draws its predicted rings on every detector image | planned | thread "Audit app against the mock-up" |
| img-meta | The Image panel lists the shown image's file, time, exposure, I₀, I₂ and flight tube | built | test/workspace/WorkspacePage.test.tsx › lists the shown exposure's time, I₀ and I signal under the .prp names |
| img-transmission | The Image panel shows I₂ / I₀ | planned | held small item from the 2026-11-23 mock-up audit, pending Jonathan's call (thread "Audit app against the mock-up") |
| img-lowq | The Image panel shows "Low-q vs others" | planned | thread "Show combined curve by default" |

## Plot and peaks

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| plot-zoom | Scrolling zooms q | built | test/workspace/WorkspacePage.test.tsx › drops a kept zoom that falls outside the next sample's q range, and ⌘/Ctrl+0 resets the zoom |
| plot-zoom-anchor | The zoom keeps the q under the cursor in place | built | test/plot/interaction.test.ts › keeps the q under the cursor at the same place on the axis |
| plot-reset-zoom | A key resets the zoom (mock-up Z) | superseded | ⌘/Ctrl+0 (Jonathan 2026-11-21 keymap) |
| plot-readout | The cursor reads out q and d | built | test/workspace/WorkspacePlot.test.tsx › shows q and d under the cursor and nothing once it leaves |
| plot-peak-click | A click near an auto peak turns it off, near a manual one removes it, elsewhere adds one | built | test/workspace/WorkspacePage.test.tsx › a click near an auto peak disables it, near a manual peak removes it, and elsewhere adds one |
| plot-peak-markers | Peak markers sit in a row above the trace; off peaks are outlined | built | test/workspace/WorkspacePlot.test.tsx › draws the trace multiplied by q² in the q²I view and keeps the peak markers in a row above the plot |
| plot-comb-hover | A candidate's reflections are drawn only while it is hovered | superseded | Jonathan 2026-11-23, PR 500: the called phase's reflections stay drawn, faded; see plot-comb-faded-call |
| plot-comb-faded-call | The called phase's reflections stay drawn faded, and a hovered candidate's are drawn at full strength | built | test/workspace/WorkspacePage.test.tsx › the plot keeps the called phase's reflections faded, and draws a hovered candidate's at full strength instead |
| plot-log-linear | L toggles log and linear intensity | superseded | views live in the palette and L opens labels (Jonathan 2026-11-21, PR 390) |
| plot-peak-edit-call | A peak edit keeps each ticked phase on the new candidate of that phase | superseded | the call is re-attached by phase and lattice within 5%, and a drop is announced (docs/event-log.md, MEMBER_REATTACH_RELTOL) |
| plot-peak-edit-undo | Undoing a peak edit also restores the call it dropped | planned | thread "Audit app against the mock-up" (follow-up to PR 507) |

## Comparisons

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| cmp-pick | C opens a picker; picks are drawn dashed, one colour each | built | test/workspace/WorkspacePage.test.tsx › C opens a picker, and Enter adds the sample scaled to the current trace without subtracting it |
| cmp-remove | A chip's × removes the comparison | built | test/workspace/WorkspacePage.test.tsx › Ctrl-clicking the subtracted chip stops subtracting, the palette hides comparisons, and × removes one |
| cmp-chip-click | Clicking a chip makes it the subtraction background | superseded | a click shows or hides the chip, ⌘-click subtracts (compare-sidebar decision, 2026-11-21) |
| cmp-strip | A difference strip sits under the plot while a background is set | superseded | Raw/Minus toggle replaces the strip (compare-sidebar decision, 2026-11-21) |
| cmp-scale | Comparisons are scaled to the current trace at high q | superseded | scaled by .prp beam monitors, high q as the fallback (buffer-choice decision, PR 325) |
| cmp-hide-all | R hides or shows every comparison | superseded | palette only (Jonathan 2026-11-21 keymap) |
| cmp-carry | Comparisons are kept as you move between samples | superseded | comparisons belong to the sample (per-sample-state decision, PR 422) |
| cmp-mean | A comparison draws the other sample's kept mean | planned | thread "Show combined curve by default" |

## Call

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| call-keys | 1–9 tick or untick a candidate, saved at once | built | test/workspace/WorkspacePage.test.tsx › ranks candidates by score and number keys toggle them on the shown exposure |
| call-row | A candidate row shows its key, phase chip, lattice and score, or "Custom" | built | test/workspace/WorkspacePage.test.tsx › an added custom index is listed as Custom and ticked, though it has no score |
| call-match-count | A candidate row says how many reflections matched | planned | held small item from the 2026-11-23 mock-up audit (thread "Audit app against the mock-up") |
| call-empty | No candidates says "Broad bumps only, no Bragg peaks." or "No peaks." | planned | held small item from the 2026-11-23 mock-up audit (thread "Audit app against the mock-up") |
| call-no-phase | Unticking every phase leaves the sample reviewed with no phase | superseded | 0 marks no phase (form factor) (Jonathan 2026-11-20) |
| call-index-chart | The indexing chart shows the hovered candidate's fit, else the ticked ones | built | test/workspace/WorkspacePage.test.tsx › hovering a candidate shows its fit in place of the ticked ones |
| call-index-legend | The chart legend gives each fit's lattice, R² and reflections on peaks | built | test/workspace/WorkspacePage.test.tsx › the indexing plot shows the assigned candidate's fit against the peaks |

## Custom index

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| cix-open | I opens a custom index on the top candidate | built | test/workspace/WorkspacePage.test.tsx › I opens a custom index on the top candidate; arrows nudge the lattice and Enter adds it |
| cix-controls | Symmetry buttons, a lattice slider and a number field | built | test/workspace/WorkspacePage.test.tsx › a lattice typed outside the range is flagged and blocks Add until it is clamped on blur |
| cix-landed | It says how many reflections land on a peak | built | test/workspace/WorkspacePage.test.tsx › the indexing plot follows the custom index while it is open |
| cix-nudge | ← → nudge the lattice, Shift for a larger step | built | test/workspace/WorkspacePage.test.tsx › arrows nudge the lattice after a symmetry click, while focus is on that button |
| cix-snap | Clicking near a peak puts the first order on it | built | test/workspace/WorkspacePage.test.tsx › while a custom index is open, a plot click puts the first order on the nearest peak and edits no peaks |
| cix-enter-after-snap | Enter adds right after a click snapped the lattice, whatever had focus | planned | PR 505 |
| cix-esc | Esc cancels without writing | built | test/workspace/WorkspacePage.test.tsx › Escape cancels a custom index without writing |
| cix-ticks | Adding it ticks it into the call | built | test/workspace/WorkspacePage.test.tsx › an added custom index is listed as Custom and ticked, though it has no score |

## Palette, undo and help

| id | Behaviour | Status | Evidence |
|---|---|---|---|
| pal-open | ⌘K opens the palette; typing a sample's name and Enter goes there | built | test/workspace/WorkspacePage.test.tsx › ⌘K opens the command palette and Enter goes to the matching sample |
| pal-same-action | A palette command runs the same action as its key | built | test/workspace/WorkspacePage.test.tsx › Ctrl+K toggles the palette, and a command runs the same action as its key |
| pal-theme | The palette switches the theme | built | test/workspace/WorkspacePage.test.tsx › the palette's theme page switches the workspace to dark and back, remembers it, and can follow the system again |
| undo-key | A key undoes the last drop, call or peak edit, and a toast names it (mock-up U) | superseded | ⌘/Ctrl+Z, reading the server first (Jonathan 2026-11-21 keymap); test/workspace/WorkspacePage.test.tsx covers it |
| undo-nothing | Undo with nothing to undo says so | built | test/workspace/WorkspacePage.test.tsx › says so when there is nothing to undo |
| help-list | ? opens the shortcut list | built | test/workspace/WorkspacePage.test.tsx › ? opens the shortcut list, the header button opens it too, and Escape closes it |
| help-keys | Shortcuts do nothing while the list or a dialog is open | built | test/workspace/WorkspacePage.test.tsx › shortcuts do nothing while the shortcut list is open |
