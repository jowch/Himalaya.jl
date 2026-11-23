# Product decisions

This file lists the product decisions Jonathan has made, one row each, and says where each one stands in the code. Read it before you build or review a feature: if a change contradicts a row, the row changes in the same PR, with the new decision's date. `npm run lint:contracts` (part of `npm run build`) checks every row, so a decision can't be recorded and then quietly go unbuilt, and a built one can't lose its test without a visible change here.

The mock-up's own behaviours are listed in [mockup-contract.md](mockup-contract.md). This file holds the decisions that sit above or beside the mock-up. Where a decision replaced a mock-up behaviour, that row in the contract is `superseded` and cites the decision.

## Statuses

| Status | Meaning | Evidence the lint requires |
|---|---|---|
| `built` | The app does it, and a test proves it. | `path › test title`, a path from the repository root to a Vitest or Playwright file holding `it`/`test` with that title, or a Julia test file holding `@testset` with that title. |
| `untested` | The app does it, but no test pins it yet. | What shows it works. Turning a row `built` means writing its test. |
| `superseded` | A later decision replaced it. Keep the row so the history stays readable. | The replacing decision: a date, a PR number or a doc path. |
| `planned` | Decided but not built, and someone owns it. | The owner: a PR number, an issue (`#n`) or a thread name. |
| `contradicted` | The app does the opposite today. A known bug, not a gap. | The owner: a PR number, an issue (`#n`) or a thread name, and where the app contradicts it. |

Every row also needs the date it was decided, as `YYYY-MM-DD`.

## Who the app is for

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| user-one-sample | The app serves an experienced user triaging a fresh beamtime fast, one sample at a time. | 2026-11-19 | untested | The workspace opens on the first unreviewed sample (docs/mockup-contract.md list rows) |
| mockup-one-to-one | The workspace matches the approved mock-up 1:1, and app-only features are replayed on top. | 2026-11-21 | untested | docs/mockup-contract.md tracks each mock-up behaviour row by row |
| compare-zero-setup | Comparing with another sample needs no setup: pick it and it is drawn. | 2026-11-19 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › C opens a picker, and Enter adds the sample scaled to the current trace without subtracting it |

## Keys

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| keys-standard | Standard keys only, and nothing advances to the next sample by itself. | 2026-11-21 | untested | No handler moves to another sample after a call; the mock-up contract's list rows cover the keys that do |
| keys-next-unreviewed | Z goes to the next unreviewed sample. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › Z jumps to the next unreviewed sample, like Gradescope's next ungraded |
| keys-reset-zoom | ⌘/Ctrl+0 resets the zoom. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › drops a kept zoom that falls outside the next sample's q range, and ⌘/Ctrl+0 resets the zoom |
| keys-undo | ⌘/Ctrl+Z undoes and ⌘/Ctrl+Shift+Z redoes. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › Ctrl+Shift+Z redoes what was undone |
| keys-label | L opens the label picker. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › L opens the label picker without typing the L, and Enter confirms a suggestion |
| keys-space | Space opens the enlarged detector view. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › Space enlarges the shown frame's detector image, and Space or Esc closes it |
| keys-plot-view | V and Shift+V step through the plot views, from a switcher in the plot's top-right corner. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › V steps forward through the plot views and Shift+V steps back, wrapping |
| keys-lists | In lists Tab and Shift+Tab move, Enter commits and Esc backs out. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/Palette.test.tsx › Tab, Shift+Tab and the arrows move the highlight and wrap, with focus kept in the field |
| keys-lists-filter | The sample filter follows the list rule too: Tab and Shift+Tab move the highlight. | 2026-11-21 | planned | thread "Find decisions that were never built": the filter has Enter for the first match and Esc, but no Tab highlight |
| keys-lists-symmetry | The custom index's symmetry choice follows the list rule. | 2026-11-21 | planned | thread "Find decisions that were never built": the symmetries are a plain button group |
| keys-keycap-every-control | A control with a shortcut shows its key as a keycap, except Drop (keys-drop-no-keycap). | 2026-11-20 | untested | Buttons with shortcuts render the Kbd primitive; no test checks every one |
| keys-drop-no-keycap | The Drop button carries no keycap, which would be noise; its X shortcut is listed in the ? shortcut list. | 2026-11-23 | untested | The shortcut list in packages/HimalayaUI/frontend/src/workspace/WorkspacePage.tsx has the X row; the cited ? test checks only that the list opens |
| keys-tab-assignment | In the workspace Tab and Shift+Tab cycle only the assignment choices, and Esc frees Tab until the next sample. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › Tab and Shift+Tab cycle only the phase choices, wrapping |
| keys-drop-moves | Dropping an image moves to the nearest kept one. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › X drops the shown frame and shows the nearest kept frame, so a phase key calls there |

## Undo

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| undo-reads-server | Undo reads the server first and never blindly replays an edit someone else has since changed. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/editOps.test.tsx › reports nothing sent when the edit is already undone on the server |
| undo-experiment-meta | Renaming the experiment and editing its description are undoable. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › undoes and redoes renaming the experiment and changing its description from the palette |
| undo-setup-edits | Hiding and renaming a setup are undoable. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/setupUndo.test.ts › undoes a hide and a rename when the server still shows the edit |
| undo-peak-restores-call | Undoing a peak edit also restores a call the edit dropped. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › undoing a peak edit that dropped the call ticks the phase again on the re-fitted candidate |

## Wording

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| word-image | The UI calls one detector image of a sample an "image". "Exposure" stays only for exposure time and in code. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/eslint-ui-copy.test.ts › enforces the house vocabulary in user-facing copy only |
| word-frame | The UI calls one detector image a "frame". | 2026-11-21 | superseded | word-image, 2026-11-22 (PR 480) |
| word-label | "Tag" is renamed "Label". | 2026-11-21 | built | packages/HimalayaUI/frontend/test/eslint-ui-copy.test.ts › enforces the house vocabulary in user-facing copy only |
| word-assignment | The tools column's Call section is titled "Assignment", with Custom index in its header. Other uses of "call" stay. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › the Custom index button sits in the Assignment header and stays there while the custom index it opens replaces the list |

## Look and layout

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| look-closed-primitives | Appearance lives in the src/ui primitives; other code passes placement-only classes. | 2026-05-29 | built | packages/HimalayaUI/frontend/test/check-design.test.ts › flags a pixel font size with no line height in the same class string, even in src/ui |
| look-type | Source Sans 3, no monospace, four sizes, weights 400 and 600, whole-pixel line heights. | 2026-11-21 | untested | packages/HimalayaUI/frontend/src/styles.css font tokens |
| look-themes | Light mode is the target; dark mode is a charcoal blue-violet with the same tokens for the system and explicit choice. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/theme.test.ts › states the same dark tokens for the system setting and the explicit choice |
| look-phase-status-hues | Phase hues never reuse a status hue; error is crimson so it stays apart from the terracotta accent. | 2026-11-21 | untested | packages/HimalayaUI/frontend/src/styles.css colour tokens (PR 405) |
| look-beam-off-warning | A beam-off image is flagged like a flare: an outlined tile with the reason, and a dot on the sample row. | 2026-11-21 | planned | thread "Find decisions that were never built": only the plot caption shows it today; it belongs with the suspect-image warning of thread "Show combined curve by default" |
| look-selection | A ticked phase choice is one flat tint, with no outline. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/CandidateList.test.tsx › fills a ticked row instead of outlining it, so selection never reads as keyboard focus |
| look-selection-list | The current sample row is one flat tint, with no accent bar. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/SampleList.test.tsx › tints the current row flat, like an Assignment choice, with no accent bar |
| look-no-mouse-ring | A ring shows only for keyboard focus: a row or image picked with the mouse shows none, even after the arrow keys move the selection. | 2026-11-23 | built | packages/HimalayaUI/frontend/e2e/workspace-focus.spec.ts › a clicked sample row shows no focus ring after the arrow keys move the selection |
| layout-no-shift | Nothing in the workspace may shift between samples or as data loads. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › lists every exposure field, with a dash for the ones the frame lacks, so the section keeps its height |
| layout-tools-tall | On a tall window the tools column runs Call, then Indexing, then the image. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › puts the Assignment list first in the tools column, then the indexing plot, then the exposure |
| layout-tools-short | On a short window (under 860 px) it runs Call, then the image, then Indexing, with the chart pinned; Call and the image fold as an accordion. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › folds Assignment and Image as an accordion on a short window, and opening a custom index from the folded Assignment header shows it in place of the list |
| plot-fixed | The plot has a fixed size; scrolling zooms q and a click adds or disables a peak. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › a click near an auto peak disables it, near a manual peak removes it, and elsewhere adds one |
| plot-hover | The plot's cursor is a crosshair; the readout snaps to a hovered peak; the called phase's reflections stay faded when nothing is hovered. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePlot.test.tsx › snaps the crosshair and readout to a peak's q while the pointer is on its marker, and not in the plot below it |
| tooltip-plate | Hover hints use the Tooltip primitive in the plate look; a DOM title is banned. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/ui/Tooltip.test.tsx › does not open when a mouse click focuses the trigger |
| tooltip-timing | A tooltip opens after 400 ms, at once while another was just open, and at once on keyboard focus. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/ui/Tooltip.test.tsx › opens a neighbour instantly while warm, and waits again once cold |

## Per-sample state

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| state-per-sample | Phase calls, comparisons, the subtraction background and labels belong to the sample; switching samples carries nothing over. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › comparisons and the subtraction belong to the sample: moving on carries nothing, and coming back restores them |
| state-call-survives-drop | A phase call belongs to the sample and survives dropping the image it was made on. | 2026-11-21 | built | packages/HimalayaUI/test/test_sample_calls.jl › a phase call survives dropping the frame it was made on (#350) |
| state-recents-local | A local recents cache only ranks picker suggestions and stays in this browser. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › the label picker lists labels added recently first, and the choice stays on this browser only |
| state-compare-cap | A sample holds at most five comparisons, one per comparison colour. | 2026-11-22 | built | packages/HimalayaUI/test/test_routes_compares.jl › a sample holds at most five comparisons, one per comparison colour |
| state-compare-live | Comparisons are shared live; chips other people add are hidden by default and fetch nothing until shown. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › a comparison someone else added arrives hidden, and yours start shown |

## Comparisons and subtraction

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| cmp-sidebar-setups | The sample list is grouped by setup. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/SampleList.test.tsx › heads each run of a setup's samples with its name |
| cmp-sidebar-click | ⌘/Ctrl-clicking a sidebar row toggles it as a comparison without leaving the sample. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › Ctrl/Cmd-clicking a sidebar row toggles it as a comparison without leaving the sample |
| cmp-chip-gestures | A click on a chip shows or hides it; ⌘/Ctrl-click makes it the subtraction background. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › a plain click on a chip hides or shows just that trace and keeps the subtraction |
| cmp-raw-minus | The Raw/Minus toggle shows only while a background is set. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › Ctrl-click subtracts in the main plot, Raw brings the trace and comparisons back, and a dimmed chip goes back to Raw |
| cmp-pseudolog | The difference is drawn on an asinh pseudolog axis in the log views. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/difference.test.ts › is linear in Linear and q²I, and pseudolog in the log views |
| cmp-chip-setup | A chip from another setup carries a small setup icon, with the full name in its tooltip; numbered badges are rejected. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/ui/CompareBar.test.tsx › names the chip and its setup and says how to hide and subtract |
| cmp-chips-visible | Chips are the plot's legend, so none may hide or scroll away: the bar reserves two rows and wraps. | 2026-11-22 | untested | Compare bar layout (PR 477) |
| buffer-no-default | No buffer is picked automatically; labelled buffers rank first and the user picks. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › the compare picker lists buffers first, marks which are only suggested, and picks nothing by itself |
| buffer-monitor-scale | Subtraction scales by the .prp beam monitors, with a high-q match as the fallback. | 2026-11-20 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › scales by the transmitted beam when both traces carry monitors, even where high-q matching would disagree |
| flight-path-label | The flight path is an auto-suggested label. | 2026-11-21 | built | packages/HimalayaUI/test/test_label_suggestions.jl › suggests the flight path when every frame of a sample agrees |

## Labels

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| label-colours | Labels use eight colours, kept distinct within labels; they may overlap phase or status hues. | 2026-11-21 | built | packages/HimalayaUI/test/test_routes_label_colors.jl › the eight colours match the frontend's list, in order |
| label-system | Buffer (blue), Empty capillary (stone) and Calibrant (ochre) are fixed system labels, added per sample. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/labels.test.ts › a stored colour never overrides a system label's colour |
| label-suggestions-stay-removed | A label suggestion the user removed is never suggested again. | 2026-11-21 | built | packages/HimalayaUI/test/test_routes_labels.jl › a removed suggestion is rebuilt from the log as removed |

## Palette

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| palette-widgets | Every ⌘K command is a list, a single field or an inline form; Esc goes back one level. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/Palette.test.tsx › turns a Create row into an inline form: Tab moves between fields, arrows change a choice, Enter submits |
| palette-settings | Settings sit behind one palette entry, each edited in its own form. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/workspacePalette.test.ts › lists each setup on the Settings page and saves a rename or a hide from its form |
| palette-custom-index | The custom index modal moves into the palette as a widget. | 2026-11-21 | planned | thread "Design palette with inline widgets" |

## Images and combining them

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| combine-include-toggle | Each image has an include toggle; the combined trace is a time-weighted sum, or the median. | 2026-11-20 | superseded | combine-select, 2026-11-23 |
| combine-select | What you select is what is analysed. Images are selected by a checkbox, Shift-click or ⌘-click; one image is that image, several are their monitor-normalised average with its own peaks, keyed by the set. Every kept image starts selected, with Select all and Select none. A plain click only views. The call stays with the sample. | 2026-11-23 | built | packages/HimalayaUI/test/test_combine.jl › refresh_selection! averages every selected image, one curve per set |
| combine-keep-one | Dropping the last kept image is refused. | 2026-11-23 | superseded | combine-drop-last, 2026-11-23 |
| combine-drop-last | The last kept image can be dropped. Nothing is analysed, and the sample row says why ("All dropped" or "None selected") instead of blocking the drop. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/SampleList.test.tsx › says on the row why nothing is analyzed |
| combine-restore-unselected | A dropped image that is restored comes back unselected, unless nothing else is selected. Undo of a drop puts back the selection it had. | 2026-11-23 | built | packages/HimalayaUI/test/test_combine.jl › selection, drop and restore routes |
| combine-beam-off-deselected | Beam-off images start deselected, so they are averaged in only when ticked. | 2026-11-23 | built | packages/HimalayaUI/test/test_ingest_reconcile.jl › a beam-off image starts left out of the selection, at ingest and when its trace lands late |
| combine-no-cross-setup | Images never combine across setups. | 2026-11-22 | built | packages/HimalayaUI/test/test_combine.jl › a bridged sample averages only its own images |

## Bridges and beamtimes

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| bridge-is-merge | A bridge says two samples are the same sample measured in two setups. It is not a comparison and takes no comparison slot. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › the same sample takes no comparison slot, so a fifth comparison goes in past a bridge |
| bridge-note | A bridged sample shows "Same as X in SETUP" and its mate's images under their own setup heading. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › shows a same sample's frames under their own heading after this sample's, and picking one shows its image |
| bridge-hidden-setup | A bridge to a sample in a hidden setup stays stored but is not drawn. | 2026-11-21 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › draws no bridge to a sample whose setup is hidden, and draws it again once the setup is shown |
| bridge-shared-call | Bridged samples share one call and their labels, except flight paths. | 2026-11-22 | built | packages/HimalayaUI/test/test_bridge_calls.jl › bridged samples share one call |
| bridge-keep-call | Opening a bridged pair whose calls differ asks once which call to keep. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › opening a bridged sample whose call state differs from its mate's asks once which call to keep |
| bridge-export | The export keeps one row per sample and names its bridge mate. | 2026-11-22 | built | packages/HimalayaUI/test/test_routes_bridges.jl › the export names the other bridged sample and its setup, never the row itself |
| bridge-merged-into | The old merged_into_id column is dropped; samples it emptied stay as rows. | 2026-11-22 | built | packages/HimalayaUI/test/test_ingestion_structural_events.jl › the retired samples.merged_into_id column is dropped and its rows kept |
| beamtime-one-workspace | One workspace per beamtime folder, with one setup per subfolder. | 2026-11-21 | built | packages/HimalayaUI/test/test_beamtime.jl › POST /api/experiments with setups makes one beamtime |
| beamtime-rescan-adds | A rescan adds a new setup folder to the beamtime. | 2026-11-21 | built | packages/HimalayaUI/test/test_beamtime.jl › a rescan adds a new setup folder to the beamtime |
| beamtime-merge-cli | `himalaya merge` folds per-setup experiments into one beamtime, as a dry run unless --apply. | 2026-11-21 | built | packages/HimalayaUI/test/test_beamtime.jl › merge folds per-setup experiments into one beamtime |

## Phases

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| phase-curvature | Bicontinuous cubic calls show the mean Gaussian curvature in nm⁻². | 2026-11-22 | superseded | 2026-11-23: phase-curvature-legend and phase-curvature-bonnet |
| phase-curvature-legend | Each indexing plot legend line reads phase, lattice, κ in nm⁻² (bicontinuous cubics only) and R², on one line. It drops the reflections-on-peaks count, which both plots already show. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › the indexing plot shows the assigned candidate's fit against the peaks |
| phase-curvature-bonnet | In the Assignment list only a Bonnet group's header shows κ. Single candidate rows and group members show none. | 2026-11-23 | built | packages/HimalayaUI/frontend/test/workspace/CandidateList.test.tsx › leaves the curvature off single candidate rows, cubic or not |
| phase-bonnet-groups | Bonnet-related cubics that share peaks are offered as one group, ticked whole or one member at a time. | 2026-11-22 | built | packages/HimalayaUI/frontend/test/workspace/WorkspacePage.test.tsx › ticks a Bonnet group whole with its number key and undoes it as one edit |

## Export

| id | Decision | Decided | Status | Evidence |
|---|---|---|---|---|
| export-call | JSON and CSV export read the sample call; an uncalled sample exports call_state "uncalled". | 2026-11-22 | built | packages/HimalayaUI/test/test_routes_export.jl › an uncalled sample exports call_state uncalled with its candidates and labels, as JSON and CSV |
| export-call-members | An image's indices hold only the call's members, and its auto_phases keeps the pipeline's own pick for bin/verify. | 2026-11-22 | built | packages/HimalayaUI/test/test_routes_export.jl › an indexed sample's indices hold only its call's members, and auto_phases keeps the auto pick |
| export-candidates | A separate candidates field lists every index, called or not. | 2026-11-22 | built | packages/HimalayaUI/test/test_routes_export.jl › a form-factor sample exports no indices and lists every index as a candidate |
