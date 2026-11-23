# Selecting images to analyze

Decided with Jonathan on 2026-11-23. This replaces the "representative frame"
(the newest kept image) and the combined-by-default curve that was tried
before it.

## The rule

What you select in a sample's image strip is what gets analyzed.

- **One selected image:** the plot, peaks, candidates and auto call are that image's.
- **Two or more:** they are the average of the selected images, a derived
  `exposures` row (`kind = 'combined'`) whose `exposure_sources` list exactly
  those images. It is drawn bold over one faint line per selected image.
- **None selected:** nothing is analyzed. The plot says "Select an image to analyze."

Nothing is ever averaged without the user choosing it. Averaging over images
that differ (radiation damage, a phase change during the series, flare on one
image) measures nothing, so the selection is the user's call, made cheap.

## Selection

| Situation | Behaviour |
| --- | --- |
| A sample's images are first ingested | Every kept image is selected |
| Select all / Select none (strip toolbar) | Every kept image / no image |
| Click a tile | Only looks: Image panel and detector view. The selection doesn't change |
| Checkbox in a tile's top-left corner, or ⌘/Ctrl-click | Toggles that image |
| Shift-click | Adds the range from the last toggled tile |
| Drop an image | It leaves the selection, and a dropped image can't be selected |
| Drop the last kept image | Allowed: nothing is analyzed (a leaked capillary can leave no usable image). The sample row says "All dropped" |
| Restore an image | It comes back unselected, unless nothing else is selected |
| A rescan adds an image | It joins only if every other kept image is selected |

Storage is `exposures.deselected` (0 = in the selection), which is
`selected: boolean` on `GET /api/samples/{id}/exposures`. `PUT
/api/samples/{id}/selection {exposure_ids}` replaces the selection in one
undoable `set_selection` event. `analyzed_exposure_id` on `GET /api/samples`
names the image the selection resolves to, and `analysis_gap` says why there
is none: `all_dropped` or `none_selected` (the sample row shows "All dropped" or
"None selected"). Undo of a drop sends `selected` with the restore, so the image
comes back as it was.

## Peaks and calls follow the set

Each distinct set of selected images gets its own averaged row, which is never
deleted. Going back to a set used before brings back its peaks and edits.
Calls stay with the sample (#350): members on any image or average count.
Export lists the analyzed average and any earlier average that holds a call
member, not every set tried. At boot, `backfill_averages!` builds the average a
selection names but doesn't have yet (a DB from before averages existed).

## Physics

The pooled estimate of intensity per monitor count is ΣIᵢ / ΣTᵢ. Rescaled to
one image (times the mean monitor ΣTᵢ / N) it is ΣIᵢ / N, the plain mean of
raw counts, with σ = √(Σσᵢ²) / N. That is unbiased under beam decay.
Normalising each image first and taking an unweighted mean would not be. The
average's monitors are the mean of the sources' flux, so subtraction scales it
like one image. Images in one average must share a q grid.

## Migration

`exposures_deselected_v1` adds the column. A sample with phase members or
peak edits on an image starts with only the newest such image selected, so
reviewed samples open where they were left. Others start fully selected.

## Follow-ups

- Images flagged beam-off start deselected (exception to "every kept image",
  waiting on PR 517's `beam_off` flag).
- A disagreement flag on tiles whose image differs from the rest of the
  selection. Nothing is deselected automatically.
- Peak-model training from reviewed averages.
