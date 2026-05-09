import { useMemo, useRef, useState } from "react";
import type { RefObject } from "react";
import { Skeleton } from "boneyard-js/react";
import {
  useRecentlyPickedExposures,
  useSampleTags,
  usePickerSamples,
} from "../queries";
import { useCurrentUserId } from "../hooks/useCurrentUserId";
import { SamplePickerRow } from "./SamplePickerRow";
import { HintText } from "./ui";
import type { PickerSampleRow } from "../api";

const FIXTURE_ROWS: PickerSampleRow[] = [
  {
    sample: { id: 1, experiment_id: 0, name: "Sample A", display_name: null,
              notes: "Cubic phase replicate · 24°C", tags: [] },
    indexing_exposure_id: 11,
    all_exposures: [{ id: 11, sample_id: 1, filename: "JC068P_E257_S1418_tot", selected: true }],
  },
  {
    sample: { id: 2, experiment_id: 0, name: "Sample B", display_name: null,
              notes: "Hexagonal phase reference · 37°C", tags: [] },
    indexing_exposure_id: 21,
    all_exposures: [{ id: 21, sample_id: 2, filename: "JC068P_E258_S1418_tot", selected: true }],
  },
  {
    sample: { id: 3, experiment_id: 0, name: "Form factor", display_name: null, notes: null, tags: [] },
    indexing_exposure_id: 31,
    all_exposures: [{ id: 31, sample_id: 3, filename: "JC068P_E259_S1418_tot", selected: true }],
  },
];

// Skeleton fixture covers ONLY the result list area — the filter strip
// (search input + tag chips) renders outside the Skeleton so it's
// interactive immediately on cold load (regression fix from PR #96 review:
// search input was unfocusable when ref pointed inside the Skeleton's
// fixture/fallback DOM during cold fetch).
const COMPARISON_PICKER_BODY_FIXTURE = (
  <div className="flex-1 overflow-hidden">
    <div className="py-2">
      <div className="px-4 text-xs font-medium text-fg-muted uppercase tracking-wide pb-1">
        All samples
      </div>
      <ul className="flex flex-col">
        {FIXTURE_ROWS.map((r) => (
          <li key={r.sample.id}>
            <SamplePickerRow
              row={r}
              checked={false}
              onCheckedChange={() => {}}
              overrideExposureId={undefined}
              onOverrideChange={() => {}}
            />
          </li>
        ))}
      </ul>
    </div>
  </div>
);

export type Pick = {
  sample_id: number;
  exposure_id: number;
  source: "default" | "override";
};

interface Props {
  experimentId: number | undefined;
  /** Controlled picks list. */
  picks: Pick[];
  onPicksChange: (next: Pick[]) => void;
  /**
   * Immediate-mode callback. When set, the body fires `onPick` on each toggle
   * and does NOT mutate the `picks`/`onPicksChange` controlled state — the
   * caller is responsible for committing each pick directly (e.g. `addMember`).
   * `picks` and `onPicksChange` are still required for checked-state display in
   * the modal shell; pass `[]` + a no-op when using immediate mode.
   */
  onPick?: (pick: Pick) => void;
  /** Set of exposure ids already in the active draft — rendered as locked. */
  alreadyAddedExposureIds: Set<number>;
  /** Optional ref the parent threads down to focus the search input. */
  searchInputRef?: RefObject<HTMLInputElement>;
  // Filter chip state (search / tag-chip) is internal to the body.
}

export function ComparisonPickerBody({
  experimentId, picks, onPicksChange, onPick, alreadyAddedExposureIds,
  searchInputRef,
}: Props): JSX.Element {
  const userId = useCurrentUserId();

  // Data sources.
  const pickerQ  = usePickerSamples(experimentId);
  const recentsQ = useRecentlyPickedExposures(userId, 20);
  const tagsQ    = useSampleTags(experimentId);
  // `useExperiments` is intentionally not invoked — the picker scope is the
  // active experiment only in PR1. Cross-experiment widening can return in
  // a follow-on if user demand surfaces.

  const fallbackRef = useRef<HTMLInputElement>(null);
  const inputRef = searchInputRef ?? fallbackRef;

  const [search, setSearch]             = useState("");
  const [selectedTags, setSelectedTags] = useState<Set<string>>(new Set());

  // Build the rendered row list (filtered + sorted) from the picker-samples query.
  const allRows = pickerQ.data ?? [];
  const filteredRows = useMemo<PickerSampleRow[]>(() => {
    return allRows.filter((r) => {
      if (selectedTags.size > 0) {
        const match = r.sample.tags.some((t) => selectedTags.has(`${t.key}:${t.value}`));
        if (!match) return false;
      }
      if (search.trim() !== "") {
        const needle = search.toLowerCase();
        const haystack = [
          r.sample.name ?? "", r.sample.display_name ?? "", r.sample.notes ?? "",
          ...r.sample.tags.map((t) => t.value),
          ...r.all_exposures.map((e) => e.filename ?? ""),
        ].map((s) => s.toLowerCase());
        if (!haystack.some((s) => s.includes(needle))) return false;
      }
      return true;
    }).sort((a, b) =>
      (a.sample.name ?? "").localeCompare(b.sample.name ?? ""),
    );
  }, [allRows, selectedTags, search]);

  // Recents → derive sample-id ordering, dedupe to one row per sample, then
  // exclude from main list. useMemo keeps server state out of Zustand.
  const recentSamples = useMemo<PickerSampleRow[]>(() => {
    const recentIds = recentsQ.data ?? [];
    if (recentIds.length === 0) return [];
    const sampleByExposureId = new Map<number, PickerSampleRow>();
    for (const r of allRows) {
      for (const e of r.all_exposures) sampleByExposureId.set(e.id, r);
    }
    const seen = new Set<number>();
    const out: PickerSampleRow[] = [];
    for (const eid of recentIds) {
      const r = sampleByExposureId.get(eid);
      if (!r || seen.has(r.sample.id)) continue;
      seen.add(r.sample.id);
      out.push(r);
    }
    return out;
  }, [recentsQ.data, allRows]);

  const recentSampleIds = useMemo(
    () => new Set(recentSamples.map((r) => r.sample.id)),
    [recentSamples],
  );
  const mainListRows = useMemo(
    () => filteredRows.filter((r) => !recentSampleIds.has(r.sample.id)),
    [filteredRows, recentSampleIds],
  );

  // Pick-set lookup. In immediate mode (`onPick` set), derive from draft
  // state via alreadyAddedExposureIds × picker data so SamplePickerRow's
  // checkbox + override radios reflect the actual committed state. In modal
  // mode, derive from the controlled `picks` list.
  const pickedSampleIds = useMemo(() => {
    if (onPick) {
      const out = new Set<number>();
      for (const r of allRows) {
        if (r.all_exposures.some((e) => alreadyAddedExposureIds.has(e.id))) {
          out.add(r.sample.id);
        }
      }
      return out;
    }
    return new Set(picks.map((p) => p.sample_id));
  }, [onPick, picks, allRows, alreadyAddedExposureIds]);

  const overrideBySampleId = useMemo(() => {
    if (onPick) {
      // Immediate mode: a member exists in the draft for some exposure of
      // this sample. If that exposure is NOT the indexing one, surface as an
      // override so the corresponding radio shows checked.
      const m = new Map<number, number>();
      for (const r of allRows) {
        const matching = r.all_exposures.find((e) => alreadyAddedExposureIds.has(e.id));
        if (matching && matching.id !== r.indexing_exposure_id) {
          m.set(r.sample.id, matching.id);
        }
      }
      return m;
    }
    const m = new Map<number, number>();
    for (const p of picks) {
      if (p.source === "override") m.set(p.sample_id, p.exposure_id);
    }
    return m;
  }, [onPick, picks, allRows, alreadyAddedExposureIds]);

  const togglePickFor = (row: PickerSampleRow, next: boolean): void => {
    if (row.indexing_exposure_id === null) return;
    const pick: Pick = {
      sample_id: row.sample.id,
      exposure_id: row.indexing_exposure_id,
      source: "default",
    };
    if (onPick) {
      if (next) onPick(pick);
      return;
    }
    if (next) onPicksChange([...picks, pick]);
    else onPicksChange(picks.filter((p) => p.sample_id !== row.sample.id));
  };

  const overridePickFor = (row: PickerSampleRow, exposureId: number): void => {
    const next: Pick = {
      sample_id: row.sample.id,
      exposure_id: exposureId,
      source: exposureId === row.indexing_exposure_id ? "default" : "override",
    };
    if (onPick) { onPick(next); return; }
    const i = picks.findIndex((p) => p.sample_id === row.sample.id);
    if (i < 0) onPicksChange([...picks, next]);
    else onPicksChange(picks.map((p, j) => (j === i ? next : p)));
  };

  // Tag filter helper.
  const toggleTag = (key: string, value: string): void => {
    const id = `${key}:${value}`;
    setSelectedTags((prev) => {
      const m = new Set(prev);
      if (m.has(id)) m.delete(id);
      else m.add(id);
      return m;
    });
  };

  return (
    <div className="flex flex-col flex-1 min-h-0">
      {/* Filters (search + tag chips). Rendered OUTSIDE the Skeleton so the
          search input is mounted on first paint — keeps `inputRef.current?.focus()`
          working during cold load (PR #96 review fix). */}
      <div className="px-4 py-2 border-b border-border space-y-2">
        <input
          ref={inputRef}
          data-testid="comparison-picker-search"
          type="search"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          placeholder="Search samples, notes, filenames…"
          className="w-full bg-transparent border border-border rounded px-2 py-1 text-sm"
          spellCheck={false}
        />
        {(tagsQ.data ?? []).length > 0 && (
          <div data-testid="picker-filter-tag" className="flex flex-wrap gap-1 items-center">
            <span className="text-xs text-fg-muted">Tags:</span>
            {(tagsQ.data ?? []).map((t) => {
              const id = `${t.key}:${t.value}`;
              const active = selectedTags.has(id);
              return (
                <button
                  key={id}
                  type="button"
                  data-testid={`picker-tag-option-${id}`}
                  data-active={active ? "true" : undefined}
                  onClick={() => toggleTag(t.key, t.value)}
                  className={
                    "text-xs px-2 py-0.5 rounded-full border " +
                    (active
                      ? "bg-accent/15 border-accent text-accent"
                      : "border-border text-fg-muted hover:text-fg")
                  }
                >
                  {t.key}:{t.value}
                </button>
              );
            })}
          </div>
        )}
      </div>

      <Skeleton
        name="comparison-picker-body"
        className="flex-1 min-h-0 flex flex-col"
        loading={pickerQ.isLoading}
        fixture={COMPARISON_PICKER_BODY_FIXTURE}
        fallback={<div className="flex-1 flex items-center justify-center"><HintText>Loading samples…</HintText></div>}
      >
      <div className="flex-1 min-h-0 overflow-y-auto">
        {recentSamples.length > 0 && (
          <section data-testid="comparison-picker-recents" className="border-b border-border/40 py-2">
            <div className="px-4 text-xs font-medium text-fg-muted uppercase tracking-wide pb-1">
              Recently used
            </div>
            <ul role="listbox" aria-label="Recently used samples" className="flex flex-col">
              {recentSamples.map((r) => {
                // alreadyAdded = "any exposure of this sample is in the draft"
                // — the indexing-only check missed override-added rows in
                // immediate mode (PR #97 review): clicking the checkbox to
                // un-check an override-added row was a silent no-op since
                // togglePickFor(_, false) does nothing in onPick mode.
                // Locking via the some() check covers both modes.
                const alreadyAdded = r.all_exposures.some(
                  (e) => alreadyAddedExposureIds.has(e.id),
                );
                return (
                  <li key={`recent-${r.sample.id}`} data-testid="picker-row">
                    <SamplePickerRow
                      row={r}
                      checked={pickedSampleIds.has(r.sample.id)}
                      onCheckedChange={(next) => togglePickFor(r, next)}
                      overrideExposureId={overrideBySampleId.get(r.sample.id)}
                      onOverrideChange={(eid) => overridePickFor(r, eid)}
                      alreadyAdded={alreadyAdded}
                    />
                  </li>
                );
              })}
            </ul>
          </section>
        )}

        {mainListRows.length === 0 ? (
          <div data-testid="comparison-picker-empty"
               className="px-4 py-8 text-center text-fg-muted text-sm">
            No samples match. Try clearing filters.
          </div>
        ) : (
          <section className="py-2">
            <div className="px-4 text-xs font-medium text-fg-muted uppercase tracking-wide pb-1">
              All samples
            </div>
            <ul role="listbox" aria-label="Samples" className="flex flex-col">
              {mainListRows.map((r) => {
                // alreadyAdded = "any exposure of this sample is in the draft"
                // — the indexing-only check missed override-added rows in
                // immediate mode (PR #97 review): clicking the checkbox to
                // un-check an override-added row was a silent no-op since
                // togglePickFor(_, false) does nothing in onPick mode.
                // Locking via the some() check covers both modes.
                const alreadyAdded = r.all_exposures.some(
                  (e) => alreadyAddedExposureIds.has(e.id),
                );
                return (
                  <li key={r.sample.id} data-testid="picker-row">
                    <SamplePickerRow
                      row={r}
                      checked={pickedSampleIds.has(r.sample.id)}
                      onCheckedChange={(next) => togglePickFor(r, next)}
                      overrideExposureId={overrideBySampleId.get(r.sample.id)}
                      onOverrideChange={(eid) => overridePickFor(r, eid)}
                      alreadyAdded={alreadyAdded}
                    />
                  </li>
                );
              })}
            </ul>
          </section>
        )}
      </div>
      </Skeleton>
    </div>
  );
}
