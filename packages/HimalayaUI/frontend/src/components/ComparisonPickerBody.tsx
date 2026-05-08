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
import type { PickerSampleRow } from "../api";

const COMPARISON_PICKER_BODY_FIXTURE = (
  <div className="flex flex-col flex-1 min-h-0">
    <div className="px-4 py-2 border-b border-border h-[44px]" />
    <div className="flex-1" />
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
  /** Set of exposure ids already in the active draft — rendered as locked. */
  alreadyAddedExposureIds: Set<number>;
  /** Optional ref the parent threads down to focus the search input. */
  searchInputRef?: RefObject<HTMLInputElement>;
  // Filter chip state (search / tag-chip) is internal to the body.
  // PR2 adds an `onPick` prop for immediate-commit shells.
}

export function ComparisonPickerBody({
  experimentId, picks, onPicksChange, alreadyAddedExposureIds,
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
          r.sample.name ?? "", r.sample.label ?? "", r.sample.notes ?? "",
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

  // Pick-set lookup.
  const pickedSampleIds = useMemo(
    () => new Set(picks.map((p) => p.sample_id)),
    [picks],
  );
  const overrideBySampleId = useMemo(() => {
    const m = new Map<number, number>();
    for (const p of picks) {
      if (p.source === "override") m.set(p.sample_id, p.exposure_id);
    }
    return m;
  }, [picks]);

  const togglePickFor = (row: PickerSampleRow, next: boolean): void => {
    if (row.indexing_exposure_id === null) return;
    const pick: Pick = {
      sample_id: row.sample.id,
      exposure_id: row.indexing_exposure_id,
      source: "default",
    };
    if (next) onPicksChange([...picks, pick]);
    else onPicksChange(picks.filter((p) => p.sample_id !== row.sample.id));
  };

  const overridePickFor = (row: PickerSampleRow, exposureId: number): void => {
    const next: Pick = {
      sample_id: row.sample.id,
      exposure_id: exposureId,
      source: exposureId === row.indexing_exposure_id ? "default" : "override",
    };
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
    <Skeleton
      name="comparison-picker-body"
      className="flex flex-col flex-1 min-h-0"
      loading={pickerQ.isLoading}
      fixture={COMPARISON_PICKER_BODY_FIXTURE}
    >
    <div className="flex flex-col flex-1 min-h-0">
      {/* Filters (search + tag chips). */}
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

      <div className="flex-1 min-h-0 overflow-y-auto">
        {recentSamples.length > 0 && (
          <section data-testid="comparison-picker-recents" className="border-b border-border/40 py-2">
            <div className="px-4 text-xs font-medium text-fg-muted uppercase tracking-wide pb-1">
              Recently used
            </div>
            <ul role="listbox" aria-label="Recently used samples" className="flex flex-col">
              {recentSamples.map((r) => {
                const alreadyAdded = r.indexing_exposure_id !== null
                  && alreadyAddedExposureIds.has(r.indexing_exposure_id);
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
                const alreadyAdded = r.indexing_exposure_id !== null
                  && alreadyAddedExposureIds.has(r.indexing_exposure_id);
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
    </div>
    </Skeleton>
  );
}
