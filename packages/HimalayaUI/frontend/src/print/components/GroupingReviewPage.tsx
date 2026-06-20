import { useMemo, useState, type JSX } from "react";
import type { Load } from "../../api";
import { useLoads } from "../../queries";
import { LoadFold } from "./LoadFold";
import { GroupingBulkBar } from "./GroupingBulkBar";
import { SearchInput } from "../ui/SearchInput";
import { SegmentedControl } from "../ui/SegmentedControl";
import { EmptyState } from "../ui/EmptyState";
import { Kicker } from "../ui/Kicker";
import { matchSample } from "../../lib/matchSample";

type Filter = "attn" | "all";

export interface GroupingReviewPageProps {
  experimentId: number;
  onBack: () => void;
  className?: string;
}

export function GroupingReviewPage({ experimentId, onBack, className }: GroupingReviewPageProps): JSX.Element {
  const { data: loads = [], isLoading } = useLoads(experimentId);
  const [filter, setFilter] = useState<Filter>("attn");
  const [search, setSearch] = useState("");
  // ORDERED selection (first-selected = bulk-merge survivor -- Task 15). Membership
  // checks use `.includes`; the LoadFold/SampleFold `selected` prop takes a Set,
  // so derive one. Selection PERSISTS across filter changes (never cleared here).
  const [selection, setSelection] = useState<number[]>([]);
  const selectedSet = useMemo(() => new Set(selection), [selection]);
  const [openLoads, setOpenLoads] = useState<Set<number>>(new Set());
  const [openSamples, setOpenSamples] = useState<Set<number>>(new Set());

  const q = search.trim();
  const flaggedTotal = useMemo(
    () => loads.reduce((a, l) => a + l.samples.filter((s) => s.flag).length, 0),
    [loads],
  );
  const totalSamples = useMemo(() => loads.reduce((a, l) => a + l.samples.length, 0), [loads]);

  // Visible (load, samples) pairs. A searched query overrides the attn filter
  // (search across all loads); a selected sample stays visible even when it
  // would be filtered out, so the persistent selection is verifiable.
  const visible = useMemo(() => {
    const base: Load[] = !q && filter === "attn" ? loads.filter((l) => l.samples.some((s) => s.flag)) : loads;
    return base
      .map((l) => ({
        load: l,
        samples: q
          ? l.samples.filter((s) => matchSample(s, q) || selectedSet.has(s.sample_id))
          : l.samples,
      }))
      .filter((x) => x.samples.length > 0);
  }, [loads, filter, q, selectedSet]);

  const toggleSelect = (id: number) =>
    setSelection((prev) => (prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]));

  const toggleSet = (id: number, set: (u: (p: Set<number>) => Set<number>) => void) =>
    set((p) => { const n = new Set(p); n.has(id) ? n.delete(id) : n.add(id); return n; });

  const noop = () => {};

  return (
    <div className={`mx-auto max-w-[1180px] px-10 pb-32 pt-8${className ? ` ${className}` : ""}`}>
      <button type="button" className="mb-4 inline-flex items-center gap-1.5 text-xs font-semibold text-accent" onClick={onBack}>
        {"←"} Back to corpus
      </button>
      <div className="flex items-end justify-between gap-8">
        <div>
          <Kicker>Grouping review</Kicker>
          <h1 className="text-display text-ink">Check the grouping</h1>
          <p className="mt-2 max-w-[66ch] text-sm text-ink-soft">
            Confirm every sample loaded, has all its exposures, and is split where it should be.
          </p>
        </div>
        <div className="shrink-0 text-right">
          <div className="text-display text-ink"><b className="text-accent">{flaggedTotal}</b> flagged</div>
          <div className="text-xs font-bold uppercase text-ink-faint">of {totalSamples} samples</div>
        </div>
      </div>

      <div className="my-4 flex items-center gap-2">
        <SegmentedControl<Filter>
          value={filter}
          onChange={setFilter}
          options={[{ value: "attn", label: "Needs review" }, { value: "all", label: "All loads" }]}
          aria-label="Grouping filter"
        />
        <div className="ml-auto w-80">
          <SearchInput
            value={search}
            onChange={setSearch}
            ariaLabel="Filter samples"
            placeholder="Filter samples -- name or glob, e.g. HA8* or JC C0?"
          />
        </div>
      </div>

      {isLoading ? null : visible.length === 0 ? (
        <EmptyState title={`No samples match "${q}"`} />
      ) : (
        visible.map(({ load, samples }) => (
          <LoadFold
            key={load.load_id}
            load={load}
            open={openLoads.has(load.load_id) || !!q || load.samples.some((s) => s.flag) || filter === "all"}
            visibleSamples={samples}
            openSamples={openSamples}
            selected={selectedSet}
            onToggleLoad={(id) => toggleSet(id, setOpenLoads)}
            onToggleSampleOpen={(id) => toggleSet(id, setOpenSamples)}
            onToggleSelect={toggleSelect}
            onRename={noop}
            onSplit={noop}
            onMerge={noop}
            onDismissFlag={noop}
            onMoveExposure={noop}
            thumbSrcFor={() => null}
          />
        ))
      )}

      {selection.length > 0 ? (
        <GroupingBulkBar
          count={selection.length}
          noun="sample"
          primaryLabel="Merge"
          primaryEnabled={selection.length >= 2}
          onPrimary={noop}
          onClear={() => setSelection([])}
        />
      ) : null}
    </div>
  );
}
