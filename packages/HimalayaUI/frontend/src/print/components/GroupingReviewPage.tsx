import { useMemo, useState, type JSX } from "react";
import type { Load, LoadSample } from "../../api";
import {
  useLoads, useExperiment,
  useMergeSamples, useRenameSample, useMoveExposure,
  useSplitSample, useDismissGroupingFlag, useUndoDismissGroupingFlag,
} from "../../queries";
import { useAppState } from "../../state";
import { LoadFold } from "./LoadFold";
import { GroupingBulkBar } from "./GroupingBulkBar";
import { SearchInput } from "../ui/SearchInput";
import { SegmentedControl } from "../ui/SegmentedControl";
import { EmptyState } from "../ui/EmptyState";
import { Kicker } from "../ui/Kicker";
import { ProgressBar } from "../ui/ProgressBar";
import { ModalShell } from "../ui/ModalShell";
import { Button } from "../ui/Button";
import { Menu } from "../ui/Menu";
import { matchSample } from "../../lib/matchSample";
import { effectiveIngestStatus } from "../../lib/ingestStatus";
import { showToast } from "../../lib/toast";
import { useUndoStack } from "../../hooks/useUndoStack";

type Filter = "attn" | "all";

interface ConfirmState {
  kind: "merge";
  loserId: number;
  survivorId: number;
  survivorLabel: string;
}

interface BulkMergeConfirmState {
  survivorId: number;
  survivorLabel: string;
  loserIds: number[];
}

/** State for the Move picker Menu: which exposure is being moved and which
 *  same-load sibling samples are available as destinations. */
interface MovePickerState {
  exposureId: number;
  /** Current owner (excluded from the destination list). */
  fromSampleId: number;
  /** All sibling samples in the same load (current owner excluded). */
  siblings: LoadSample[];
  /** The Move button element — used to position the relative wrapper. */
  anchorEl: HTMLElement;
}

interface UndoEntry {
  label: string;
  undo: () => void;
}

export interface GroupingReviewPageProps {
  experimentId: number;
  onBack: () => void;
  /** "Confirm groups" target: where to go once the scan finishes and the user
   *  accepts the grouping. Defaults to onBack (the corpus). */
  onConfirm?: () => void;
  className?: string;
}

export function GroupingReviewPage({ experimentId, onBack, onConfirm, className }: GroupingReviewPageProps): JSX.Element {
  const { data: loads = [], isLoading } = useLoads(experimentId);
  const exp = useExperiment(experimentId);
  const inFlight = useAppState((s) => s.ingestInFlight?.[experimentId]);
  // Scanning = an initial-scan frame is in flight (the combined scan + review
  // surface, p1-grouping). Loads unfold live as the SSE invalidates the query.
  // Reconcile the live SSE overlay with the persisted resting truth: a terminal
  // persisted state (complete/failed) overrides a stale "scanning" overlay left
  // behind by a dropped `ingest_complete` frame (8c). With useExperiment's
  // scoped refetchInterval refreshing ingest_status while a scan is active, the
  // surface self-heals to the post-scan review without a manual reload.
  const scanning =
    effectiveIngestStatus(inFlight?.status, exp.data?.ingest_status) === "scanning";
  const [filter, setFilter] = useState<Filter>("attn");
  const [search, setSearch] = useState("");
  // ORDERED selection (first-selected = bulk-merge survivor -- Task 15). Membership
  // checks use `.includes`; the LoadFold/SampleFold `selected` prop takes a Set,
  // so derive one. Selection PERSISTS across filter changes (never cleared here).
  const [selection, setSelection] = useState<number[]>([]);
  const selectedSet = useMemo(() => new Set(selection), [selection]);
  const [openLoads, setOpenLoads] = useState<Set<number>>(new Set());
  const [openSamples, setOpenSamples] = useState<Set<number>>(new Set());

  // Confirm modal state (merge only for now; split confirm added if needed)
  const [confirm, setConfirm] = useState<ConfirmState | null>(null);
  // Bulk-merge confirm modal state
  const [bulkMergeConfirm, setBulkMergeConfirm] = useState<BulkMergeConfirmState | null>(null);
  // Move picker state: open when an exposure's Move button is activated.
  const [movePicker, setMovePicker] = useState<MovePickerState | null>(null);

  // Mutation hooks — one instance per experiment, entity ids ride in mutate() input
  const { mutate: mergeMutate } = useMergeSamples(experimentId);
  const { mutate: renameMutate } = useRenameSample(experimentId);
  const { mutate: moveMutate } = useMoveExposure(experimentId);
  const { mutate: splitMutate } = useSplitSample(experimentId);
  const { mutate: dismissMutate } = useDismissGroupingFlag(experimentId);
  const { mutate: undoDismissMutate } = useUndoDismissGroupingFlag(experimentId);

  // Undo stack for reversible single-entity ops (rename / move)
  const undoStack = useUndoStack<UndoEntry>();

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
    // During scanning the surface shows EVERY load as it lands (clean ones
    // collapsed); the attn filter only narrows the post-scan standalone review.
    const base: Load[] = !q && filter === "attn" && !scanning
      ? loads.filter((l) => l.samples.some((s) => s.flag))
      : loads;
    return base
      .map((l) => ({
        load: l,
        samples: q
          ? l.samples.filter((s) => matchSample(s, q) || selectedSet.has(s.sample_id))
          : l.samples,
      }))
      .filter((x) => x.samples.length > 0);
  }, [loads, filter, q, selectedSet, scanning]);

  const toggleSelect = (id: number) =>
    setSelection((prev) => (prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]));

  const toggleSet = (id: number, set: (u: (p: Set<number>) => Set<number>) => void) =>
    set((p) => { const n = new Set(p); n.has(id) ? n.delete(id) : n.add(id); return n; });

  // --- Edit callbacks ---

  const handleRename = (sampleId: number, newName: string) => {
    // Find current name for undo
    let prevName: string | undefined;
    for (const l of loads) {
      const s = l.samples.find((x) => x.sample_id === sampleId);
      if (s) { prevName = s.name; break; }
    }
    renameMutate({ sampleId, name: newName });
    if (prevName !== undefined) {
      const prev = prevName;
      undoStack.push({
        label: "rename",
        undo: () => renameMutate({ sampleId, name: prev }),
      });
    }
    showToast("Renamed", "info", {
      label: "Undo",
      onClick: () => { const e = undoStack.pop(); if (e) e.undo(); },
    });
  };

  const handleSplit = (sampleId: number) => {
    // Find the sample and its split flag
    let sample: LoadSample | undefined;
    for (const l of loads) {
      const s = l.samples.find((x) => x.sample_id === sampleId);
      if (s) { sample = s; break; }
    }
    if (!sample) return;
    const flag = sample.flag?.kind === "split" ? sample.flag : null;
    const splitIdx = flag?.split_at_index ?? Math.floor(sample.exposures.length / 2);
    const splitExposureIds = sample.exposures.slice(splitIdx).map((e) => e.id);
    // Use the sample name with a suffix for the split-off portion
    const splitName = `${sample.name} (split)`;
    splitMutate({ sampleId, exposureIds: splitExposureIds, name: splitName });
    showToast("Split sample", "info");
  };

  const handleMerge = (loserId: number, survivorId: number) => {
    // Find the survivor label for the confirm dialog
    let survivorLabel = String(survivorId);
    for (const l of loads) {
      const s = l.samples.find((x) => x.sample_id === survivorId);
      if (s) { survivorLabel = s.name; break; }
    }
    setConfirm({ kind: "merge", loserId, survivorId, survivorLabel });
  };

  const handleConfirmMerge = () => {
    if (!confirm || confirm.kind !== "merge") return;
    mergeMutate({ loserId: confirm.loserId, survivorId: confirm.survivorId });
    showToast(`Merged into ${confirm.survivorLabel}`, "info");
    // Merge is multi-row — no undo action in v1 (spec §9.3). Flag to Jonathan.
    setConfirm(null);
  };

  const handleDismissFlag = (sampleId: number) => {
    // Capture the original flag for the undo optimistic restore
    let originalFlag: LoadSample["flag"] = null;
    for (const l of loads) {
      const s = l.samples.find((x) => x.sample_id === sampleId);
      if (s) { originalFlag = s.flag; break; }
    }
    if (!originalFlag) return;

    dismissMutate({
      sampleId,
      flagKind: originalFlag.kind,
      ...(originalFlag.kind === "merge" ? { mergeWithSampleId: originalFlag.merge_with_sample_id } : {}),
    });

    // Undo: re-show the flag via the inverse backend route
    showToast("Flag dismissed", "info", {
      label: "Undo",
      onClick: () => {
        undoDismissMutate({ sampleId, originalFlag });
      },
    });
  };

  const handleMoveExposure = (sampleId: number, exposureId: number, anchorEl: HTMLElement) => {
    // Find the load that owns this sample and collect its siblings (every
    // sample in the same load except the current owner).
    let siblings: LoadSample[] = [];
    for (const l of loads) {
      if (l.samples.some((s) => s.sample_id === sampleId)) {
        siblings = l.samples.filter((s) => s.sample_id !== sampleId);
        break;
      }
    }
    setMovePicker({ exposureId, fromSampleId: sampleId, siblings, anchorEl });
  };

  const handlePickDestination = (destSampleId: number) => {
    if (!movePicker) return;
    const { exposureId, fromSampleId } = movePicker;
    setMovePicker(null);
    moveMutate({ exposureId, sampleId: destSampleId });
    undoStack.push({
      label: "move",
      undo: () => moveMutate({ exposureId, sampleId: fromSampleId }),
    });
    showToast("Moved exposure", "info", {
      label: "Undo",
      onClick: () => { const e = undoStack.pop(); if (e) e.undo(); },
    });
  };

  const openBulkMergeConfirm = () => {
    if (selection.length < 2) return;
    const [survivorId, ...loserIds] = selection;
    // Find survivor label
    let survivorLabel = String(survivorId);
    for (const l of loads) {
      const s = l.samples.find((x) => x.sample_id === survivorId);
      if (s) { survivorLabel = s.name; break; }
    }
    setBulkMergeConfirm({ survivorId, survivorLabel, loserIds });
  };

  const handleConfirmBulkMerge = () => {
    if (!bulkMergeConfirm) return;
    const { survivorId, survivorLabel, loserIds } = bulkMergeConfirm;
    loserIds.forEach((loserId) => mergeMutate({ loserId, survivorId }));
    setSelection([]);
    setBulkMergeConfirm(null);
    showToast(`Merged ${loserIds.length + 1} samples into ${survivorLabel}`, "info");
  };

  return (
    <div className={`mx-auto max-w-[1180px] px-10 pb-32 pt-8${className ? ` ${className}` : ""}`}>
      {!scanning && (
        <button type="button" className="mb-4 inline-flex items-center gap-1.5 text-xs font-semibold text-accent" onClick={onBack}>
          {"←"} Back to corpus
        </button>
      )}

      {scanning ? (
        // Combined scan + grouping-review header (p1-grouping): live progress as
        // loads land. processed/total ride the SSE ingestInFlight; loads/samples
        // come from the live-invalidated query.
        <div data-testid="grouping-scanning-header">
          <Kicker tone="accent">Scanning · review groups as they land</Kicker>
          <h1 className="text-display text-ink">{exp.data?.name ?? `Experiment ${experimentId}`}</h1>
          <div className="mt-3 flex items-center gap-4">
            <span className="text-meta text-ink-soft whitespace-nowrap">
              Parsing exposures… {inFlight?.processed ?? 0} / ~{inFlight?.total ?? 0}
              {" · "}{loads.length} loads · {totalSamples} samples
            </span>
            <div className="flex-1 min-w-[8rem]">
              <ProgressBar
                value={inFlight?.processed ?? 0}
                total={Math.max(inFlight?.total ?? 1, 1)}
                label="Scan progress"
              />
            </div>
            {flaggedTotal > 0 && (
              <span className="text-meta font-semibold text-print-accent whitespace-nowrap" data-testid="grouping-flag-count">
                {flaggedTotal} {flaggedTotal === 1 ? "flag" : "flags"} to review
              </span>
            )}
          </div>
        </div>
      ) : (
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
      )}

      {/* Filter/search is for the post-scan standalone review; during a live
          scan the surface just shows every load as it lands. */}
      <div className={`my-4 flex items-center gap-2${scanning ? " hidden" : ""}`}>
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
            placeholder="Filter samples by name or glob, e.g. HA8* or JC C0?"
          />
        </div>
      </div>

      {isLoading ? null : visible.length === 0 ? (
        // Three distinct empty cases. Collapsing them all to `No samples match
        // ""` (the old bug) read as an error even in the all-clear case.
        q ? (
          <EmptyState title={`No samples match "${q}"`} />
        ) : totalSamples === 0 ? (
          <EmptyState
            title="No loads yet"
            body="Rescan this experiment to group its exposures into samples."
          />
        ) : filter === "attn" ? (
          <EmptyState
            title="Nothing needs review"
            body="Every sample looks correctly grouped. Switch to All loads to see them all."
          />
        ) : (
          <EmptyState title="No samples in this experiment" />
        )
      ) : (
        visible.map(({ load, samples }) => (
          <LoadFold
            key={load.load_id}
            load={load}
            open={openLoads.has(load.load_id) || !!q || load.samples.some((s) => s.flag) || (filter === "all" && !scanning)}
            visibleSamples={samples}
            openSamples={openSamples}
            selected={selectedSet}
            onToggleLoad={(id) => toggleSet(id, setOpenLoads)}
            onToggleSampleOpen={(id) => toggleSet(id, setOpenSamples)}
            onToggleSelect={toggleSelect}
            onRename={handleRename}
            onSplit={handleSplit}
            onMerge={handleMerge}
            onDismissFlag={handleDismissFlag}
            onMoveExposure={handleMoveExposure}
            thumbSrcFor={() => null}
          />
        ))
      )}

      {/* Still-landing loads: a faint placeholder while the scan parses more. */}
      {scanning && (inFlight?.processed ?? 0) < (inFlight?.total ?? 0) && (
        <div
          data-testid="grouping-unfolding"
          className="mt-3 rounded-md border border-dashed border-hair px-5 py-4 text-meta text-ink-faint"
        >
          unfolding…
        </div>
      )}

      {/* Merge confirm modal (single-entity) */}
      <ModalShell
        open={confirm?.kind === "merge"}
        onClose={() => setConfirm(null)}
        size="sm"
        aria-label="Confirm merge"
        testId="merge-confirm"
      >
        <div className="p-6">
          <p className="text-sm text-ink">
            Merge this sample into <b className="text-ink">{confirm?.survivorLabel}</b>?
            All exposures will be moved to the surviving sample. This cannot be undone automatically.
          </p>
          <div className="mt-4 flex justify-end gap-2">
            <Button variant="ghost" onClick={() => setConfirm(null)}>Cancel</Button>
            <Button variant="outline" onClick={handleConfirmMerge}>Merge</Button>
          </div>
        </div>
      </ModalShell>

      {/* Bulk-merge confirm modal */}
      <ModalShell
        open={bulkMergeConfirm !== null}
        onClose={() => setBulkMergeConfirm(null)}
        size="sm"
        aria-label="Confirm bulk merge"
        testId="bulk-merge-confirm"
      >
        <div className="p-6">
          <p className="text-sm text-ink">
            Merge {bulkMergeConfirm ? bulkMergeConfirm.loserIds.length + 1 : 0} samples
            into <b className="text-ink">{bulkMergeConfirm?.survivorLabel}</b>?
            All exposures will be moved to the surviving sample. This cannot be undone automatically.
          </p>
          <div className="mt-4 flex justify-end gap-2">
            <Button variant="ghost" onClick={() => setBulkMergeConfirm(null)}>Cancel</Button>
            <Button variant="outline" onClick={handleConfirmBulkMerge}>Merge</Button>
          </div>
        </div>
      </ModalShell>

      {/* Move picker: a floating Menu of same-load sibling samples. Rendered
          as a fixed overlay anchored to the trigger element's bounding rect so
          the absolute Menu finds its relative parent without disrupting layout.
          Close on Escape (Menu handles this) or outside-click (overlay backdrop). */}
      {movePicker !== null ? (() => {
        const rect = movePicker.anchorEl.getBoundingClientRect();
        const pickerOptions = movePicker.siblings.map((s) => ({
          value: String(s.sample_id),
          label: s.name,
        }));
        return (
          <>
            {/* Invisible backdrop to catch outside-clicks. Above the Confirm
                footer (z-30) so a picker opened near the bottom isn't covered. */}
            <div
              aria-hidden="true"
              style={{ position: "fixed", inset: 0, zIndex: 39 }}
              onClick={() => setMovePicker(null)}
            />
            {/* Anchor wrapper: positioned at the trigger element */}
            <div
              style={{
                position: "fixed",
                top: rect.bottom,
                left: rect.left,
                zIndex: 40,
              }}
            >
              <div className="relative">
                <Menu<string>
                  open
                  options={pickerOptions}
                  onSelect={(value) => handlePickDestination(Number(value))}
                  onClose={() => setMovePicker(null)}
                  aria-label="Move exposure to sample"
                  className="right-0 top-0"
                />
              </div>
            </div>
          </>
        );
      })() : null}

      {selection.length > 0 && !bulkMergeConfirm ? (
        <GroupingBulkBar
          count={selection.length}
          noun="sample"
          primaryLabel="Merge"
          primaryEnabled={selection.length >= 2}
          onPrimary={openBulkMergeConfirm}
          onClear={() => setSelection([])}
        />
      ) : null}

      {/* Confirm-groups footer (p1-grouping): the surface's exit. Disabled while
          the scan is in flight — settled loads are reviewable immediately, but
          Confirm waits for the scan to finish (later loads can still raise flags).
          Ingest is additive, so confirming never re-touches a settled load. */}
      <footer
        data-testid="grouping-footer"
        className="fixed bottom-0 left-0 right-0 z-30 flex items-center justify-between gap-4 border-t border-hair bg-paper px-8 py-3"
      >
        <span className="text-meta text-ink-soft">
          {scanning
            ? "Review flags as loads land. Confirm unlocks when the scan finishes. Later loads can still raise flags."
            : "Confirm the grouping to head to the corpus. You can always return to review."}
        </span>
        <Button
          variant="accent"
          data-testid="grouping-confirm"
          disabled={scanning}
          onClick={() => (onConfirm ?? onBack)()}
        >
          {scanning ? "Confirm groups · scanning…" : "Confirm groups"}
        </Button>
      </footer>
    </div>
  );
}
