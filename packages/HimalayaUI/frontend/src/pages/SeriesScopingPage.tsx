import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Link, useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import type { Trace } from "../api";
import {
  useCorpusSampleTags,
  useCorpusPickerSamples,
  useScopeSeries,
  useMemberTraces,
  useMemberIndices,
} from "../queries";
import { proposeOrdering, type OrderingRow } from "../lib/scoping/proposeOrdering";
import { splitProposal, humanizeKey } from "../lib/scoping/splitProposal";
import { parseSortKey } from "../lib/scoping/parseSortKey";
import { dominantPhase } from "../lib/scoping/dominantPhase";
import { ScopingAutogroupCard } from "../components/ScopingAutogroupCard";
import { ScopingOrderField } from "../components/ScopingOrderField";
import { ScopingRow } from "../components/ScopingRow";
import { ScopingLooseMatches } from "../components/ScopingLooseMatches";
import { Card } from "../components/ui";
import { PhaseStrip } from "../components/ui/PhaseStrip";
import { Kicker } from "../components/ui/Kicker";
import { ScopingFoot } from "../components/ScopingFoot";
import { ScopingConfirmModal } from "../components/ScopingConfirmModal";

// Static skeleton shape for boneyard's headless capture (docs/boneyard.md).
const SCOPING_FIXTURE = (
  <div className="space-y-2">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="flex items-center gap-3 rounded border border-hair p-3">
        <div className="h-4 w-4 rounded bg-paper-sunk" />
        <div className="h-4 w-1/3 rounded bg-paper-sunk" />
        <div className="h-4 w-1/4 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

interface UndoEntry {
  rows: OrderingRow[];
  loose: OrderingRow[];
  label: string;
}

/**
 * SeriesScopingPage — the machine-proposes / human-confirms scoping worksheet
 * at /series/new (#174 / I3.4 / R7 #230). Per master plan §6.4 this is the
 * confirm-and-build GATE that *writes* the structured (key,value) sample_tags —
 * NOT series creation.
 *
 * It surfaces `proposeOrdering` as a worksheet plate (series-scoping.html): an
 * autogroup summary, a parsed "Ordered by" field, member rows with trace
 * sparklines + confirm-not-fill-out values + amber "check the read" flags, a
 * "Himalaya also found" loose-match section (samples lacking the ordering key),
 * a preview phase strip, and a narrative gate footer. Nothing is durable until
 * Confirm & build; every change steps back in-session via Undo / ⌘Z. On confirm
 * it writes one batch (source='scoping') through the queue, then navigates to
 * the folio.
 */
export function SeriesScopingPage(): JSX.Element {
  const navigate = useNavigate();
  const tagsQ = useCorpusSampleTags();
  const pickerQ = useCorpusPickerSamples();
  const scopeSeries = useScopeSeries();

  const isLoading = tagsQ.isLoading || pickerQ.isLoading;
  const isError = tagsQ.isError || pickerQ.isError;

  const proposal = useMemo(
    () => proposeOrdering(tagsQ.data ?? [], pickerQ.data ?? []),
    [tagsQ.data, pickerQ.data],
  );
  const split = useMemo(() => splitProposal(proposal), [proposal]);

  // Local editable copies — seeded from the split once queries resolve, then
  // user-owned. `rows` = members; `loose` = Himalaya-also-found.
  const [rows, setRows] = useState<OrderingRow[]>([]);
  const [loose, setLoose] = useState<OrderingRow[]>([]);
  useEffect(() => {
    setRows(split.members);
    setLoose(split.looseMatches);
  }, [split.members, split.looseMatches]);

  // In-session undo stack (S-F). Each mutating action snapshots prior state.
  const [history, setHistory] = useState<UndoEntry[]>([]);
  const snapshot = useCallback(
    (label: string) => setHistory((h) => [...h, { rows, loose, label }]),
    [rows, loose],
  );
  const undo = useCallback(() => {
    setHistory((h) => {
      if (h.length === 0) return h;
      const last = h[h.length - 1]!;
      setRows(last.rows);
      setLoose(last.loose);
      return h.slice(0, -1);
    });
  }, []);

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        undo();
      }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [undo]);

  // Member ordering: render low→high by parsed sort key; unparseable last,
  // stable by sampleId. Pure surfacing — does NOT reorder the batch payload.
  const orderedRows = useMemo(() => {
    return [...rows].sort((a, b) => {
      const ka = parseSortKey(a.value), kb = parseSortKey(b.value);
      if (ka === null && kb === null) return a.sampleId - b.sampleId;
      if (ka === null) return 1;
      if (kb === null) return -1;
      return ka - kb || a.sampleId - b.sampleId;
    });
  }, [rows]);

  // Per-member trace + phase fetches keyed on the indexing exposure id.
  const pickerById = useMemo(() => {
    const m = new Map<number, number | null>();
    for (const s of pickerQ.data ?? []) m.set(s.sample.id, s.indexing_exposure_id);
    return m;
  }, [pickerQ.data]);
  const allSampleIds = useMemo(() => [...rows, ...loose].map((r) => r.sampleId), [rows, loose]);
  const exposureIds = useMemo(
    () => allSampleIds.map((id) => pickerById.get(id)).filter((e): e is number => e != null),
    [allSampleIds, pickerById],
  );
  const tracesByExposure = useMemberTraces(exposureIds);
  const indicesByExposure = useMemberIndices(exposureIds);

  // Re-key trace/phase maps from exposure id → sample id for the row components.
  const traceBySample = useMemo(() => {
    const m = new Map<number, Trace>();
    for (const id of allSampleIds) {
      const eid = pickerById.get(id);
      if (eid != null && tracesByExposure.has(eid)) m.set(id, tracesByExposure.get(eid)!);
    }
    return m;
  }, [allSampleIds, pickerById, tracesByExposure]);
  const phaseBySample = useMemo(() => {
    const m = new Map<number, string | null>();
    for (const id of allSampleIds) {
      const eid = pickerById.get(id);
      const idx = eid != null ? indicesByExposure.get(eid) : undefined;
      m.set(id, idx ? dominantPhase(idx).dominant : null);
    }
    return m;
  }, [allSampleIds, pickerById, indicesByExposure]);

  const phaseReads = useMemo(
    () => orderedRows.map((r) => {
      const eid = pickerById.get(r.sampleId);
      const idx = eid != null ? indicesByExposure.get(eid) : undefined;
      return idx ? dominantPhase(idx) : { dominant: null, coexist: null };
    }),
    [orderedRows, pickerById, indicesByExposure],
  );

  const [confirmOpen, setConfirmOpen] = useState(false);

  const flagCount = rows.filter((r) => r.flagged).length;
  // Build gate (unchanged contract): every included member must be non-flagged,
  // and there must be at least one — the batch route 400s on an empty array.
  const included = rows.filter((r) => r.include && !r.flagged);
  const canBuild =
    proposal.orderingKey !== undefined &&
    included.length > 0 &&
    rows.every((r) => !r.include || !r.flagged);

  function changeValue(sampleId: number, value: string): void {
    snapshot(`edited smp_${sampleId}`);
    setRows((prev) =>
      prev.map((r) => (r.sampleId === sampleId ? { ...r, value, flagged: value === "" } : r)));
  }
  function toggleFlag(sampleId: number): void {
    snapshot(`toggled smp_${sampleId}`);
    setRows((prev) =>
      prev.map((r) => (r.sampleId === sampleId ? { ...r, flagged: !r.flagged } : r)));
  }
  function addLoose(sampleId: number): void {
    snapshot(`added smp_${sampleId}`);
    const found = loose.find((r) => r.sampleId === sampleId);
    setLoose((prev) => prev.filter((r) => r.sampleId !== sampleId));
    if (found) setRows((prev) => [...prev, { ...found, include: true }]);
  }

  // Defer navigation until the batch write actually succeeds — mirror
  // Compare.tsx: gate the nav on isSuccess via a pending ref (#212 review).
  const pendingBuildRef = useRef(false);

  function handleConfirm(): void {
    if (proposal.orderingKey === undefined) return;
    pendingBuildRef.current = true;
    scopeSeries.mutate({
      key: proposal.orderingKey,
      tags: included.map((r) => ({ sampleId: r.sampleId, value: r.value })),
    });
    setConfirmOpen(false);
  }

  useEffect(() => {
    if (!scopeSeries.isSuccess || !pendingBuildRef.current) return;
    pendingBuildRef.current = false;
    // D1: land on the folio (exists since I3.3).
    navigate("/series");
  }, [scopeSeries.isSuccess, navigate]);

  useEffect(() => {
    if (scopeSeries.error) pendingBuildRef.current = false;
  }, [scopeSeries.error]);

  const keyLabel = humanizeKey(proposal.orderingKey);

  return (
    <div
      data-testid="scoping-page"
      className="flex flex-1 flex-col items-center overflow-auto px-10 py-10"
    >
      <div className="flex w-full max-w-[760px] items-center justify-end pb-3">
        <button
          type="button"
          data-testid="scoping-discard"
          onClick={() => navigate("/series")}
          className="px-[4px] py-[7px] text-xs font-semibold text-ink-faint hover:text-ink"
        >
          Discard
        </button>
      </div>

      <Card
        elevated
        data-testid="scoping-plate"
        className="w-full max-w-[760px] px-8 py-7"
      >
        <Kicker tone="accent" className="mb-1.5">New series</Kicker>
        <h1 data-testid="scoping-title" className="text-display text-ink">
          {proposal.orderingKey ? `Series by ${keyLabel}` : "New series"}
        </h1>

        {scopeSeries.error ? (
          <div
            data-testid="scoping-error-banner"
            role="alert"
            className="mt-4 rounded border border-print-accent bg-paper-sunk px-4 py-2 text-sm text-print-accent"
          >
            Could not write the scoping tags. Nothing was saved. Adjust and try Confirm &amp; build
            again.
          </div>
        ) : null}

        {isError ? (
          <div data-testid="scoping-error" className="px-4 py-8 text-sm text-ink-soft">
            Could not load corpus tags. Try reloading the page.
          </div>
        ) : (
          <Skeleton
            name="scoping"
            className="w-full"
            loading={isLoading}
            stagger={50}
            transition={200}
            fixture={SCOPING_FIXTURE}
          >
            {proposal.orderingKey === undefined ? (
              <div
                data-testid="scoping-empty"
                className="flex flex-col items-center gap-3 px-4 py-12 text-center text-sm text-ink-faint"
              >
                <p>
                  {rows.length === 0 && loose.length === 0
                    ? "No samples in the corpus to scope."
                    : "No shared ordering variable yet. Tag these samples to propose a series."}
                </p>
                <Link
                  to="/samples"
                  data-testid="scoping-empty-cta"
                  className="font-medium text-print-accent hover:underline focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
                >
                  Open the contact sheet
                </Link>
              </div>
            ) : (
              <>
                <ScopingAutogroupCard memberCount={rows.length} keyLabel={keyLabel} flagCount={flagCount} />
                <ScopingOrderField keyLabel={keyLabel} />

                <div className="mb-1 mt-6 flex items-baseline justify-between">
                  <Kicker as="span" tone="faint">The series</Kicker>
                  <span className="flex items-baseline gap-3.5">
                    {history.length > 0 ? (
                      <button
                        type="button"
                        data-testid="scoping-undo"
                        onClick={undo}
                        title={`Step back: ${history[history.length - 1]!.label}`}
                        className="text-sm font-semibold text-print-accent hover:underline"
                      >
                        ↺ Undo last change
                      </button>
                    ) : null}
                    <span className="font-mono text-xs text-ink-faint">
                      {rows.length} samples · low to high
                    </span>
                  </span>
                </div>

                <div>
                  {orderedRows.map((r) => (
                    <ScopingRow
                      key={r.sampleId}
                      row={r}
                      trace={traceBySample.get(r.sampleId)}
                      phase={phaseBySample.get(r.sampleId) ?? null}
                      onChangeValue={(v) => changeValue(r.sampleId, v)}
                      onToggleFlag={() => toggleFlag(r.sampleId)}
                    />
                  ))}
                </div>

                <ScopingLooseMatches
                  rows={loose}
                  traces={traceBySample}
                  phases={phaseBySample}
                  onAdd={addLoose}
                />
                <div className="mt-5">
                  {/* interim inline kicker — replace with <Kicker> when it lands */}
                  <div className="mb-1.5 text-xs font-bold uppercase tracking-wider text-ink-faint">
                    Preview: phase across the series
                  </div>
                  <PhaseStrip
                    size="sm"
                    emptyLabel="Members not yet indexed; phase preview unavailable."
                    segments={phaseReads.map((r) => ({ phase: r.dominant, coexistWith: r.coexist }))}
                  />
                </div>
                <ScopingFoot
                  flagCount={flagCount}
                  memberCount={rows.length}
                  keyLabel={keyLabel}
                  canBuild={canBuild}
                  onBuild={() => setConfirmOpen(true)}
                />
              </>
            )}
          </Skeleton>
        )}
      </Card>

      <ScopingConfirmModal
        open={confirmOpen}
        orderingKey={proposal.orderingKey}
        count={included.length}
        onConfirm={handleConfirm}
        onClose={() => setConfirmOpen(false)}
      />
    </div>
  );
}
