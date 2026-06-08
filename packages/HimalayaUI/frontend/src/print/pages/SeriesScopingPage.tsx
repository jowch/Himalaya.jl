import { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { ScopePlate } from "../components/ScopePlate";
import { ScopeSampleRow } from "../components/ScopeSampleRow";
import { ScopeCandidateRow } from "../components/ScopeCandidateRow";
import { useDragReorder, reorder } from "../components/useDragReorder";
import { EmptyState, Button } from "../ui";
import type { Trace } from "../../api";
import {
  useCorpusSampleTags,
  useCorpusPickerSamples,
  useScopeSeries,
  useMemberTraces,
  useMemberIndices,
} from "../../queries";
import { proposeOrdering, type OrderingRow } from "../../lib/scoping/proposeOrdering";
import { splitProposal, humanizeKey } from "../../lib/scoping/splitProposal";
import { parseSortKey } from "../../lib/scoping/parseSortKey";
import { dominantPhase } from "../../lib/scoping/dominantPhase";
import { buildFootState, canScopeBuild, buildScopePayload, toPreviewSegments } from "./scopingDerive";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

// Token-only skeleton fixture (no inline appearance literals — design-guard clean).
// No scoping.bones.json capture exists yet (deferred, needs the data volume).
const SCOPING_FIXTURE = (
  <div className="space-y-2">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="flex items-center gap-3 rounded border border-hair p-3">
        <div className="h-4 w-4 rounded bg-paper-sunk" />
        <div className="h-4 w-1/3 rounded bg-paper-sunk" />
        <div className="ml-auto h-4 w-1/5 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

type HistoryEntry =
  | { type: "flag"; id: number; prev: boolean; label: string }
  | { type: "add"; id: number; label: string };

/**
 * SeriesScopingPage (greenfield) — the machine-proposes / human-confirms
 * scoping worksheet at /series/new. The confirm-and-build GATE that *writes*
 * the structured (key,value) sample_tags — NOT series creation (that is the
 * builder). Assembled from src/print composites + the series-scoping mockup;
 * carried logic only (proposeOrdering/splitProposal/dominantPhase + the
 * useScopeSeries batch write), no legacy presentation. See the plan's scope
 * table for the six honesty decisions (display-only order field, no confirm
 * modal, display-only drag-reorder, honest grouping copy, flag-not-edit).
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
  const keyLabel = humanizeKey(proposal.orderingKey);

  // Default value-sorted member order (low → high; unparseable last, stable by id).
  const seededOrder = useMemo(
    () =>
      [...split.members]
        .sort((a, b) => {
          const ka = parseSortKey(a.value);
          const kb = parseSortKey(b.value);
          if (ka === null && kb === null) return a.sampleId - b.sampleId;
          if (ka === null) return 1;
          if (kb === null) return -1;
          return ka - kb || a.sampleId - b.sampleId;
        })
        .map((r) => r.sampleId),
    [split.members],
  );

  // Local, user-owned copies seeded once the proposal resolves.
  const [rows, setRows] = useState<OrderingRow[]>([]);
  const [loose, setLoose] = useState<OrderingRow[]>([]);
  const [order, setOrder] = useState<number[]>([]);
  const [history, setHistory] = useState<HistoryEntry[]>([]);
  useEffect(() => {
    setRows(split.members);
    setLoose(split.looseMatches);
    setOrder(seededOrder);
    setHistory([]);
  }, [split.members, split.looseMatches, seededOrder]);

  // Display-only manual reorder (scope decision #5): rewrites `order` + the
  // preview; never touches the written (key,value) payload.
  const { dragItemProps, dropEdge } = useDragReorder((from, to) =>
    setOrder((o) => reorder(o, from, to)),
  );

  const byId = useMemo(() => new Map(rows.map((r) => [r.sampleId, r])), [rows]);
  const sorted = useMemo(
    () => order.map((id) => byId.get(id)).filter((r): r is OrderingRow => r != null),
    [order, byId],
  );

  // Clicking a value confirms / re-opens the machine's read (flag-not-edit,
  // scope decision #2). Recorded so Undo / ⌘Z steps it back.
  const toggleFlag = (id: number): void => {
    const m = rows.find((r) => r.sampleId === id);
    if (!m) return;
    setHistory((h) => [...h, { type: "flag", id, prev: m.flagged, label: `smp_${id}` }]);
    setRows((cur) => cur.map((r) => (r.sampleId === id ? { ...r, flagged: !r.flagged } : r)));
  };

  // Folding a candidate in: include it and APPEND to the displayed order
  // (manual order wins — never re-sort).
  const addCandidate = (id: number): void => {
    const c = loose.find((r) => r.sampleId === id);
    if (!c) return;
    setLoose((cur) => cur.filter((r) => r.sampleId !== id));
    setRows((cur) => [...cur, { ...c, include: true }]);
    setOrder((o) => [...o, id]);
    setHistory((h) => [...h, { type: "add", id, label: `smp_${id}` }]);
  };

  const undo = (): void => {
    setHistory((h) => {
      const e = h[h.length - 1];
      if (!e) return h;
      if (e.type === "flag") {
        setRows((cur) => cur.map((r) => (r.sampleId === e.id ? { ...r, flagged: e.prev } : r)));
      } else {
        setRows((cur) => {
          const m = cur.find((r) => r.sampleId === e.id);
          if (m) setLoose((ls) => [...ls, { ...m, include: false }]);
          return cur.filter((r) => r.sampleId !== e.id);
        });
        setOrder((o) => o.filter((id) => id !== e.id));
      }
      return h.slice(0, -1);
    });
  };

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        undo();
      }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  });

  // Trace + dominant-phase maps, keyed exposure → sample (carried wiring).
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

  // Preview strip follows the displayed order.
  const preview = useMemo(
    () =>
      toPreviewSegments(
        sorted.map((r) => {
          const eid = pickerById.get(r.sampleId);
          const idx = eid != null ? indicesByExposure.get(eid) : undefined;
          return idx ? dominantPhase(idx) : { dominant: null, coexist: null };
        }),
      ),
    [sorted, pickerById, indicesByExposure],
  );

  const flagCount = rows.filter((r) => r.flagged).length;
  const footState = buildFootState(flagCount, rows.length);
  const canBuild = canScopeBuild(rows, proposal.orderingKey);
  const lastLabel = history.length ? history[history.length - 1]!.label : undefined;

  // Defer navigation until the batch write actually succeeds (pending ref).
  const pendingBuildRef = useRef(false);
  const handleBuild = (): void => {
    if (proposal.orderingKey === undefined) return;
    pendingBuildRef.current = true;
    scopeSeries.mutate({ key: proposal.orderingKey, tags: buildScopePayload(rows) });
  };
  useEffect(() => {
    if (!scopeSeries.isSuccess || !pendingBuildRef.current) return;
    pendingBuildRef.current = false;
    navigate("/series");
  }, [scopeSeries.isSuccess, navigate]);
  useEffect(() => {
    if (scopeSeries.error) pendingBuildRef.current = false;
  }, [scopeSeries.error]);

  // ── State 1: corpus load failed (distinct from an empty result). ──────────
  if (isError) {
    return (
      <PageFrame width="scoping" className="px-10 py-10">
        <EmptyState
          title="Couldn't load the corpus"
          body="The sample tags failed to load. Try reloading the page."
        />
      </PageFrame>
    );
  }

  return (
    <PageFrame width="scoping" className="px-10 py-10">
      <div data-testid="scoping-page">
        {/* Discard (mockup top-right, outside the plate). */}
        <div className="flex items-center justify-end pb-3">
          <button
            type="button"
            data-testid="scoping-discard"
            onClick={() => navigate("/series")}
            className="px-1 py-1.5 text-caption font-semibold text-ink-faint hover:text-ink"
          >
            Discard
          </button>
        </div>

        {/* State 4: the batch write failed. */}
        {scopeSeries.error ? (
          <div
            data-testid="scoping-error-banner"
            role="alert"
            className="mb-4 rounded border border-print-accent bg-paper-sunk px-4 py-2 text-meta text-print-accent"
          >
            Could not write the scoping tags. Nothing was saved. Adjust and try Confirm &amp; build again.
          </div>
        ) : null}

        <Skeleton
          name="scoping"
          className="block w-full"
          loading={isLoading}
          stagger={50}
          transition={200}
          fixture={SCOPING_FIXTURE}
          fallback={<div className="p-8 text-meta text-ink-faint">Loading the worksheet…</div>}
        >
          {proposal.orderingKey === undefined ? (
            /* State 3: nothing shares an ordering variable yet. */
            <EmptyState
              title={rows.length === 0 && loose.length === 0 ? "Nothing to scope yet" : "No shared ordering variable"}
              body={
                <div className="flex flex-col items-center gap-4">
                  <span>
                    {rows.length === 0 && loose.length === 0
                      ? "No samples in the corpus to scope."
                      : "These samples share no ordering variable yet. Tag them on the contact sheet to propose a series."}
                  </span>
                  <Button variant="outline" onClick={() => navigate("/samples")}>
                    Open the contact sheet
                  </Button>
                </div>
              }
            />
          ) : (
            <ScopePlate
              seriesName={`Series by ${keyLabel}`}
              grouping={
                <>
                  Himalaya grouped <strong>{rows.length} samples</strong> by their{" "}
                  <strong>{keyLabel}</strong>, read from the sample names.
                  {flagCount > 0 ? (
                    <>
                      {" "}
                      {rows.length - flagCount} parsed cleanly, {flagCount}{" "}
                      {flagCount === 1 ? "needs" : "need"} a look.
                    </>
                  ) : (
                    <> All {rows.length} parsed cleanly.</>
                  )}
                </>
              }
              orderedBy={keyLabel}
              orderNote="Read from the sample names."
              count={`${rows.length} samples · low to high`}
              {...(history.length
                ? { onUndo: undo, ...(lastLabel ? { undoLabel: `Step back: ${lastLabel}` } : {}) }
                : {})}
              rows={sorted.map((r, i) => {
                const dprops = dragItemProps(i);
                const edge = dropEdge(i);
                return (
                  <div
                    key={r.sampleId}
                    {...dprops}
                    className={`relative cursor-grab${dprops["data-dragging"] ? " opacity-50" : ""}`}
                  >
                    {edge && (
                      <span
                        aria-hidden="true"
                        className={`pointer-events-none absolute left-0 right-0 z-10 h-0.5 rounded-full bg-accent ${edge === "top" ? "-top-px" : "-bottom-px"}`}
                      />
                    )}
                    <ScopeSampleRow
                      name={r.sampleName}
                      sampleId={`smp_${r.sampleId}`}
                      trace={traceBySample.get(r.sampleId) ?? EMPTY_TRACE}
                      phase={phaseBySample.get(r.sampleId) ?? null}
                      value={r.value}
                      {...(r.flagged ? { flagged: true } : {})}
                      onToggleFlag={() => toggleFlag(r.sampleId)}
                    />
                  </div>
                );
              })}
              candidates={
                loose.length ? (
                  loose.map((c) => (
                    <ScopeCandidateRow
                      key={c.sampleId}
                      name={c.sampleName}
                      why={
                        <>
                          lacks the <strong className="text-accent font-semibold">{keyLabel}</strong> — add
                          it if it belongs.
                        </>
                      }
                      trace={traceBySample.get(c.sampleId) ?? EMPTY_TRACE}
                      phase={phaseBySample.get(c.sampleId) ?? null}
                      onAdd={() => addCandidate(c.sampleId)}
                    />
                  ))
                ) : (
                  <div className="text-meta text-ink-faint italic">
                    Nothing else in the corpus matches this grouping.
                  </div>
                )
              }
              preview={preview}
              footState={footState}
              footNote={
                <>
                  Confirming records the {keyLabel} on every sample — the next series that needs it already
                  knows.
                </>
              }
              {...(!canBuild ? { buildDisabled: true } : {})}
              onBuild={handleBuild}
            />
          )}
        </Skeleton>
      </div>
    </PageFrame>
  );
}
