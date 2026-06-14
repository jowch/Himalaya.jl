import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { ReactNode } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { ScopePlate } from "../components/ScopePlate";
import { ScopeSampleRow } from "../components/ScopeSampleRow";
import { useDragReorder, reorder } from "../components/useDragReorder";
import { Sparkline } from "../plot/Sparkline";
import { EmptyState, Button, Card, Dot, Field, Input, Kicker, PhaseChip } from "../ui";
import { ColdAssignPanel } from "../components/ColdAssignPanel";
import type { SampleTagPair, Trace } from "../../api";
import {
  useCorpusSampleTags,
  useCorpusPickerSamples,
  useScopeSeries,
  useCreateSeries,
  useMemberTraces,
  useMemberIndices,
} from "../../queries";
import { proposeOrdering, type OrderingRow } from "../../lib/scoping/proposeOrdering";
import { splitProposal, humanizeKey } from "../../lib/scoping/splitProposal";
import { parseSortKey } from "../../lib/scoping/parseSortKey";
import { dominantPhase } from "../../lib/scoping/dominantPhase";
import {
  buildFootState,
  canScopeBuild,
  buildScopePayload,
  toPreviewSegments,
  isKept,
  filterPickerBySeed,
  buildColdAssignRows,
  canColdBuild,
  buildColdScopePayload,
  type ColdAssignRow,
} from "./scopingDerive";
import { readNewSeriesSeed } from "../../lib/series/newSeriesNav";
import { suppressGlobalKeys } from "../../lib/keys";
import { useReorderShortcuts } from "../shell/useReorderShortcuts";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

// The sentinel "Define your own..." dropdown entry (locked decision A):
// selecting it inline-converts the warm worksheet to the ColdAssignPanel flow
// (seeded visits carry the WHOLE selection, members first then loose; direct
// visits carry the current members only), committing through the SAME two-op
// scope->create chain. value === label so Menu's strict activeValue equality
// highlights it while custom mode is active.
const DEFINE_YOUR_OWN = "Define your own…";

// How many loose-candidate exemplars the worksheet renders. The worksheet is a
// review surface; the contact sheet is the corpus browser. On a whole-corpus
// visit the loose list can run to 100+ samples — rendering them all buries the
// members and pushes the build action far below the fold, and discovery beyond
// a few exemplars belongs on the contact sheet. The section note below the
// rows states the hidden remainder and links there.
const SCOPE_CANDIDATE_PREVIEW_COUNT = 3;

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

// Undo history (⌘Z / the plate's "Undo last change"). Both worksheet edit
// gestures push a typed entry so each is recoverable — skip and reorder reach
// parity (no asymmetric freedom): a `flag` entry restores one member's skip
// state; a `reorder` entry restores the WHOLE prior display order.
type HistoryEntry =
  | { type: "flag"; id: number; prev: boolean; label: string }
  | { type: "reorder"; prevOrder: number[]; label: string }
  // SC-VALUECORRECT: an inline value correction is recoverable too — `prev` is
  // the pre-edit read so skip/reorder/edit reach undo parity.
  | { type: "value"; id: number; prev: string; label: string }
  // SCOPE-LOOSEADD: pulling a loose candidate into the members is recoverable —
  // `prevLoose` is the original value-less candidate so undo restores it whole
  // to the loose list (carried in the entry so undo needs no live `rows` read,
  // keeping every history updater pure — no StrictMode double-apply).
  | { type: "add"; id: number; prevLoose: OrderingRow; label: string };

/**
 * ColdAssignSection: the shared assign body (ColdAssignPanel + the gated
 * confirm foot) used by BOTH the cold path (no proposable key + a seed) and
 * the warm worksheet's "Define your own..." custom mode. Presentational; the
 * page owns the state and the build chain (`onBuild` is handleBuild in both
 * uses, so each commits through the same two-op scope->create chain).
 */
function ColdAssignSection({
  rows,
  variableKey,
  onKeyChange,
  onValueChange,
  canBuildNow,
  onBuild,
  buildBusy,
  intro,
  discardSlot,
}: {
  rows: ColdAssignRow[];
  variableKey: string;
  onKeyChange: (key: string) => void;
  onValueChange: (sampleId: number, value: string) => void;
  canBuildNow: boolean;
  onBuild: () => void;
  /** Discard affordance, placed left of "Confirm & build" (SC-POLISH2). */
  discardSlot?: ReactNode;
  /** The scope→create chain is in flight: the foot button flips to the
   *  progressive register ("Building…") with `aria-busy`, disabled (no
   *  double-submit). Same contract as ScopePlate's buildBusy. */
  buildBusy: boolean;
  /** Pass-through to ColdAssignPanel: undefined keeps the default cold-corpus
   *  intro; an explicit null suppresses it (custom mode, where that copy
   *  would be false). */
  intro?: ReactNode;
}): JSX.Element {
  return (
    <>
      <ColdAssignPanel
        rows={rows}
        variableKey={variableKey}
        onKeyChange={onKeyChange}
        onValueChange={onValueChange}
        {...(intro !== undefined ? { intro } : {})}
      />
      <div className="mt-6 pt-4 border-t border-hair flex items-center justify-between gap-5">
        <div className="flex flex-col gap-1">
          <div
            className={`flex items-center gap-2 text-meta font-semibold ${
              canBuildNow ? "text-ink" : "text-accent"
            }`}
          >
            <Dot tone={canBuildNow ? "success" : "accent"} aria-hidden />
            {canBuildNow
              ? `${rows.length} value${rows.length === 1 ? "" : "s"} ready to commit`
              : "Name the variable and assign every value to build"}
          </div>
          <div className="text-caption text-ink-soft max-w-[42ch]">
            Confirming records the variable on every sample. The next series that needs it
            already knows.
          </div>
        </div>
        <div className="flex items-center gap-2 flex-shrink-0">
          {discardSlot}
          <Button
            variant="solid"
            {...(!canBuildNow || buildBusy ? { disabled: true } : {})}
            {...(buildBusy ? { "aria-busy": true } : {})}
            onClick={onBuild}
          >
            {buildBusy ? "Building…" : <>Confirm &amp; build →</>}
          </Button>
        </div>
      </div>
    </>
  );
}

/**
 * SeriesScopingPage (greenfield) — the machine-proposes / human-confirms
 * scoping worksheet at /series/new. The confirm-and-build GATE that *writes*
 * the structured (key,value) sample_tags AND creates the series (M-A): a
 * TWO-OP queue chain (mirrors the builder's save→commit) — Op A writes the
 * ordering tags (useScopeSeries), Op B creates the series (useCreateSeries,
 * POST /api/series) — sequenced by a ref state machine; on success the page
 * navigates to the new `/series/:id` builder. Assembled from src/print
 * composites + the series-scoping mockup; carried logic only
 * (proposeOrdering/splitProposal/dominantPhase + the scope→create chain), no
 * legacy presentation.
 *
 * Honest commit-gate model (Option A): in the carried data a member ALWAYS has a
 * value (loose matches — value==="" — split out as informational candidates), so
 * members are never machine-flagged. Scoping's only durable effect is writing
 * Himalaya's parsed (key,value) read onto each member; every read is committed
 * unless the user SKIPS it (the value control toggles skip). Candidates carry no
 * value, so they cannot be committed — they're discovery only, with no add path.
 */
export function SeriesScopingPage(): JSX.Element {
  const navigate = useNavigate();
  const location = useLocation();
  const tagsQ = useCorpusSampleTags();
  const pickerQ = useCorpusPickerSamples();
  const scopeSeries = useScopeSeries();
  const createSeries = useCreateSeries();

  const isLoading = tagsQ.isLoading || pickerQ.isLoading;
  const isError = tagsQ.isError || pickerQ.isError;

  // Seeded selection (from the contact-sheet picker → /series/new): when present
  // we scope BOTH the proposal rows AND the ordering-variable vocabulary to just
  // those samples; when null (a direct visit) the full corpus is used — the
  // existing whole-corpus behaviour.
  const seed = useMemo(() => readNewSeriesSeed(location), [location]);
  const seededPicker = useMemo(
    () => filterPickerBySeed(pickerQ.data ?? [], seed),
    [pickerQ.data, seed],
  );

  // SC-SEEDDEAD: the ordering-variable vocabulary scopes to the seed. A seeded
  // visit derives it from the seeded samples' OWN tags, deduped to first match
  // per (sample, key) — retries can duplicate sample_tags rows and readers take
  // the first match — so frequency = how many seeded samples carry the key; a
  // direct visit keeps the corpus-wide distinct pair list. Without this, a seed
  // of samples lacking the corpus's dominant key dead-ends in a 0-member
  // worksheet whose foot gate can never be satisfied.
  const scopedTags = useMemo<SampleTagPair[]>(
    () =>
      seed === null
        ? tagsQ.data ?? []
        : seededPicker.flatMap((s) => {
            const seen = new Set<string>();
            return s.sample.tags
              .filter((t) => {
                if (seen.has(t.key)) return false;
                seen.add(t.key);
                return true;
              })
              .map((t) => ({ key: t.key, value: t.value }));
          }),
    [seed, tagsQ.data, seededPicker],
  );

  // The real ordering-variable vocabulary: the distinct tag keys of the scoped
  // universe (the seed's own tags, or the whole corpus on a direct visit), in
  // proposeOrdering's deterministic order (frequency desc, lexicographic tie),
  // so the dropdown lists exactly what the scoped samples expose — no declared
  // schema, and no corpus key the seed can't satisfy.
  const orderKeys = useMemo(() => {
    const freq = new Map<string, number>();
    for (const t of scopedTags) freq.set(t.key, (freq.get(t.key) ?? 0) + 1);
    return [...freq.entries()]
      .sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1))
      .map(([k]) => k);
  }, [scopedTags]);

  // Human override of the machine's frequency winner. Null → the machine picks.
  const [overrideKey, setOverrideKey] = useState<string | null>(null);
  // A vocabulary change can't strand a stale override: drop it once it names a
  // key the scoped universe (seed's tags, or the corpus on a direct visit) no
  // longer exposes.
  useEffect(() => {
    if (overrideKey !== null && !orderKeys.includes(overrideKey)) setOverrideKey(null);
  }, [overrideKey, orderKeys]);

  const proposal = useMemo(
    () => proposeOrdering(scopedTags, seededPicker, overrideKey ?? undefined),
    [scopedTags, seededPicker, overrideKey],
  );
  const split = useMemo(() => splitProposal(proposal), [proposal]);
  const keyLabel = humanizeKey(proposal.orderingKey);

  // Ordering-variable dropdown wiring. The Menu highlights its active option by
  // strict value-equality with the displayed value (`orderedBy` === keyLabel),
  // so each option's `value` IS its humanized label; map the selected label back
  // to the raw key for the override. Deduped by humanized label so two raw keys
  // can't collide on one option (last raw key wins the label). A real key whose
  // humanized label collides with the sentinel is filtered out so it can never
  // shadow the "Define your own..." routing.
  const labelToKey = useMemo(
    () =>
      new Map(
        orderKeys
          .map((k) => [humanizeKey(k), k] as const)
          .filter(([label]) => label !== DEFINE_YOUR_OWN),
      ),
    [orderKeys],
  );
  // The dropdown always renders (even with a single key): "Define your own..."
  // is a real alternative on every corpus, listed last.
  const orderOptions = useMemo(
    () => [
      ...[...labelToKey.keys()].map((label) => ({ value: label, label })),
      { value: DEFINE_YOUR_OWN, label: DEFINE_YOUR_OWN },
    ],
    [labelToKey],
  );

  // Custom ("define your own") mode is page state, not URL state. Only
  // reachable from the warm worksheet; the cold path renders independently.
  const [customMode, setCustomMode] = useState(false);

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

  // Local, user-owned copies seeded once the proposal resolves. The reseed on a
  // `split` change is intentional: TanStack structural sharing keeps the query
  // `data` referentially stable unless the corpus content actually changed, in
  // which case reseeding the worksheet to the new proposal is correct.
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

  // ── Cold path (no variable proposable FROM THE SEED, user arrived with one) ──
  const [coldKey, setColdKey] = useState("");
  const [coldRows, setColdRows] = useState<ColdAssignRow[]>([]);
  // SC-POLISH2: a two-step guard so Discard never destroys in-progress work
  // silently. The first click on a dirty worksheet arms this inline confirm;
  // a clean worksheet discards immediately (nothing to lose).
  const [confirmingDiscard, setConfirmingDiscard] = useState(false);
  // Seed cold rows once when the page loads into the cold path.
  useEffect(() => {
    if (!isLoading && proposal.orderingKey === undefined && seed !== null && rows.length === 0) {
      const names = new Map<number, string>();
      for (const s of pickerQ.data ?? [])
        names.set(s.sample.id, s.sample.display_name ?? s.sample.name ?? "");
      setColdRows(
        buildColdAssignRows(
          (seed ?? []).map((id) => ({ sampleId: id, sampleName: names.get(id) ?? `smp_${id}` })),
        ),
      );
      setColdKey("");
    }
  }, [isLoading, proposal.orderingKey, seed, rows.length, pickerQ.data]);

  // Display-only manual reorder (scope decision #5): rewrites `order` + the
  // preview; never touches the written (key,value) payload. `applyReorder` is
  // the SINGLE order mutation — the pointer (drag) and keyboard (grip-button
  // arrows) paths both converge on it, so recording the undo entry HERE covers
  // both. A no-op move (boundary press, identical indices) must not push a dead
  // history entry — `reorder` returns an equal-content array, so compare and
  // bail. The change-detection + both setState calls live OUTSIDE any updater:
  // pushing `setHistory` from inside the `setOrder` updater would be an impure
  // updater that StrictMode double-invokes (a single reorder would mint two
  // history entries). `order` from the render closure is the pre-mutation array,
  // captured for the undo restore.
  const applyReorder = useCallback(
    (from: number, to: number): void => {
      const next = reorder(order, from, to);
      const changed =
        next.length !== order.length || next.some((id, i) => id !== order[i]);
      if (!changed) return;
      setHistory((h) => [...h, { type: "reorder", prevOrder: order, label: "reorder" }]);
      setOrder(next);
    },
    [order],
  );
  const { dragItemProps, dropEdge } = useDragReorder(applyReorder);

  const byId = useMemo(() => new Map(rows.map((r) => [r.sampleId, r])), [rows]);
  const sorted = useMemo(
    () => order.map((id) => byId.get(id)).filter((r): r is OrderingRow => r != null),
    [order, byId],
  );

  // ── Keyboard reorder (SC-KBD, WCAG 2.1.1) ─────────────────────────────────
  // Announcements land in a persistent polite live region (rendered
  // unconditionally below, never mounted-on-demand, so SRs reliably pick up
  // text changes). The toggled trailing space mirrors LiveRegion.tsx: two
  // identical consecutive announcements (e.g. ArrowUp on the first row twice)
  // must still mutate textContent or the SR silently skips the repeat.
  const [reorderMsg, setReorderMsg] = useState("");
  const reorderFlip = useRef(false);
  const announceReorder = useCallback((msg: string): void => {
    reorderFlip.current = !reorderFlip.current;
    setReorderMsg(reorderFlip.current ? msg : `${msg} `);
  }, []);

  // Arrow-key move from row i's grip button. A boundary press is an order
  // no-op but still announces — a silent keypress is dead air for an SR user.
  const moveRow = (i: number, delta: -1 | 1): void => {
    const r = sorted[i];
    if (!r) return;
    if (i === 0 && delta === -1) {
      announceReorder(`${r.sampleName} is already first.`);
      return;
    }
    if (i === sorted.length - 1 && delta === 1) {
      announceReorder(`${r.sampleName} is already last.`);
      return;
    }
    applyReorder(i, i + delta);
    announceReorder(`Moved ${r.sampleName} to position ${i + delta + 1} of ${sorted.length}.`);
  };

  // Unified Alt+↑/↓ reorder power-gesture (shared registry). Moves the focused
  // worksheet row from any control inside it, mirroring the grip's own arrows;
  // moveRow owns the clamp + SR announcement. The keyed row wrapper carries the
  // data-reorder-index the hook reads.
  useReorderShortcuts({
    rowSelector: "[data-reorder-row]",
    move: (i, delta) => moveRow(i, delta),
  });

  // Routing for the ordering dropdown. The sentinel enters custom mode, seeding
  // the assign rows conditionally on the seed (re-entry deliberately re-seeds:
  // previous custom edits are discarded). A SEEDED visit carries the WHOLE
  // selection — members in displayed order first, then the loose ones (the user
  // deliberately picked them all; custom mode must not silently drop any). A
  // direct visit stays members-only: a whole-corpus DYO must not fan in 100+
  // loose candidates. Any existing label exits custom mode and applies the
  // override (existing behaviour).
  const onOrderSelect = useCallback(
    (label: string): void => {
      if (label === DEFINE_YOUR_OWN) {
        const source = seed !== null ? [...sorted, ...loose] : sorted;
        setColdRows(
          buildColdAssignRows(
            source.map((r) => ({ sampleId: r.sampleId, sampleName: r.sampleName })),
          ),
        );
        setColdKey("");
        setCustomMode(true);
        return;
      }
      setCustomMode(false);
      setOverrideKey(labelToKey.get(label) ?? null);
    },
    [labelToKey, sorted, loose, seed],
  );

  const onColdValueChange = useCallback(
    (id: number, v: string): void =>
      setColdRows((cur) => cur.map((r) => (r.sampleId === id ? { ...r, value: v } : r))),
    [],
  );

  // Clicking a member's value toggles "skip this sample from the write" — the
  // honest commit-gate (Option A). A skipped read is excluded from the batch,
  // never blocks the build. Recorded so Undo / ⌘Z steps it back.
  const toggleFlag = (id: number): void => {
    const m = rows.find((r) => r.sampleId === id);
    if (!m) return;
    setHistory((h) => [...h, { type: "flag", id, prev: m.flagged, label: `smp_${id}` }]);
    setRows((cur) => cur.map((r) => (r.sampleId === id ? { ...r, flagged: !r.flagged } : r)));
  };

  // SC-VALUECORRECT: correct a misread auto-read in place. The ScopeSampleRow
  // control already guards changed + non-empty, but re-check here so a stray
  // no-op caller can never push an empty history entry. Recorded for Undo/⌘Z.
  const editValue = (id: number, value: string): void => {
    const m = rows.find((r) => r.sampleId === id);
    if (!m || value === m.value) return;
    setHistory((h) => [...h, { type: "value", id, prev: m.value, label: `smp_${id}` }]);
    setRows((cur) => cur.map((r) => (r.sampleId === id ? { ...r, value } : r)));
  };

  // SCOPE-LOOSEADD: per-loose-candidate draft of the value the user is naming.
  // A loose candidate lacks a value for the ordering key, so naming one is the
  // act that makes it committable — the inline value field is the labelling
  // gesture the "Himalaya also found" section used to lack.
  const [looseDraft, setLooseDraft] = useState<Record<number, string>>({});
  // Pull a loose candidate into the members WITH the named value. The caller
  // (the Add button / Enter) only fires with a non-empty value; re-check here so
  // a stray no-op can never seat a `value:""` member. The candidate becomes a
  // real, included, unflagged member; its id joins the display order at the end.
  const addLoose = (id: number, value: string): void => {
    const v = value.trim();
    if (v === "") return;
    const c = loose.find((r) => r.sampleId === id);
    if (!c) return;
    setHistory((h) => [...h, { type: "add", id, prevLoose: c, label: `smp_${id}` }]);
    setLoose((cur) => cur.filter((r) => r.sampleId !== id));
    // include:true is load-bearing — loose matches carry include:false, and the
    // commit gate (isKept) drops a non-included row, so a spread that kept the
    // candidate's false would seat a member that silently never commits.
    setRows((cur) => [...cur, { ...c, value: v, flagged: false, include: true }]);
    setOrder((cur) => [...cur, id]);
    setLooseDraft((d) => {
      const next = { ...d };
      delete next[id];
      return next;
    });
  };

  const undo = useCallback((): void => {
    setHistory((h) => {
      const e = h[h.length - 1];
      if (!e) return h;
      if (e.type === "flag") {
        setRows((cur) => cur.map((r) => (r.sampleId === e.id ? { ...r, flagged: e.prev } : r)));
      } else if (e.type === "value") {
        // value: restore the pre-edit read.
        setRows((cur) => cur.map((r) => (r.sampleId === e.id ? { ...r, value: e.prev } : r)));
      } else if (e.type === "add") {
        // add: drop the seated member + its order slot, restore the candidate to
        // the FRONT of the loose pool — undo should pop it back where you grabbed
        // it (visible above the preview cap), not bury it past the hidden tail.
        setRows((cur) => cur.filter((r) => r.sampleId !== e.id));
        setOrder((cur) => cur.filter((oid) => oid !== e.id));
        setLoose((cur) => [e.prevLoose, ...cur]);
      } else {
        // reorder: restore the whole prior display order.
        setOrder(e.prevOrder);
      }
      return h.slice(0, -1);
    });
  }, []);

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        // Guard BEFORE preventDefault — ⌘Z inside an input must stay the
        // browser's native text undo, not the page's skip-undo.
        if (suppressGlobalKeys(e)) return;
        e.preventDefault();
        undo();
      }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [undo]);

  // Trace + dominant-phase maps, keyed exposure → sample (carried wiring).
  const pickerById = useMemo(() => {
    const m = new Map<number, number | null>();
    for (const s of pickerQ.data ?? []) m.set(s.sample.id, s.indexing_exposure_id);
    return m;
  }, [pickerQ.data]);
  // Only the candidate exemplars that actually RENDER (the preview cap) join
  // the trace/index fan-out — the hidden remainder must not trigger fetches.
  // EmptyState gating and the remainder count keep using the FULL loose list.
  const visibleLoose = useMemo(() => loose.slice(0, SCOPE_CANDIDATE_PREVIEW_COUNT), [loose]);
  const hiddenLooseCount = loose.length - visibleLoose.length;
  const allSampleIds = useMemo(
    () => [...rows, ...visibleLoose].map((r) => r.sampleId),
    [rows, visibleLoose],
  );
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

  // Preview strip: the strip's job is to preview the series figure that will
  // be BUILT, so it shows exactly the commit membership (isKept — the same
  // predicate as the write payload), in displayed order. Skipped members must
  // not paint a cell the foot says is excluded.
  const preview = useMemo(
    () =>
      toPreviewSegments(
        sorted.filter(isKept).map((r) => {
          const eid = pickerById.get(r.sampleId);
          const idx = eid != null ? indicesByExposure.get(eid) : undefined;
          return idx ? dominantPhase(idx) : { dominant: null, coexist: null };
        }),
      ),
    [sorted, pickerById, indicesByExposure],
  );

  const skippedCount = rows.filter((r) => r.flagged).length;
  const keptCount = rows.filter(isKept).length;

  // SC-COUNTHONEST: the count caption must not claim "low to high" once the user
  // has manually reordered (drag or keyboard). The canonical default IS the
  // value-sorted `seededOrder`; while the live `order` still equals it the
  // ordering is genuinely low→high, otherwise it is a custom order. Compared by
  // content (same ids, same positions) so a reseed that produces an equal array
  // still reads as the default. (Lengths can differ transiently mid-reseed.)
  const isCanonicalOrder =
    order.length === seededOrder.length && order.every((id, i) => id === seededOrder[i]);
  // SC-PROV: when some reads are skipped, the total-count caption ("N samples")
  // and the foot ("M values ready · K skipped") differ by the skip count. Annotate
  // the skip inline so the two reconcile (N total = M committed + K skipped) instead
  // of silently diverging — the reader shouldn't have to do the subtraction.
  const orderCaption = `${rows.length} samples${skippedCount > 0 ? ` (${skippedCount} skipped)` : ""} · ${isCanonicalOrder ? "low to high" : "custom order"}`;
  const footState = buildFootState(keptCount, skippedCount);
  const canBuild = canScopeBuild(rows, proposal.orderingKey);
  const lastLabel = history.length ? history[history.length - 1]!.label : undefined;

  // SC-POLISH2: Discard moved from the page top into each plate foot beside
  // "Confirm & build", and guarded against silently destroying in-progress
  // work. "Work" = any warm gesture recorded in history (reorder/skip) OR a
  // typed cold/custom value. A clean worksheet skips the confirm.
  const hasUnsavedWork =
    history.length > 0 || coldRows.some((r) => (r.value ?? "").trim() !== "");
  const onDiscardClick = (): void => {
    if (hasUnsavedWork) setConfirmingDiscard(true);
    else navigate("/series");
  };
  const discardControl = confirmingDiscard ? (
    <span
      data-testid="scoping-discard-confirm"
      className="flex items-center gap-2"
      role="group"
      aria-label="Confirm discard"
    >
      <span className="text-caption text-ink-soft">Discard your work?</span>
      <Button variant="ghost" onClick={() => setConfirmingDiscard(false)}>
        Keep editing
      </Button>
      <Button
        variant="ghost"
        data-testid="scoping-discard-yes"
        onClick={() => navigate("/series")}
      >
        Discard
      </Button>
    </span>
  ) : (
    <Button variant="ghost" data-testid="scoping-discard" onClick={onDiscardClick}>
      Discard
    </Button>
  );

  // Cold-corpus build gate (controls-don't-lie: enabled only once the user has
  // named the variable AND filled in every sample's value).
  const canColdBuildNow = canColdBuild(coldKey, coldRows);
  const isColdPath = proposal.orderingKey === undefined && seed !== null;

  // ── Confirm & build: a TWO-OP chain (mirrors the builder's save→commit). ──
  // Confirm & build can't be ONE queue op: a compound op that writes the tags
  // then creates the series would have its tag-write `add_tag` SSE frame
  // resolve the op's deferred BEFORE the in-flight create returned, so
  // `mutation.data` would be the tag confirmation, not the Series (the queue
  // resolves the deferred with whichever lands first — see useQueueMutation).
  // So it's two SEPARATE single-write ops sequenced by a ref state machine:
  //   Op A  scopeSeries.mutate({key, tags})   — the batch sample_tags write
  //   Op B  createSeries.mutate({title, …})   — POST /api/series
  // Only the single-write create reliably surfaces the full Series in
  // `createSeries.data` (its only own-op frame is `series_created`).
  //
  // LOAD-BEARING (mirrors SeriesBuilderPage): the page never calls
  // scopeSeries.reset()/createSeries.reset(). A lingering isSuccess===true
  // between Confirm runs is INERT because every chain effect is `stage`-ref
  // gated — they fire only while stage.current is the matching phase, and
  // handleBuild re-arms by flipping stage to "tagging" itself. Do NOT add a
  // reset() (it would race the SSE-vs-HTTP deferred resolution) or drop the
  // gate (a stale isSuccess would re-fire the next op on an unrelated render).
  const stage = useRef<"idle" | "tagging" | "creating">("idle");
  // The create body, stashed when Op A fires so Op B can read it on success.
  const createBodyRef = useRef<{
    title: string;
    samples: { sample_id: number; position: number }[];
    ordering_variable: string;
  } | null>(null);

  const handleBuild = (): void => {
    if (stage.current !== "idle") return;
    if (isColdPath || customMode) {
      if (!canColdBuildNow) return;
      const key = coldKey.trim();
      // Cold/custom members: every assigned sample, in the worksheet order.
      createBodyRef.current = {
        title: `Series by ${key}`,
        samples: coldRows.map((r, i) => ({ sample_id: r.sampleId, position: i })),
        ordering_variable: key,
      };
      stage.current = "tagging";
      scopeSeries.mutate({ key, tags: buildColdScopePayload(coldRows) });
    } else {
      if (proposal.orderingKey === undefined) return;
      // Warm members: the SCOPED, kept set (skips excluded — isKept, the same
      // predicate as buildScopePayload) in the displayed low→high order, so
      // the series recipe matches the tags actually written.
      const keptInOrder = sorted.filter(isKept);
      createBodyRef.current = {
        title: `Series by ${keyLabel}`,
        samples: keptInOrder.map((r, i) => ({ sample_id: r.sampleId, position: i })),
        ordering_variable: proposal.orderingKey,
      };
      stage.current = "tagging";
      // SC-TAGORDER: write tags from the SAME displayed/reordered set as the
      // position-create above (`sorted`, which buildScopePayload re-filters by
      // isKept → keptInOrder), so tag-write order and series position agree.
      scopeSeries.mutate({ key: proposal.orderingKey, tags: buildScopePayload(sorted) });
    }
  };

  // Op A landed (tags written) → fire Op B (create the series) from the
  // stashed body.
  useEffect(() => {
    if (stage.current !== "tagging" || !scopeSeries.isSuccess) return;
    const body = createBodyRef.current;
    if (!body) {
      stage.current = "idle";
      return;
    }
    stage.current = "creating";
    createSeries.mutate(body);
  }, [scopeSeries.isSuccess, createSeries]);

  // Op B landed (series created) → navigate to the new builder. Read the new
  // series `id` DIRECTLY, never gated on isFullSeries: under the SSE-wins race
  // `createSeries.data` is the PARTIAL shape from saveSeriesMutator's
  // synthesizeFromSse (`{...base, ...payload, id: entity_id}` — no
  // members/state, so isFullSeries is false) — but it DOES carry the new
  // series id. On the HTTP-wins path it's the full Series; the id is present in
  // BOTH shapes. (Live, Op A's tag write warms the /api/events stream so the
  // `series_created` frame usually wins the race; mocked tests drain SSE so
  // HTTP wins — only the id-direct read works for both.)
  useEffect(() => {
    if (stage.current !== "creating" || !createSeries.isSuccess) return;
    stage.current = "idle";
    const created = createSeries.data as { id?: unknown } | undefined;
    const newId = typeof created?.id === "number" ? created.id : undefined;
    navigate(newId !== undefined ? `/series/${newId}` : "/series");
  }, [createSeries.isSuccess, createSeries.data, navigate]);

  // Either op errored → reset so the user can retry.
  useEffect(() => {
    if (scopeSeries.error || createSeries.error) stage.current = "idle";
  }, [scopeSeries.error, createSeries.error]);

  // Busy = the chain is in flight, derived from the SAME stage/mutation
  // sources the chain effects read (never a parallel flag). Each stage clause
  // carries its own truthful-exit term (the error that ends that stage) so the
  // busy register reverts IN THE SAME RENDER the failure banner surfaces — the
  // error exit flips the stage REF in an effect, which triggers no re-render
  // of its own, so a plain `stage.current !== "idle"` would leave a lying
  // "Building…" stuck next to the banner. Scoping the error term per stage
  // also keeps a STALE createSeries.error from a prior attempt (see
  // chainErrorCopy) from suppressing busy during a retry's tagging phase.
  const buildBusy =
    scopeSeries.isPending ||
    createSeries.isPending ||
    (stage.current === "tagging" && scopeSeries.error == null) ||
    (stage.current === "creating" && createSeries.error == null);

  // Chain-failure banner copy, differentiated by WHICH op failed — the two
  // modes leave different data committed and the banner must not lie:
  //   Op A failed → no tags from THIS attempt, no series. (Deliberately not
  //   "Nothing was saved": a PRIOR attempt's tags may already persist.)
  //   Op B failed → Op A already committed, so the ordering tags ARE durably
  //   written; only the series create is missing. A retry re-inserts the same
  //   (key, value) rows (sample_tags has no unique constraint); readers take
  //   the first match, so identical values are a no-op in effect.
  // scopeSeries.error is checked FIRST: handleBuild re-fires scopeSeries.mutate
  // (clearing its error), so when both are set the scope error is the more
  // recent attempt and a stale createSeries.error may linger from a prior one.
  const chainErrorCopy = scopeSeries.error
    ? "Could not save the ordering tags. The series was not created. Adjust and try Confirm & build again."
    : createSeries.error
      ? "The ordering tags were saved, but the series was not created. Try Confirm & build again."
      : null;

  // SC-RETRYBANNER: once a retry is in flight, suppress the prior failure banner
  // (and its "try again" copy) — it would otherwise keep shouting next to the
  // "Building…" control. A stale createSeries.error in particular lingers through
  // the retry's tagging phase (handleBuild only clears scopeSeries.error), so
  // gate on the SAME single in-flight truth the busy button reads. On resolution
  // the banner returns only if a FRESH error surfaces.
  const showErrorBanner = chainErrorCopy !== null && !buildBusy;

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
        {/* Discard now lives in each plate foot beside "Confirm & build"
            (SC-POLISH2), with a confirm-before-destroy guard — see
            `discardControl`. */}

        {/* SC-KBD: keyboard-reorder live region. Persistent (never conditional)
            so screen readers track its text mutations; sr-only is the project's
            visually-hidden mechanism (see ui/LiveRegion.tsx). */}
        <div
          data-testid="reorder-announcement"
          role="status"
          aria-live="polite"
          aria-atomic="true"
          className="sr-only"
        >
          {reorderMsg}
        </div>

        {/* State 4: a chain op failed. Two modes (see chainErrorCopy): an Op-A
            (tag write) failure means no series and no tags from this attempt; an
            Op-B (create) failure means the tags ARE committed — the copy admits
            it. Op A wins precedence because a stale create error can linger. */}
        {showErrorBanner ? (
          <div
            data-testid="scoping-error-banner"
            role="alert"
            className="mb-4 rounded border border-print-accent bg-paper-sunk px-4 py-2 text-meta text-print-accent"
          >
            {chainErrorCopy}
          </div>
        ) : null}

        <Skeleton
          name="scoping"
          className="block w-full"
          loading={isLoading}
          stagger={50}
          transition={200}
          fixture={SCOPING_FIXTURE}
          fallback={<div className="p-8 text-meta text-ink-soft">Loading the worksheet…</div>}
        >
          {proposal.orderingKey === undefined ? (
            seed !== null ? (
              /* Cold path: the user arrived from the contact sheet with a
                 deliberate selection, but no tag key can be proposed from the
                 seeded samples' own tags (SC-SEEDDEAD: the corpus may still
                 carry keys — they don't bind here). Rather than dead-end, let
                 them name the ordering variable and assign each sample's
                 value — then commit through the SAME scope→create chain the
                 warm path uses. SC-COLD/SC-COLDHEAD: the worksheet carries the
                 SAME plate identity as the warm/custom plates ("New series"
                 kicker + serif h1 following the typed key), so ColdAssignPanel's
                 h2 section label nests under a real h1 with no level skip.
                 No order Field here: there is no proposable key to offer, and
                 a dropdown with nothing to select would lie. */
              <Card
                elevated
                data-testid="cold-scope-plate"
                className="w-full max-w-[760px] px-8 pt-7 pb-6"
              >
                <Kicker tone="accent">New series</Kicker>
                <h1 className="text-display text-ink mt-1.5">
                  Series by {coldKey.trim() || "…"}
                </h1>
                <div className="mt-6">
                  <ColdAssignSection
                    rows={coldRows}
                    variableKey={coldKey}
                    onKeyChange={setColdKey}
                    onValueChange={onColdValueChange}
                    canBuildNow={canColdBuildNow}
                    onBuild={handleBuild}
                    buildBusy={buildBusy}
                    discardSlot={discardControl}
                  />
                </div>
              </Card>
            ) : (
              /* State 3: nothing shares an ordering variable yet (no seed). */
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
            )
          ) : customMode ? (
            /* Custom ("Define your own...") mode: the warm worksheet swaps for
               the cold-assign flow, seeded with the whole selection on a seeded
               visit (members first, then loose) or the current members on a
               direct visit; the same order Field is the escape hatch back to a
               proposed key. Commits run the SAME two-op scope->create chain via
               handleBuild. */
            <Card
              elevated
              data-testid="custom-scope-plate"
              className="w-full max-w-[760px] px-8 pt-7 pb-6"
            >
              <Kicker tone="accent">New series</Kicker>
              <h1 className="text-display text-ink mt-1.5">
                Series by {coldKey.trim() || "…"}
              </h1>
              <Kicker as="h2" tone="soft" className="mt-5 mb-2">
                Ordered by
              </Kicker>
              <Field
                testId="order-field"
                srLabel="Ordered by"
                menuLabel="Ordered by"
                value={DEFINE_YOUR_OWN}
                options={orderOptions}
                onSelect={onOrderSelect}
              />
              <div className="text-caption text-ink-soft mt-1.5">
                Name the variable below and assign each sample's value.
              </div>
              <div className="mt-6">
                {/* intro={null}: the default cold-corpus paragraph ("These
                    samples share no ordering variable yet...") would be false
                    here, and the caption above already instructs. */}
                <ColdAssignSection
                  rows={coldRows}
                  variableKey={coldKey}
                  onKeyChange={setColdKey}
                  onValueChange={onColdValueChange}
                  canBuildNow={canColdBuildNow}
                  onBuild={handleBuild}
                  buildBusy={buildBusy}
                  intro={null}
                  discardSlot={discardControl}
                />
              </div>
            </Card>
          ) : (
            <ScopePlate
              seriesName={`Series by ${keyLabel}`}
              grouping={
                <>
                  Himalaya grouped <strong>{rows.length} samples</strong> by their stored{" "}
                  <strong>{keyLabel}</strong> value. Confirm the order and build.
                  {/* SC-SKIPDISC: name the gesture, not just the possibility —
                      the skip toggle lives on the value itself.
                      SC-PROVREAD: the value is a stored tag (possibly
                      hand-entered), not necessarily parsed from the name. */}{" "}
                  Click a value to skip that sample, or the pencil to correct a misread.
                </>
              }
              orderedBy={keyLabel}
              orderOptions={orderOptions}
              onOrderSelect={onOrderSelect}
              orderNote={`Each value is the sample's stored ${keyLabel} tag. Switch the ordering variable above, or define your own.`}
              count={orderCaption}
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
                    data-reorder-row
                    data-reorder-index={i}
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
                      onMoveBy={(delta) => moveRow(i, delta)}
                      onEditValue={(v) => editValue(r.sampleId, v)}
                    />
                  </div>
                );
              })}
              candidates={
                loose.length ? (
                  /* SCOPE-LOOSEADD: loose matches lack a value for this key, so
                     each row now lets you NAME one and pull the sample in (the
                     same labelling gesture the cold path uses). The add is gated
                     on a non-empty value, so it never lies or seats a value:""
                     member. Capped to a few exemplars
                     (SCOPE_CANDIDATE_PREVIEW_COUNT); the section note below owns
                     the remainder count + the contact-sheet door for reaching the
                     hidden ones. */
                  <div data-testid="scope-candidates" className="space-y-2">
                    {visibleLoose.map((c) => (
                      <div
                        key={c.sampleId}
                        data-testid="scope-candidate"
                        className="flex items-center gap-3 px-2 py-2.5 border border-dashed border-hair-strong rounded"
                      >
                        <Sparkline
                          trace={traceBySample.get(c.sampleId) ?? EMPTY_TRACE}
                          {...(phaseBySample.get(c.sampleId) != null
                            ? { phase: phaseBySample.get(c.sampleId)! }
                            : {})}
                          className="opacity-70"
                        />
                        <div className="flex-1 min-w-0">
                          {/* Same identity vocabulary as the member rows
                              (ScopeSampleRow's mono smp_ caption) — display
                              names duplicate in the corpus, and without the id
                              a candidate can read as a member that is somehow
                              also out of the series. Inline (not stacked) so
                              the dimmed row stays two lines tall. */}
                          <div className="flex items-center gap-2 min-w-0">
                            {/* SC-CANDNAME-DIM: the name is the row's primary
                                identifier, so it reads at text-ink (out-ranking
                                its own smp_/why metadata). The candidate stays
                                visibly secondary via the smaller text-meta size,
                                the dashed border, and the dimmed sparkline, not
                                by dimming the name's colour to match its caption. */}
                            <span className="text-meta font-semibold text-ink">
                              {c.sampleName}
                            </span>
                            <span className="text-caption font-mono text-ink-soft">
                              smp_{c.sampleId}
                            </span>
                          </div>
                          <div className="text-caption text-ink-soft">
                            lacks the{" "}
                            <strong className="text-ink font-semibold">{keyLabel}</strong>
                          </div>
                        </div>
                        {/* SC-MISC (WCAG 1.4.1): the candidate's phase was a
                            colour-only signal on the aria-hidden sparkline.
                            Loose candidates never reach the preview strip that
                            labels member phases, so name the phase here as text
                            (PhaseChip's always-on second channel survives
                            grayscale + reaches a screen reader). Omitted when the
                            candidate is unindexed (no phase to name). */}
                        {phaseBySample.get(c.sampleId) != null && (
                          <PhaseChip
                            phase={phaseBySample.get(c.sampleId)!}
                            variant="tint"
                            size="sm"
                            className="flex-shrink-0"
                          />
                        )}
                        {/* SCOPE-LOOSEADD: name the missing value, then pull the
                            candidate into the series. Add is gated on a non-empty
                            value — a value-less add would seat a value:"" member
                            the build gate drops (controls-don't-lie). Enter in
                            the field is the same commit. */}
                        <Input
                          value={looseDraft[c.sampleId] ?? ""}
                          onValueChange={(v) =>
                            setLooseDraft((d) => ({ ...d, [c.sampleId]: v }))
                          }
                          placeholder={`${keyLabel}…`}
                          aria-label={`Value for ${c.sampleName}`}
                          inputSize="sm"
                          mono
                          className="w-24 flex-shrink-0"
                          onKeyDown={(e: React.KeyboardEvent) => {
                            if (e.key === "Enter") {
                              e.preventDefault();
                              addLoose(c.sampleId, looseDraft[c.sampleId] ?? "");
                            }
                          }}
                        />
                        <Button
                          variant="outline"
                          className="flex-shrink-0"
                          disabled={(looseDraft[c.sampleId] ?? "").trim() === ""}
                          onClick={() => addLoose(c.sampleId, looseDraft[c.sampleId] ?? "")}
                        >
                          + Add to series
                        </Button>
                      </div>
                    ))}
                    <div data-testid="scope-candidates-note" className="text-caption text-ink-soft pt-1">
                      {hiddenLooseCount > 0 ? (
                        <>
                          …and {hiddenLooseCount} more lack{hiddenLooseCount === 1 ? "s" : ""} the{" "}
                          <strong className="text-ink font-semibold">{keyLabel}</strong>.{" "}
                        </>
                      ) : null}
                      Tag a sample with the{" "}
                      <strong className="text-ink font-semibold">{keyLabel}</strong> on the
                      contact sheet if it belongs here.{" "}
                      <button
                        type="button"
                        onClick={() => navigate("/samples")}
                        className="text-caption font-semibold text-accent hover:underline"
                      >
                        Open the contact sheet
                      </button>
                    </div>
                  </div>
                ) : (
                  <div className="text-meta text-ink-soft italic">
                    Nothing else in the corpus matches this grouping.
                  </div>
                )
              }
              preview={preview}
              footState={footState}
              footNote={
                <>
                  Confirming records the {keyLabel} on every kept sample. The next series that needs it
                  already knows.
                </>
              }
              {...(!canBuild ? { buildDisabled: true } : {})}
              buildBusy={buildBusy}
              onBuild={handleBuild}
              discardSlot={discardControl}
            />
          )}
        </Skeleton>
      </div>
    </PageFrame>
  );
}
