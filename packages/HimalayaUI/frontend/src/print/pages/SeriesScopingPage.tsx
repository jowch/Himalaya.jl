import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { ReactNode } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { ScopePlate } from "../components/ScopePlate";
import { ScopeSampleRow } from "../components/ScopeSampleRow";
import { useDragReorder, reorder } from "../components/useDragReorder";
import { Sparkline } from "../plot/Sparkline";
import { EmptyState, Button, Card, Dot, Field, Kicker } from "../ui";
import { ColdAssignPanel } from "../components/ColdAssignPanel";
import type { Trace } from "../../api";
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
  filterPickerBySeed,
  buildColdAssignRows,
  canColdBuild,
  buildColdScopePayload,
  type ColdAssignRow,
} from "./scopingDerive";
import { readNewSeriesSeed } from "../../lib/series/newSeriesNav";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

// The sentinel "Define your own..." dropdown entry (locked decision A):
// selecting it inline-converts the warm worksheet to the ColdAssignPanel flow,
// seeded with the current members, committing through the SAME two-op
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

type HistoryEntry = { type: "flag"; id: number; prev: boolean; label: string };

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
  intro,
}: {
  rows: ColdAssignRow[];
  variableKey: string;
  onKeyChange: (key: string) => void;
  onValueChange: (sampleId: number, value: string) => void;
  canBuildNow: boolean;
  onBuild: () => void;
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
            Confirming records the variable on every sample — the next series that needs it
            already knows.
          </div>
        </div>
        <Button
          variant="solid"
          {...(!canBuildNow ? { disabled: true } : {})}
          onClick={onBuild}
        >
          Confirm &amp; build →
        </Button>
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
  // we scope the proposal to just those samples; when null (a direct visit) the
  // full corpus is used — the existing whole-corpus behaviour.
  const seed = useMemo(() => readNewSeriesSeed(location), [location]);
  const seededPicker = useMemo(
    () => filterPickerBySeed(pickerQ.data ?? [], seed),
    [pickerQ.data, seed],
  );

  // The corpus's real ordering-variable vocabulary: the distinct tag keys, in
  // proposeOrdering's deterministic order (frequency desc, lexicographic tie),
  // so the dropdown lists exactly what the corpus exposes — no declared schema.
  const orderKeys = useMemo(() => {
    const freq = new Map<string, number>();
    for (const t of tagsQ.data ?? []) freq.set(t.key, (freq.get(t.key) ?? 0) + 1);
    return [...freq.entries()]
      .sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1))
      .map(([k]) => k);
  }, [tagsQ.data]);

  // Human override of the machine's frequency winner. Null → the machine picks.
  const [overrideKey, setOverrideKey] = useState<string | null>(null);
  // A corpus change can't strand a stale override: drop it once it names a key
  // the corpus no longer exposes.
  useEffect(() => {
    if (overrideKey !== null && !orderKeys.includes(overrideKey)) setOverrideKey(null);
  }, [overrideKey, orderKeys]);

  const proposal = useMemo(
    () => proposeOrdering(tagsQ.data ?? [], seededPicker, overrideKey ?? undefined),
    [tagsQ.data, seededPicker, overrideKey],
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

  // ── Cold-corpus path (no proposable variable, user arrived with a seed) ──
  const [coldKey, setColdKey] = useState("");
  const [coldRows, setColdRows] = useState<ColdAssignRow[]>([]);
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
  // preview; never touches the written (key,value) payload.
  const { dragItemProps, dropEdge } = useDragReorder((from, to) =>
    setOrder((o) => reorder(o, from, to)),
  );

  const byId = useMemo(() => new Map(rows.map((r) => [r.sampleId, r])), [rows]);
  const sorted = useMemo(
    () => order.map((id) => byId.get(id)).filter((r): r is OrderingRow => r != null),
    [order, byId],
  );

  // Routing for the ordering dropdown. The sentinel enters custom mode, seeding
  // the assign rows from the current members in their displayed order (re-entry
  // deliberately re-seeds: previous custom edits are discarded). Any existing
  // label exits custom mode and applies the override (existing behaviour).
  const onOrderSelect = useCallback(
    (label: string): void => {
      if (label === DEFINE_YOUR_OWN) {
        setColdRows(
          buildColdAssignRows(
            sorted.map((r) => ({ sampleId: r.sampleId, sampleName: r.sampleName })),
          ),
        );
        setColdKey("");
        setCustomMode(true);
        return;
      }
      setCustomMode(false);
      setOverrideKey(labelToKey.get(label) ?? null);
    },
    [labelToKey, sorted],
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

  const undo = useCallback((): void => {
    setHistory((h) => {
      const e = h[h.length - 1];
      if (!e) return h;
      setRows((cur) => cur.map((r) => (r.sampleId === e.id ? { ...r, flagged: e.prev } : r)));
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

  const skippedCount = rows.filter((r) => r.flagged).length;
  const keptCount = rows.filter((r) => r.include && !r.flagged && r.value !== "").length;
  const footState = buildFootState(keptCount, skippedCount);
  const canBuild = canScopeBuild(rows, proposal.orderingKey);
  const lastLabel = history.length ? history[history.length - 1]!.label : undefined;

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
      // Warm members: the SCOPED, kept set (skips excluded — same predicate as
      // buildScopePayload) in the displayed low→high order, so the series
      // recipe matches the tags actually written.
      const keptInOrder = sorted.filter((r) => r.include && !r.flagged && r.value !== "");
      createBodyRef.current = {
        title: `Series by ${keyLabel}`,
        samples: keptInOrder.map((r, i) => ({ sample_id: r.sampleId, position: i })),
        ordering_variable: proposal.orderingKey,
      };
      stage.current = "tagging";
      scopeSeries.mutate({ key: proposal.orderingKey, tags: buildScopePayload(rows) });
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
            className="px-1 py-1.5 text-caption font-semibold text-ink-soft hover:text-ink"
          >
            Discard
          </button>
        </div>

        {/* State 4: a chain op failed. Two modes (see chainErrorCopy): an Op-A
            (tag write) failure means no series and no tags from this attempt; an
            Op-B (create) failure means the tags ARE committed — the copy admits
            it. Op A wins precedence because a stale create error can linger. */}
        {chainErrorCopy !== null ? (
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
                 deliberate selection, but no tag key can be proposed. Rather
                 than dead-end, let them name the ordering variable and assign
                 each sample's value — then commit through the SAME
                 scope→create chain the warm path uses. */
              <Card border="strong" padding="lg" data-testid="cold-scope-plate" className="w-full">
                <ColdAssignSection
                  rows={coldRows}
                  variableKey={coldKey}
                  onKeyChange={setColdKey}
                  onValueChange={onColdValueChange}
                  canBuildNow={canColdBuildNow}
                  onBuild={handleBuild}
                />
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
               the cold-assign flow, seeded with the current members; the same
               order Field is the escape hatch back to a proposed key. Commits
               run the SAME two-op scope->create chain via handleBuild. */
            <Card
              elevated
              data-testid="custom-scope-plate"
              className="w-full max-w-[760px] px-8 pt-7 pb-6"
            >
              <Kicker tone="accent">New series</Kicker>
              <h1 className="text-display text-ink mt-1.5">
                Series by {coldKey.trim() || "…"}
              </h1>
              <Kicker as="h2" tone="faint" className="mt-5 mb-2">
                Ordered by
              </Kicker>
              <Field
                testId="order-field"
                srLabel="Ordered by"
                value={DEFINE_YOUR_OWN}
                options={orderOptions}
                onSelect={onOrderSelect}
              />
              <div className="text-caption text-ink-soft mt-1.5">
                Name the variable below and assign each sample's value.
              </div>
              <div className="mt-6">
                {/* intro={null}: the default cold-corpus paragraph ("These
                    samples share no tag key yet...") would be false here, and
                    the caption above already instructs. */}
                <ColdAssignSection
                  rows={coldRows}
                  variableKey={coldKey}
                  onKeyChange={setColdKey}
                  onValueChange={onColdValueChange}
                  canBuildNow={canColdBuildNow}
                  onBuild={handleBuild}
                  intro={null}
                />
              </div>
            </Card>
          ) : (
            <ScopePlate
              seriesName={`Series by ${keyLabel}`}
              grouping={
                <>
                  Himalaya grouped <strong>{rows.length} samples</strong> by their{" "}
                  <strong>{keyLabel}</strong>, read from the sample names. Confirm the reads and build —
                  skip any it misread.
                </>
              }
              orderedBy={keyLabel}
              orderOptions={orderOptions}
              onOrderSelect={onOrderSelect}
              orderNote="Read from the sample names. Switch the ordering variable above, or define your own."
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
                  /* Informational discovery only (Option A): loose matches carry
                     no value for this key, so they CANNOT be committed — there is
                     no add action. A plain sparkline + name + why list, not the
                     ScopeCandidateRow composite (whose "+ Add to series" button
                     always renders; a dead button would lie). Capped to a few
                     exemplars (SCOPE_CANDIDATE_PREVIEW_COUNT); the section note
                     below owns the remainder count + the tag instruction once,
                     instead of repeating it per row. */
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
                          <div className="text-meta font-semibold text-ink-soft">{c.sampleName}</div>
                          <div className="text-caption text-ink-soft">
                            lacks the{" "}
                            <strong className="text-accent font-semibold">{keyLabel}</strong>
                          </div>
                        </div>
                      </div>
                    ))}
                    <div data-testid="scope-candidates-note" className="text-caption text-ink-soft pt-1">
                      {hiddenLooseCount > 0 ? (
                        <>
                          …and {hiddenLooseCount} more lack{hiddenLooseCount === 1 ? "s" : ""} the{" "}
                          <strong className="text-accent font-semibold">{keyLabel}</strong>.{" "}
                        </>
                      ) : null}
                      Tag a sample with the{" "}
                      <strong className="text-accent font-semibold">{keyLabel}</strong> on the
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
                  Confirming records the {keyLabel} on every kept sample — the next series that needs it
                  already knows.
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
