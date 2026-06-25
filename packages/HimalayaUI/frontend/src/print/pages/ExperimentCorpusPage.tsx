import { useEffect, useMemo, useRef, useState } from "react";
import { Link, Navigate, useNavigate, useParams } from "react-router-dom";
import { useQuery } from "@tanstack/react-query";
import { useAppState } from "../../state";
import {
  useExperiment,
  useLoads,
  useTriggerScan,
  useCorpusSamples,
  useCorpusExposures,
  useSetExposureStatusBatch,
  useSelectExposure,
} from "../../queries";
import * as api from "../../api";
import { Button } from "../ui/Button";
import { IconButton } from "../ui/IconButton";
import { KbKey } from "../ui/KbKey";
import { ProgressBar } from "../ui/ProgressBar";
import { Dock } from "../ui/Dock";
import { DockUpLink } from "../ui/DockUpLink";
import { ScanFailedPage } from "./ScanFailedPage";
import { SheetTable } from "../components/SheetTable";
import { SampleTableRow } from "../components/SampleTableRow";
import { toSampleRowModel } from "./samplesAdapters";
import { navigateToNewSeries } from "../../lib/series/newSeriesNav";
import { useShortcuts } from "../shell/useShortcuts";
import { isNativeInteractiveTarget } from "../../lib/keys";
import { showToast } from "../../lib/toast";
import { effectiveIngestStatus } from "../../lib/ingestStatus";

/** Distinct samples a cull selection spans, counted EXACTLY the way the
 *  Drop/Keep/Restore batch routes it (unmappable ids are skipped). The CullBar
 *  disclosure and the toast receipt both go through this so the promise and the
 *  action can never count differently. (Ported from SamplesPage.) */
function selectionSpread(
  selected: ReadonlySet<number>,
  ownerOf: ReadonlyMap<number, number>,
): number {
  const owners = new Set<number>();
  for (const id of selected) {
    const sampleId = ownerOf.get(id);
    if (sampleId !== undefined) owners.add(sampleId);
  }
  return owners.size;
}

type Verdict = "accepted" | "rejected" | null;

/** Drop/Keep are TOGGLES: applying a verdict a frame already carries clears it
 *  to unscreened (null). Restore (target null) is never a toggle. */
function verdictNext(current: Verdict, target: Verdict): Verdict {
  return target !== null && current === target ? null : target;
}

/**
 * ExperimentCorpusPage — the experiment's Corpus home (the index route under
 * /experiments/:id, app-shell-unification IA). The "daily-loop" surface: a
 * contact sheet of this experiment's samples with a cull/keep workflow, scoped
 * to one experiment.
 *
 * §6.2 state machine (takeover states early-return before the sheet/dock):
 *   - scanning   (inFlight.status==="scanning")      → GroupingReviewPage (live-unfold)
 *   - rescanning (inFlight.status==="analyzing")     → inline ProgressBar
 *   - failed     (!processing && status==="failed")  → ScanFailedPage  [T4.2]
 *   - otherwise  → grouping-review banner (when flags) + real SheetTable + Dock
 *
 * Two data sources, joined only by intent:
 *   - useCorpusSamples() filtered to experiment_id===expId → the rich rows
 *     (thumbnails / kept / phase via useCorpusExposures + toSampleRowModel).
 *   - useLoads(expId) → the grouping-flag count that drives the review banner.
 */
export function ExperimentCorpusPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const navigate = useNavigate();

  const inFlight = useAppState((s) => s.ingestInFlight?.[expId]);
  const loads = useLoads(expId);
  const exp = useExperiment(expId);
  const triggerScan = useTriggerScan(expId);

  // --- State machine derivation ---
  // Effective status reconciles the live overlay with the persisted truth so a
  // dropped SSE terminal frame (8c) self-heals once the polled row turns terminal.
  // The surface CHOICE (initial-scan grouping takeover vs inline rescan progress)
  // can't ride on persisted "scanning" alone — both routes set it — so it keys on
  // last_scanned_at, which is null only until the first scan completes. That keeps
  // a reload mid-rescan (overlay lost, persisted "scanning") on the inline path
  // instead of wrongly re-entering the full grouping takeover.
  const eff = effectiveIngestStatus(inFlight?.status, exp.data?.ingest_status);
  const processing = eff === "scanning" || eff === "analyzing";
  // The overlay phase is authoritative when present: "analyzing" ⇒ a rescan
  // (inline ProgressBar), "scanning" ⇒ an initial scan (grouping takeover). Only
  // when the overlay is absent (a reload mid-scan) does persisted "scanning"
  // remain, which BOTH routes set — disambiguate with last_scanned_at, null only
  // until the first scan completes, so a reload mid-rescan stays on the inline
  // path instead of wrongly re-entering the full takeover.
  const isInitialScan = eff === "scanning" && exp.data?.last_scanned_at == null;
  const scanning = processing && isInitialScan;
  const rescanning = processing && !isInitialScan;
  const failed = !processing && eff === "failed";

  // Manifest query — unconditional (hooks rule), enabled only in the failed
  // branch when data_dir is known. Provides real unmatched + parsedCount for
  // ScanFailedPage. exactOptionalPropertyTypes: conditional spread keeps absent
  // keys absent so the backend falls back to its default globs.
  const manifestPatterns = useMemo(() => ({
    ...(exp.data?.image_pattern != null && { image: exp.data.image_pattern }),
    ...(exp.data?.metadata_pattern != null && { metadata: exp.data.metadata_pattern }),
    ...(exp.data?.integration_pattern != null && { integration: exp.data.integration_pattern }),
  }), [exp.data?.image_pattern, exp.data?.metadata_pattern, exp.data?.integration_pattern]);

  // B1 (§5.5): pass the stored leaf analysis_dir so integration (.dat) is matched
  // against the analysis subtree where it actually lives. Without it, .dat is
  // matched against data_dir (where it never is) and every integration trace
  // reports `unmatched`, so the pattern-test can never clear. The experiment's
  // analysis_dir was persisted from the resolved leaf at create time.
  const manifestQuery = useQuery({
    queryKey: ["manifest", exp.data?.data_dir ?? "", manifestPatterns, exp.data?.analysis_dir ?? ""],
    queryFn: () => api.fetchManifest(
      exp.data!.data_dir, manifestPatterns, undefined, exp.data!.analysis_dir ?? undefined,
    ),
    enabled: failed && !!exp.data?.data_dir,
  });

  // --- Grouping-review banner: flag count + breakdown from the loads roll-up ---
  const loadSamples = useMemo(
    () => (loads.data ?? []).flatMap((l) => l.samples),
    [loads.data],
  );
  const flagged = useMemo(
    () => loadSamples.filter((s) => s.flag !== null),
    [loadSamples],
  );
  const reviewCount = flagged.length;
  // Honest breakdown of WHY review is needed, mirroring the mockup's structure
  // (a "stage-position jump" is a split flag's jump_from→jump_to; a "counter
  // reset" is a merge flag). Only non-zero parts appear.
  // slotBySample: map from sample_id → slot_index for the slot chip in each row.
  const slotBySample = useMemo(() => {
    const m = new Map<number, number>();
    for (const s of loadSamples) m.set(s.sample_id, s.slot_index);
    return m;
  }, [loadSamples]);

  const reviewDetail = useMemo(() => {
    const splits = flagged.filter((s) => s.flag?.kind === "split").length;
    const merges = flagged.filter((s) => s.flag?.kind === "merge").length;
    const parts: string[] = [];
    if (splits > 0) parts.push(`${splits} stage-position ${splits === 1 ? "jump" : "jumps"}`);
    if (merges > 0) parts.push(`${merges} counter ${merges === 1 ? "reset" : "resets"}`);
    if (parts.length === 0) return "";
    const joined = parts.join(" and ");
    const verb = splits + merges === 1 ? "was" : "were";
    return `${joined[0]!.toUpperCase()}${joined.slice(1)} ${verb} ambiguous.`;
  }, [flagged]);

  // --- Rich rows: this experiment's corpus samples + their exposures ---
  const corpusQuery = useCorpusSamples();
  const scopedSamples = useMemo(
    () => (corpusQuery.data ?? []).filter((s) => s.experiment_id === expId),
    [corpusQuery.data, expId],
  );
  const corpusExposures = useCorpusExposures(scopedSamples);
  const batch = useSetExposureStatusBatch();

  // --- Page cursor (active row + detail frame), driven by ↑/↓/←/→ ---
  const [cursor, setCursor] = useState<{ sampleIndex: number; frameIndex: number }>(
    { sampleIndex: 0, frameIndex: 0 },
  );
  const activeSample = scopedSamples[cursor.sampleIndex];
  // Keep the roving ↑/↓ row visible: when the cursor row changes, pull it into
  // the shell scroller's viewport. `block: "nearest"` only moves when the row is
  // off-screen and aligns to the closest edge (no jumpy re-centering). The row
  // carries `data-cursored="true"`; scope the lookup to this page's subtree.
  const corpusRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    // `?.scrollIntoView?.()` — guard the call: jsdom (unit tests) doesn't
    // implement scrollIntoView, so optional-call it instead of throwing.
    corpusRef.current
      ?.querySelector('[data-cursored="true"]')
      ?.scrollIntoView?.({ block: "nearest" });
  }, [cursor.sampleIndex]);
  // Representative setter, keyed to the cursor's active sample. Lets "R" / the
  // dock "Mark rep" button flag which exposure represents the sample (the one
  // Focus opens) directly from the contact sheet.
  const setRepresentative = useSelectExposure(activeSample?.id ?? 0);

  function clamp(v: number, lo: number, hi: number): number {
    return v < lo ? lo : v > hi ? hi : v;
  }
  function clampSample(d: number): void {
    setCursor((c) => ({
      sampleIndex: clamp(c.sampleIndex + d, 0, Math.max(0, scopedSamples.length - 1)),
      frameIndex: 0,
    }));
  }
  function clampFrame(d: number): void {
    setCursor((c) => {
      const frames = corpusExposures.byId.get(scopedSamples[c.sampleIndex]?.id ?? -1) ?? [];
      return { ...c, frameIndex: clamp(c.frameIndex + d, 0, Math.max(0, frames.length - 1)) };
    });
  }

  // --- Selection state: exposure-grain cull + sample-grain pick ---
  const [selected, setSelected] = useState<Set<number>>(() => new Set());
  const anchorRef = useRef<{ sampleId: number; exposureId: number } | null>(null);
  const shiftRef = useRef(false);
  const [checkedSamples, setCheckedSamples] = useState<Set<number>>(() => new Set());

  function toggleSampleCheck(sampleId: number): void {
    setCheckedSamples((prev) => {
      const next = new Set(prev);
      if (next.has(sampleId)) next.delete(sampleId);
      else next.add(sampleId);
      return next;
    });
  }

  // SA-STALESELECT: navigating between experiment ids re-renders this route
  // (same path, new :id) without remounting, so a working selection from the
  // prior experiment is no longer on screen. Clear both grains on expId change;
  // on mount this is a no-op (both sets start empty).
  useEffect(() => {
    setSelected(new Set());
    setCheckedSamples(new Set());
    anchorRef.current = null;
    setCursor({ sampleIndex: 0, frameIndex: 0 });
  }, [expId]);

  // Shift-anchor tracking for contiguous range selection.
  useEffect(() => {
    function onShiftDown(e: KeyboardEvent): void { if (e.key === "Shift") shiftRef.current = true; }
    function onShiftUp(e: KeyboardEvent): void { if (e.key === "Shift") shiftRef.current = false; }
    window.addEventListener("keydown", onShiftDown);
    window.addEventListener("keyup", onShiftUp);
    return () => {
      window.removeEventListener("keydown", onShiftDown);
      window.removeEventListener("keyup", onShiftUp);
    };
  }, []);

  // exposureId → owning sampleId (route a batch mutate back to its sample).
  const ownerOf = useMemo(() => {
    const m = new Map<number, number>();
    for (const s of scopedSamples)
      for (const e of corpusExposures.byId.get(s.id) ?? []) m.set(e.id, s.id);
    return m;
  }, [scopedSamples, corpusExposures.byId]);

  function toggleSelect(sampleId: number, exposureId: number): void {
    setSelected((prev) => {
      const next = new Set(prev);
      const anchor = anchorRef.current;
      if (shiftRef.current && anchor && anchor.sampleId === sampleId) {
        const ids = (corpusExposures.byId.get(sampleId) ?? []).map((e) => e.id);
        const a = ids.indexOf(anchor.exposureId);
        const b = ids.indexOf(exposureId);
        if (a >= 0 && b >= 0) {
          const lo = Math.min(a, b);
          const hi = Math.max(a, b);
          for (let i = lo; i <= hi; i++) next.add(ids[i]!);
          return next;
        }
      }
      if (next.has(exposureId)) next.delete(exposureId);
      else next.add(exposureId);
      anchorRef.current = { sampleId, exposureId };
      return next;
    });
  }

  function batchSet(status: Verdict): void {
    const targets = [...selected].filter((id) => ownerOf.has(id));
    const statusOf = (id: number): Verdict => {
      const sId = ownerOf.get(id)!;
      return corpusExposures.byId.get(sId)?.find((e) => e.id === id)?.status ?? null;
    };
    // Toggle as a set: if every selected frame already carries this verdict,
    // clear them all to unscreened instead (mirrors the single-frame toggle).
    const next: Verdict =
      status !== null && targets.length > 0 && targets.every((id) => statusOf(id) === status)
        ? null
        : status;
    for (const id of targets) {
      batch.mutate({ sampleId: ownerOf.get(id)!, exposureId: id, status: next });
    }
    const n = targets.length;
    if (n > 0) {
      const verb = next === "rejected" ? "dropped" : next === "accepted" ? "kept" : "restored";
      const spread = selectionSpread(new Set(targets), ownerOf);
      const suffix = spread > 1 ? ` across ${spread} samples` : "";
      showToast(`${n} frame${n === 1 ? "" : "s"} ${verb}${suffix}`, "success");
    }
    setSelected(new Set());
    anchorRef.current = null;
  }

  // --- Keyboard map (shared registry): cursor + selection + cull verbs ---
  useShortcuts({
    prevSample: () => clampSample(-1),
    nextSample: () => clampSample(1),
    prevDetail: () => clampFrame(-1),
    nextDetail: () => clampFrame(1),
    openFocus: (e) => {
      // §8 invariant (b): on a native interactive target (button/link/sort
      // header) Enter activates that control natively — decline.
      if (isNativeInteractiveTarget(e)) return false;
      if (activeSample == null) return false;
      navigate(`/sample/${activeSample.id}`);
      return undefined;
    },
    openLoupe: () => {
      if (activeSample == null) return false;
      navigate(`/sample/${activeSample.id}/loupe`);
      return undefined;
    },
    toggleSelect: () => {
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      toggleSelect(s.id, frame.id);
      return undefined;
    },
    extendPrev: () => {
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      if (cursor.frameIndex > 0) {
        const frame = frames[cursor.frameIndex - 1];
        if (s != null && frame != null) toggleSelect(s.id, frame.id);
      }
      return undefined;
    },
    extendNext: () => {
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex + 1];
      if (s != null && frame != null) toggleSelect(s.id, frame.id);
      return undefined;
    },
    selectAll: () => {
      const allFrameIds: number[] = [];
      for (const sam of scopedSamples)
        for (const f of corpusExposures.byId.get(sam.id) ?? []) allFrameIds.push(f.id);
      if (allFrameIds.length === 0) return false;
      setSelected(new Set(allFrameIds));
      return undefined;
    },
    drop: () => {
      if (selected.size > 0) { batchSet("rejected"); return undefined; }
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      const next = verdictNext(frame.status ?? null, "rejected");
      batch.mutate({ sampleId: s.id, exposureId: frame.id, status: next });
      showToast(`1 frame ${next === "rejected" ? "dropped" : "restored"}`, "success");
      return undefined;
    },
    keep: () => {
      if (selected.size > 0) { batchSet("accepted"); return undefined; }
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      const next = verdictNext(frame.status ?? null, "accepted");
      batch.mutate({ sampleId: s.id, exposureId: frame.id, status: next });
      showToast(`1 frame ${next === "accepted" ? "kept" : "restored"}`, "success");
      return undefined;
    },
    restore: () => {
      if (selected.size > 0) { batchSet(null); return undefined; }
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      batch.mutate({ sampleId: s.id, exposureId: frame.id, status: null });
      showToast("1 frame restored", "success");
      return undefined;
    },
    dismiss: () => {
      if (selected.size === 0) return false;
      setSelected(new Set());
      anchorRef.current = null;
      return undefined;
    },
    representative: () => {
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      if (frame.selected) return undefined; // already the representative
      setRepresentative.mutate(frame.id, {
        onSuccess: () => showToast("Set as the representative frame", "success"),
      });
      return undefined;
    },
  });

  // LO-NEXT: hand the loupe the exact visible order so its prev/next walk
  // matches what's on screen.
  const sampleOrder = useMemo(() => scopedSamples.map((s) => s.id), [scopedSamples]);
  function loupeHref(sampleId: number, exposureId?: number): string {
    const params = new URLSearchParams({ experiment: String(expId) });
    if (exposureId !== undefined) params.set("exposure", String(exposureId));
    return `/sample/${sampleId}/loupe?${params.toString()}`;
  }

  // ── Takeover states (no sheet/dock) ──────────────────────────────────────
  if (scanning) {
    // The scanning/grouping surface is a full takeover that owns its own header
    // (p1-grouping), so it lives at /grouping OUTSIDE this experiment shell.
    // Landing on the corpus mid-scan redirects there rather than nesting it.
    return <Navigate to={`/experiments/${expId}/grouping`} replace />;
  }
  // §8 (item 8): a rescan no longer swaps the whole page to an "Analyzing…"
  // takeover (which unmounted + remounted the sheet on entry/exit — the
  // "flash"). It now renders as an inline banner ABOVE the still-mounted sheet
  // (see the rescanning banner in the return below), so the rows stay put while
  // the additive scan runs.
  if (failed) {
    return (
      <ScanFailedPage
        experimentId={expId}
        unmatched={manifestQuery.data?.unmatched ?? []}
        parsedCount={manifestQuery.data?.matched.image ?? 0}
        onIngestParsed={() => triggerScan.mutate(true)}
        dataDir={exp.data?.data_dir ?? ""}
        {...(exp.data?.analysis_dir != null ? { analysisDir: exp.data.analysis_dir } : {})}
        patterns={manifestPatterns}
      />
    );
  }

  // Dock readouts + the cull verbs' target (the cursor's active frame). Derived
  // once so the steppers' "N / M" counts and the Drop/Keep/Restore buttons read
  // the SAME cursor truth.
  const activeFrames = activeSample ? (corpusExposures.byId.get(activeSample.id) ?? []) : [];
  const activeFrame = activeFrames[cursor.frameIndex];
  const sampleTotal = scopedSamples.length;
  const samplePos = sampleTotal > 0 ? cursor.sampleIndex + 1 : 0;
  const frameTotal = activeFrames.length;
  const framePos = frameTotal > 0 ? cursor.frameIndex + 1 : 0;
  const cullActiveFrame = (status: Verdict): void => {
    if (activeSample == null || activeFrame == null) return;
    const next = verdictNext(activeFrame.status ?? null, status);
    batch.mutate({ sampleId: activeSample.id, exposureId: activeFrame.id, status: next });
    const verb = next === "rejected" ? "dropped" : next === "accepted" ? "kept" : "restored";
    showToast(`1 frame ${verb}`, "success");
  };
  const markRepresentative = (): void => {
    if (activeSample == null || activeFrame == null || activeFrame.selected) return;
    setRepresentative.mutate(activeFrame.id, {
      onSuccess: () => showToast("Set as the representative frame", "success"),
    });
  };

  // Dock-as-action-bar (items 2/5): the floating CullBar/ComposeBar are gone.
  // The dock's verbs act on the multi-frame cull selection when one exists,
  // otherwise on the cursor's active frame; the compose segment appears when
  // samples are checked. `selSpread` discloses how many samples a frame
  // selection spans (the same count the batch routes through).
  const hasSel = selected.size > 0;
  const hasChecks = checkedSamples.size > 0;
  const selSpread = selectionSpread(selected, ownerOf);
  const clearSelection = (): void => { setSelected(new Set()); anchorRef.current = null; };

  // ── Corpus sheet + cull/compose/dock ─────────────────────────────────────
  return (
    // pb-24 clears the fixed Dock (≈47px) so the last sample rows scroll above it
    // instead of hiding beneath — the sheet had no bottom clearance.
    <div ref={corpusRef} data-testid="experiment-corpus" className="flex flex-col gap-4 pb-24">
      {/* Inline rescan banner (item 8): an additive rescan reports progress here
          without unmounting the sheet/dock below — the rows stay put. */}
      {rescanning && (
        <div
          data-testid="live-ingest-slot"
          className="flex flex-col gap-2 rounded-sm border border-hair bg-paper-sunk px-4 py-3"
        >
          <p className="text-sm text-ink-soft">Analyzing exposures…</p>
          <ProgressBar
            value={inFlight ? inFlight.processed : 0}
            total={inFlight ? Math.max(inFlight.total, 1) : 1}
            label="Analysis progress"
          />
        </div>
      )}

      {/* Grouping-review banner — always present once loads are known. When
          samples need a check it's an amber call-to-action with the breakdown;
          otherwise a calm "settled, but review if something's off" that keeps the
          SAME entry point so the grouping review is always one click away. */}
      {loads.data !== undefined && (
        <div
          data-testid="grouping-review-banner"
          data-state={reviewCount > 0 ? "attention" : "clear"}
          className={`flex items-center gap-3 rounded-sm border px-4 py-3 ${
            reviewCount > 0 ? "border-warning bg-warning/10" : "border-hair bg-paper-sunk"
          }`}
        >
          <span aria-hidden className={reviewCount > 0 ? "text-warning" : "text-ink-faint"}>
            {reviewCount > 0 ? "⚠" : "✓"}
          </span>
          <p className="text-sm text-ink-soft">
            {reviewCount > 0 ? (
              <>
                <span className="font-semibold text-ink">
                  {reviewCount} {reviewCount === 1 ? "sample needs" : "samples need"} a grouping check.
                </span>{" "}
                {reviewDetail}
              </>
            ) : (
              <>
                <span className="font-semibold text-ink">Grouping looks settled.</span>{" "}
                Open the review if a sample looks split, merged, or is missing exposures.
              </>
            )}
          </p>
          <Link
            to={`/experiments/${expId}/grouping`}
            data-testid="grouping-review-link"
            className="ml-auto whitespace-nowrap text-sm font-semibold text-print-accent hover:underline"
          >
            Review grouping →
          </Link>
        </div>
      )}

      {corpusQuery.isLoading ? (
        <div className="p-8 text-sm text-ink-soft">Loading samples…</div>
      ) : (
        <SheetTable checkboxColumn>
          {scopedSamples.map((s, rowIndex) => {
            const loadedExposures = corpusExposures.byId.get(s.id);
            const m = toSampleRowModel(s, loadedExposures);
            const noExposures = loadedExposures !== undefined && loadedExposures.length === 0;
            const hasDrop = (m.dropped ?? 0) > 0;
            return (
              <SampleTableRow
                key={s.id}
                name={m.name}
                sampleId={m.sampleId}
                screened={m.screened}
                exposures={m.exposures}
                kept={m.kept}
                total={m.total}
                dropped={m.dropped}
                noExposures={noExposures}
                tags={m.tags}
                {...(m.phase !== undefined ? { phase: m.phase } : {})}
                {...(m.formFactor ? { formFactor: true } : {})}
                {...(slotBySample.has(s.id) ? { slotIndex: slotBySample.get(s.id)! } : {})}
                checked={checkedSamples.has(s.id)}
                onCheck={() => toggleSampleCheck(s.id)}
                cursored={rowIndex === cursor.sampleIndex}
                {...(rowIndex === cursor.sampleIndex && activeFrame
                  ? { cursoredExposureId: activeFrame.id }
                  : {})}
                selectedExposureIds={selected}
                onSelectExposure={(eid) => {
                  setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                  toggleSelect(s.id, eid);
                }}
                onActivateExposure={(eid) => {
                  setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                  navigate(loupeHref(s.id, eid), { state: { sampleOrder } });
                }}
                onOpenLoupe={() => {
                  setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                  navigate(loupeHref(s.id), { state: { sampleOrder } });
                }}
                {...(hasDrop ? { onRestore: () => {
                  const drops = (corpusExposures.byId.get(s.id) ?? [])
                    .filter((e) => e.status === "rejected")
                    .map((e) => e.id);
                  if (drops.length === 0) return;
                  for (const did of drops) batch.mutate({ sampleId: s.id, exposureId: did, status: null });
                  showToast(`${drops.length} frame${drops.length === 1 ? "" : "s"} restored`, "success");
                } } : {})}
              />
            );
          })}
        </SheetTable>
      )}

      {/* Contextual bottom dock (Corpus grammar §3.3, mockup b1):
          ‹ Experiments │ Sample ‹N/M› │ Frame ‹N/M› │ Drop[X] Keep[K] Restore ──→ Loupe[L] Focus
          Segments are grouped (each a flex child of the Dock's gap-2 row); a
          flex-1 spacer right-anchors the destinations. */}
      <Dock>
        <DockUpLink label="Experiments" href="/experiments" onClick={() => navigate("/experiments")} />

        <span className="mx-1 h-6 w-px bg-hair" aria-hidden />

        {/* Sample stepper — ↑/↓ axis, current / total readout (pages2 dock: faint
            label, bordered ↑/↓ boxes, faint total). */}
        <div className="flex items-center gap-1.5">
          <span className="mr-0.5 text-meta text-ink-faint">Sample</span>
          <IconButton label="Previous sample" tone="ghost" boxed disabled={cursor.sampleIndex === 0}
            onClick={() => clampSample(-1)} data-testid="dock-prev-sample">↑</IconButton>
          <span className="text-data tabular-nums text-ink text-center min-w-[3.5rem]"
            data-testid="dock-sample-count">{samplePos}<span className="font-normal text-ink-faint"> / {sampleTotal}</span></span>
          <IconButton label="Next sample" tone="ghost" boxed disabled={cursor.sampleIndex >= scopedSamples.length - 1}
            onClick={() => clampSample(1)} data-testid="dock-next-sample">↓</IconButton>
        </div>

        <span className="mx-1 h-6 w-px bg-hair" aria-hidden />

        {/* Frame stepper — ←/→ axis within the active sample */}
        <div className="flex items-center gap-1.5">
          <span className="mr-0.5 text-meta text-ink-faint">Frame</span>
          <IconButton label="Previous frame" tone="ghost" boxed disabled={cursor.frameIndex === 0}
            onClick={() => clampFrame(-1)} data-testid="dock-prev-frame">←</IconButton>
          <span className="text-data tabular-nums text-ink text-center min-w-[2.75rem]"
            data-testid="dock-frame-count">{framePos}<span className="font-normal text-ink-faint"> / {frameTotal}</span></span>
          <IconButton label="Next frame" tone="ghost" boxed
            disabled={cursor.frameIndex >= frameTotal - 1}
            onClick={() => clampFrame(1)} data-testid="dock-next-frame">→</IconButton>
        </div>

        <span className="mx-1 h-6 w-px bg-hair" aria-hidden />

        {/* Cull verbs — act on the multi-frame selection when one exists (the
            old floating CullBar's job, now folded into the dock, items 2/5),
            otherwise on the cursor's active frame. A selection readout + Clear
            appear only while a selection is live. */}
        <div className="flex items-center gap-1">
          {hasSel && (
            <span className="mr-1 inline-flex items-center gap-2 text-meta font-semibold text-ink" data-testid="dock-selection-count">
              <span className="h-2 w-2 rounded-sm bg-accent" aria-hidden />
              {selected.size} frame{selected.size === 1 ? "" : "s"}
              {selSpread > 1 ? ` · ${selSpread} samples` : ""}
            </span>
          )}
          <Button variant="outlineAccent" data-testid="dock-drop"
            onClick={() => (hasSel ? batchSet("rejected") : cullActiveFrame("rejected"))}>Drop<KbKey className="ml-1.5">X</KbKey></Button>
          <Button variant="outlineSuccess" data-testid="dock-keep"
            onClick={() => (hasSel ? batchSet("accepted") : cullActiveFrame("accepted"))}>Keep<KbKey className="ml-1.5">K</KbKey></Button>
          <Button variant="ghost" data-testid="dock-restore"
            onClick={() => (hasSel ? batchSet(null) : cullActiveFrame(null))}>Restore</Button>
          {hasSel && (
            <Button variant="ghost" data-testid="dock-clear-selection"
              onClick={clearSelection}>Clear<KbKey className="ml-1.5">Esc</KbKey></Button>
          )}
        </div>

        {/* Spacer — right-anchors the destinations */}
        <div className="flex-1" />

        {/* Compose segment — appears when samples are checked (the old floating
            ComposeBar's job, now in the dock). Sample-grain, distinct from the
            frame-grain cull selection above. */}
        {hasChecks && (
          <>
            <div className="flex items-center gap-1" data-testid="dock-compose">
              <span className="text-meta text-ink-soft">
                {checkedSamples.size} sample{checkedSamples.size === 1 ? "" : "s"}
              </span>
              <Button variant="accent" data-testid="dock-new-series"
                onClick={() => navigateToNewSeries(checkedSamples, navigate)}>+ New series</Button>
              <Button variant="ghost" data-testid="dock-clear-checks"
                onClick={() => setCheckedSamples(new Set())}>Clear</Button>
            </div>
            <span className="mx-1 h-6 w-px bg-hair" aria-hidden />
          </>
        )}

        {/* Destinations — Mark rep (flags the active frame as the sample's
            representative, the one Focus opens) + Loupe + Focus (the primary) */}
        <div className="flex items-center gap-1">
          <Button variant="ghost" data-testid="dock-mark-rep"
            disabled={activeFrame == null || (activeFrame.selected ?? false)}
            onClick={markRepresentative}
          >Mark rep<KbKey className="ml-1.5">R</KbKey></Button>
          <Button variant="outline" data-testid="dock-loupe"
            onClick={() => { if (activeSample == null) return; navigate(`/sample/${activeSample.id}/loupe`); }}
          >Loupe<KbKey className="ml-1.5">L</KbKey></Button>
          <Button variant="accent" data-testid="dock-focus"
            onClick={() => { if (activeSample == null) return; navigate(`/sample/${activeSample.id}`); }}
          >Focus<KbKey variant="frost" className="ml-1.5">↵</KbKey></Button>
        </div>
      </Dock>
    </div>
  );
}
