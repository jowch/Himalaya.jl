import { useCallback, useEffect, useMemo, useRef, useState } from "react";
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
} from "../../queries";
import * as api from "../../api";
import { Button } from "../ui/Button";
import { ProgressBar } from "../ui/ProgressBar";
import { ScanFailedPage } from "./ScanFailedPage";
import { SheetTable } from "../components/SheetTable";
import { SampleTableRow } from "../components/SampleTableRow";
import { toSampleRowModel } from "./samplesAdapters";
import { navigateToNewSeries } from "../../lib/series/newSeriesNav";
import { showToast } from "../../lib/toast";
import { useListCursor } from "../interaction/useListCursor";
import { usePageActions } from "../interaction/usePageActions";
import { core, page } from "../interaction/core";
import { effectiveIngestStatus } from "../../lib/ingestStatus";

/** Exposure screening verdict applied by the cull verbs (null = unscreened). */
type Verdict = "accepted" | "rejected" | null;

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

  // --- ID-based sample cursor (roving tabindex, Enter → navigate to Focus) ---
  const sampleIds = useMemo(() => scopedSamples.map((s) => s.id), [scopedSamples]);
  const sampleCursor = useListCursor({
    ids: sampleIds,
    onActivate: (id) => navigate(`/sample/${id}`),
    stepperLabel: "Sample",
    stepperTestIdBase: "sample",
    axis: "vertical",
  });
  const activeSample = scopedSamples.find((s) => s.id === sampleCursor.cursorId);

  // Frame axis: page-local, reset when cursor changes
  const [frameIndex, setFrameIndex] = useState(0);
  useEffect(() => { setFrameIndex(0); }, [sampleCursor.cursorId]);

  // --- Selection state ---
  // SAMPLE-grain selection lives on the cursor (`sampleCursor.selected`) — one
  // selection state, driven by the row checkbox. The EXPOSURE-grain `selected`
  // Set below is a distinct page-local concern (which thumbnails are highlighted
  // within a row); it does not feed the dock's cull verbs.
  const [selected, setSelected] = useState<Set<number>>(() => new Set());
  const anchorRef = useRef<{ sampleId: number; exposureId: number } | null>(null);
  const shiftRef = useRef(false);

  // SA-STALESELECT: navigating between experiment ids re-renders this route
  // (same path, new :id) without remounting, so a working selection from the
  // prior experiment is no longer on screen. Clear the exposure grain on expId
  // change (the sample grain on `sampleCursor.selected` self-prunes when the id
  // set swaps); on mount this is a no-op.
  useEffect(() => {
    setSelected(new Set());
    anchorRef.current = null;
    setFrameIndex(0);
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



  // --- Sample-grain cull actions (declared to the interaction registry) ---
  // All three read the ONE sample selection (`sampleCursor.selected`) and apply
  // the verdict across every exposure of each selected sample. `toggleSelect()`
  // with no arg clears nothing; we clear by toggling each selected id off after
  // the batch so the selection resets (mirrors the old set-clear).
  const cullSelected = useCallback(
    (status: Verdict, verb: string) => {
      const ids = [...sampleCursor.selected];
      for (const sampleId of ids) {
        const exps = corpusExposures.byId.get(sampleId) ?? [];
        for (const e of exps) batch.mutate({ sampleId, exposureId: e.id, status });
      }
      const n = ids.length;
      if (n > 0) showToast(`${n} sample${n === 1 ? "" : "s"} ${verb}`, "success");
      for (const sampleId of ids) sampleCursor.toggleSelect(sampleId);
    },
    [sampleCursor, corpusExposures.byId, batch],
  );
  const dropSelected = useCallback(() => cullSelected("rejected", "dropped"), [cullSelected]);
  const keepSelected = useCallback(() => cullSelected("accepted", "kept"), [cullSelected]);
  const restoreSelected = useCallback(() => cullSelected(null, "restored"), [cullSelected]);

  const mode = sampleCursor.selected.size > 0 ? "selection" : "browse";

  usePageActions({
    cursor: sampleCursor,
    // ↑/↓ drive the sample cursor; ←/→ walk the active sample's frame axis.
    // Scope-exempt (the shell keyboard layer runs it) so arrows control the
    // surface wherever focus sits, instead of scrolling the page.
    arrowHandler: (e) => {
      if (e.key === "ArrowDown") { e.preventDefault(); sampleCursor.moveBy(1); }
      else if (e.key === "ArrowUp") { e.preventDefault(); sampleCursor.moveBy(-1); }
      else if (e.key === "ArrowLeft") { e.preventDefault(); setFrameIndex((i) => Math.max(0, i - 1)); }
      else if (e.key === "ArrowRight") { e.preventDefault(); setFrameIndex((i) => Math.min(Math.max(activeFrames.length - 1, 0), i + 1)); }
    },
    actions: [
      core("back", { label: "Experiments", run: () => navigate("/experiments"), dock: true }),
      core("openFocus", {
        run: () => sampleCursor.activate(),
        dock: "primary",
        enabled: () => sampleCursor.cursorId !== null,
      }),
      core("openLoupe", {
        run: () => { if (activeSample) navigate(`/sample/${activeSample.id}/loupe`); },
        dock: true,
        enabled: () => activeSample != null,
      }),
      page("cull", {
        label: "Drop",
        keys: ["x"],
        group: "Act",
        mode: "selection",
        enabled: () => mode === "selection",
        dock: true,
        run: () => dropSelected(),
      }),
      page("keep", {
        label: "Keep",
        keys: ["k"],
        group: "Act",
        mode: "selection",
        enabled: () => mode === "selection",
        dock: true,
        run: () => keepSelected(),
      }),
      page("restore", {
        label: "Restore",
        keys: ["r"],
        group: "Act",
        mode: "selection",
        enabled: () => mode === "selection",
        dock: true,
        run: () => restoreSelected(),
      }),
    ],
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

  // Active frame for cursoredExposureId highlight in the gallery.
  const activeFrames = activeSample ? (corpusExposures.byId.get(activeSample.id) ?? []) : [];
  const activeFrame = activeFrames[frameIndex];

  // ── Corpus sheet + cull/compose/dock ─────────────────────────────────────
  return (
    // pb-24 clears the fixed Dock (≈47px) so the last sample rows scroll above it
    // instead of hiding beneath — the sheet had no bottom clearance.
    <div data-testid="experiment-corpus" className="flex flex-col gap-4 pb-24">
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

      <div
        data-testid="corpus-grid"
        data-interaction-scope
      >
        {corpusQuery.isLoading ? (
          <div className="p-8 text-sm text-ink-soft">Loading samples…</div>
        ) : (
          <SheetTable checkboxColumn>
            {scopedSamples.map((s) => {
              const { ref: rowRef, ...rowRest } = sampleCursor.rowProps(s.id);
              const loadedExposures = corpusExposures.byId.get(s.id);
              const m = toSampleRowModel(s, loadedExposures);
              const noExposures = loadedExposures !== undefined && loadedExposures.length === 0;
              const hasDrop = (m.dropped ?? 0) > 0;
              return (
                <SampleTableRow
                  key={s.id}
                  ref={rowRef}
                  {...rowRest}
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
                  checked={sampleCursor.selected.has(s.id)}
                  onCheck={() => sampleCursor.toggleSelect(s.id)}
                  cursored={sampleCursor.cursorId === s.id}
                  {...(sampleCursor.cursorId === s.id && activeFrame
                    ? { cursoredExposureId: activeFrame.id }
                    : {})}
                  selectedExposureIds={selected}
                  onSelectExposure={(eid) => {
                    sampleCursor.setCursor(s.id);
                    const frames = corpusExposures.byId.get(s.id) ?? [];
                    const fi = frames.findIndex((e) => e.id === eid);
                    if (fi >= 0) setFrameIndex(fi);
                    toggleSelect(s.id, eid);
                  }}
                  onActivateExposure={(eid) => {
                    sampleCursor.setCursor(s.id);
                    navigate(loupeHref(s.id, eid), { state: { sampleOrder } });
                  }}
                  onOpenLoupe={() => {
                    sampleCursor.setCursor(s.id);
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
      </div>

      {/* Compose segment — appears when samples are selected. Reads the ONE
          sample selection (`sampleCursor.selected`) to build a new series. */}
      {sampleCursor.selected.size > 0 && (
        <div className="flex items-center gap-2 px-4 py-2" data-testid="dock-compose">
          <span className="text-meta text-ink-soft">
            {sampleCursor.selected.size} sample{sampleCursor.selected.size === 1 ? "" : "s"}
          </span>
          <Button variant="accent" data-testid="dock-new-series"
            onClick={() => navigateToNewSeries(sampleCursor.selected, navigate)}>+ New series</Button>
          <Button variant="ghost" data-testid="dock-clear-checks"
            onClick={() => { for (const sid of [...sampleCursor.selected]) sampleCursor.toggleSelect(sid); }}>Clear</Button>
        </div>
      )}
    </div>
  );
}
