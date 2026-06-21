import { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { SheetTable } from "../components/SheetTable";
import { SampleTableRow } from "../components/SampleTableRow";
import { CullBar } from "../components/CullBar";
import { Kicker, KbLegend, ProgressBar, ComposeBar, EmptyState, Button, HintText, IconButton } from "../ui";
import { Dock } from "../ui/Dock";
import { KbdLegend } from "../shell/KbdLegend";
import {
  useCorpusSamples,
  useCorpusExposures,
  useScreenedProgress,
  useExperiments,
  useSetExposureStatusBatch,
} from "../../queries";
import { navigateToNewSeries } from "../../lib/series/newSeriesNav";
import {
  resolveExperimentFilter,
  UNKNOWN_EXPERIMENT_LABEL,
} from "../../lib/experimentFilter";
import { useShortcuts } from "../shell/useShortcuts";
import { showToast } from "../../lib/toast";
import { toSampleRowModel } from "./samplesAdapters";
import {
  sortSampleRows,
  nextSortState,
  isSortableKey,
  type SortKey,
  type SortDir,
  type SortableKey,
} from "../../lib/sample/sortSamples";
import { sampleDisplayName } from "../../lib/sample/displayName";

// Boneyard fixture — a static skeleton shaped to the SheetTable rows region so
// the headless capture CLI measures the contact-sheet body. Four placeholder
// rows of token-only bars (no inline appearance literals: design-guard clean).
const CONTACT_SHEET_FIXTURE = (
  <div className="bg-paper-sunk">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="border-b border-hair px-4 py-5">
        <div className="h-3 w-40 rounded-md bg-paper-sunk" />
        <div className="mt-3 h-12 w-full rounded-md bg-paper-sunk" />
      </div>
    ))}
  </div>
);

/**
 * Single shared predicate (SA-F3): how many distinct samples a selection
 * spans, counted EXACTLY the way the Drop/Keep/Restore handler routes the
 * batch — distinct owners among the ids `ownerOf` can map (unmappable ids
 * are skipped by the handler, so they don't count here either). Both the
 * CullBar disclosure and the toast receipt go through this one function so
 * the promise and the action can never count differently.
 */
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

/**
 * SamplesPage (greenfield) — the contact-sheet table at /samples.
 *
 * Reimplemented from src/print composites + the sample-table mockup: a
 * `PageFrame width="sheet"` body with the contact-sheet head (kicker + serif
 * experiment title + screened progress), a slotted `SheetTable` of
 * `SampleTableRow`s, a footer keyboard legend, and a floating `CullBar`.
 *
 * Carried logic only (queries / mutators / the row adapter); no legacy
 * presentation. Selection is page-owned (a Set of exposure ids); a Shift-anchor
 * ref drives the contiguous-range cull within a single sample's exposure order.
 */
export function SamplesPage(): JSX.Element {
  const navigate = useNavigate();
  const [searchParams, setSearchParams] = useSearchParams();

  // ── filter + title (carried logic) ──────────────────────────────────────────
  const corpusQuery = useCorpusSamples();
  const experimentsQuery = useExperiments();
  // SA-F5: the SAME resolver the topbar select uses — the page body and the
  // select must agree on whether the URL names a real experiment.
  const filter = resolveExperimentFilter(
    searchParams.get("experiment"),
    experimentsQuery.data,
    corpusQuery.data,
  );
  const experimentId = filter.id;

  const samples = useMemo(
    () => corpusQuery.data ?? [],
    [corpusQuery.data],
  );
  const filtered = useMemo(
    () =>
      experimentId === undefined
        ? samples
        : samples.filter((s) => s.experiment_id === experimentId),
    [samples, experimentId],
  );

  // The unknown title never repeats the raw id from the address; the loading
  // fallback (`experiment N`) only shows before the lists settle the verdict.
  const title =
    experimentId === undefined
      ? "The corpus"
      : filter.unknown
        ? UNKNOWN_EXPERIMENT_LABEL
        : (filter.name ?? `experiment ${experimentId}`);

  function clearExperimentFilter(): void {
    setSearchParams((prev) => {
      const next = new URLSearchParams(prev);
      next.delete("experiment");
      return next;
    });
  }

  const corpusExposures = useCorpusExposures(filtered);
  const { screened, total } = useScreenedProgress(filtered);
  const batch = useSetExposureStatusBatch();

  // ── sort (SA-SORT, page-owned, URL-persisted) ───────────────────────────────
  // The sort key+dir live in the address the SAME way the experiment filter does
  // (useSearchParams), so a sorted view is a shareable permalink. An unknown /
  // absent ?sort falls back to null = ingest order; ?dir defaults to asc.
  const sortKeyParam = searchParams.get("sort");
  const sortKey: SortKey =
    sortKeyParam && isSortableKey(sortKeyParam) ? sortKeyParam : null;
  const sortDir: SortDir = searchParams.get("dir") === "desc" ? "desc" : "asc";
  const sort = { key: sortKey, dir: sortDir };

  function applySort(key: SortableKey): void {
    const next = nextSortState(sort, key);
    setSearchParams((prev) => {
      const params = new URLSearchParams(prev);
      if (next.key === null) {
        params.delete("sort");
        params.delete("dir");
      } else {
        params.set("sort", next.key);
        params.set("dir", next.dir);
      }
      return params;
    });
  }

  // Apply the stable sort to the FILTERED list before it becomes rows. The
  // exposures/kept accessors read the same per-sample exposure rows the row
  // adapter does, so the column you sort by is exactly the column you see.
  const sortedSamples = useMemo(() => {
    if (sortKey === null) return filtered;
    return sortSampleRows(filtered, { key: sortKey, dir: sortDir }, {
      name: (s) => sampleDisplayName(s),
      exposures: (s) => (corpusExposures.byId.get(s.id) ?? []).length,
      kept: (s) =>
        (corpusExposures.byId.get(s.id) ?? []).filter((e) => e.status === "accepted")
          .length,
      phase: (s) => s.phase,
    });
  }, [filtered, sortKey, sortDir, corpusExposures.byId]);

  // LO-NEXT: hand the loupe the exact visible order (sorted + experiment-filtered)
  // through router state, so its prev/next walk matches what's on screen. (The
  // "exposures"/"kept" sorts can't be reconstructed from the loupe's data, so
  // passing the resolved id list is the honest source.)
  const sampleOrder = useMemo(
    () => sortedSamples.map((s) => s.id),
    [sortedSamples],
  );

  // ── page-level cursor (T2.4) ────────────────────────────────────────────────
  // Tracks the "active" row (sampleIndex) and detail frame (frameIndex) driven
  // by ↑/↓/←/→ keyboard navigation via useShortcuts. Default {0,0} = first
  // sample, first frame. Both indices are clamped (never wrap around).
  const [cursor, setCursor] = useState<{ sampleIndex: number; frameIndex: number }>(
    { sampleIndex: 0, frameIndex: 0 },
  );
  const activeSample = sortedSamples[cursor.sampleIndex];

  function clamp(v: number, lo: number, hi: number): number {
    return v < lo ? lo : v > hi ? hi : v;
  }

  function clampSample(d: number): void {
    setCursor((c) => ({
      sampleIndex: clamp(c.sampleIndex + d, 0, Math.max(0, sortedSamples.length - 1)),
      frameIndex: 0,
    }));
  }

  function clampFrame(d: number): void {
    setCursor((c) => {
      const frames = corpusExposures.byId.get(sortedSamples[c.sampleIndex]?.id ?? -1) ?? [];
      return {
        ...c,
        frameIndex: clamp(c.frameIndex + d, 0, Math.max(0, frames.length - 1)),
      };
    });
  }

  // ── selection state (page-owned) ────────────────────────────────────────────
  // Exposure-grain cull selection (drives CullBar → Drop/Restore).
  const [selected, setSelected] = useState<Set<number>>(() => new Set());
  const anchorRef = useRef<{ sampleId: number; exposureId: number } | null>(null);
  const shiftRef = useRef(false);

  // Sample-grain pick selection — DISTINCT from the exposure-grain cull above.
  // Drives the ComposeBar → "+ New series" carry; never merged with `selected`.
  const [checkedSamples, setCheckedSamples] = useState<Set<number>>(() => new Set());

  function toggleSampleCheck(sampleId: number): void {
    setCheckedSamples((prev) => {
      const next = new Set(prev);
      if (next.has(sampleId)) next.delete(sampleId);
      else next.add(sampleId);
      return next;
    });
  }

  // SA-STALESELECT: an experiment-filter change re-scopes the visible corpus, so a
  // working selection from the prior scope is no longer on screen to act on or
  // uncheck. Carrying it forward made the CullBar/ComposeBar counts, the cull
  // toast receipt, and the "+ New series" carry lie about what's actionable
  // (e.g. "1 sample" with no visible row, or a drop toast counting frames whose
  // sample left the filter). Clear both grains whenever the scope changes; on
  // mount this is a no-op (both sets start empty).
  useEffect(() => {
    setSelected(new Set());
    setCheckedSamples(new Set());
    anchorRef.current = null;
  }, [experimentId]);

  // LA-COLLIDE: the bottom Dock is always mounted and permanently claims the
  // dock lane via its own useEffect (Dock.tsx). CullBar / ComposeBar are opaque
  // transient bars that float above the Dock; the lane is always occupied while
  // this page is mounted, so no additional lane claim is needed here.

  useEffect(() => {
    function onShiftDown(e: KeyboardEvent): void {
      if (e.key === "Shift") shiftRef.current = true;
    }
    function onShiftUp(e: KeyboardEvent): void {
      if (e.key === "Shift") shiftRef.current = false;
    }
    window.addEventListener("keydown", onShiftDown);
    window.addEventListener("keyup", onShiftUp);
    return () => {
      window.removeEventListener("keydown", onShiftDown);
      window.removeEventListener("keyup", onShiftUp);
    };
  }, []);

  // exposureId → owning sampleId (for routing a batch mutate back to its sample).
  const ownerOf = useMemo(() => {
    const m = new Map<number, number>();
    for (const s of filtered)
      for (const e of corpusExposures.byId.get(s.id) ?? []) m.set(e.id, s.id);
    return m;
  }, [filtered, corpusExposures.byId]);

  function toggleSelect(sampleId: number, exposureId: number): void {
    setSelected((prev) => {
      const next = new Set(prev);
      const anchor = anchorRef.current;
      if (shiftRef.current && anchor && anchor.sampleId === sampleId) {
        // Contiguous range within THIS sample's exposure order.
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
      // Single toggle + (re)anchor.
      if (next.has(exposureId)) next.delete(exposureId);
      else next.add(exposureId);
      anchorRef.current = { sampleId, exposureId };
      return next;
    });
  }

  function batchSet(status: "accepted" | "rejected" | null): void {
    // SA-STALESELECT: only act on exposures still in the visible scope (an owner
    // resolves in `ownerOf`, which is built from `filtered`). The toast then
    // counts what was actually mutated, never the raw selection size — a
    // destructive action must never claim a larger drop/keep count than it made.
    const targets = [...selected].filter((id) => ownerOf.has(id));
    for (const id of targets) {
      batch.mutate({ sampleId: ownerOf.get(id)!, exposureId: id, status });
    }
    const n = targets.length;
    // Consequential, batch-scale change → visible toast. Drop, Keep and the
    // symmetric restore all announce so each action and its undo are equally
    // legible. Restore sends null unconditionally: it clears BOTH verdicts.
    if (n > 0) {
      const verb =
        status === "rejected" ? "dropped" : status === "accepted" ? "kept" : "restored";
      // SA-F3: the receipt matches the bar's promise — a cross-sample batch
      // discloses its spread with the same suffix, under the same predicate.
      const spread = selectionSpread(new Set(targets), ownerOf);
      const suffix = spread > 1 ? ` across ${spread} samples` : "";
      showToast(`${n} frame${n === 1 ? "" : "s"} ${verb}${suffix}`, "success");
    }
    setSelected(new Set());
    anchorRef.current = null;
  }

  // ── keyboard map (T2.4): page cursor + selection + cull verbs ──────────────
  //
  // All page-level shortcuts go through the shared registry/useShortcuts hook —
  // no hand-rolled window keydown switch. Modifier chords (Alt+↑/↓ = reorder,
  // Shift+←/→ = extend) are DISTINCT combos in the registry, never confused
  // with the bare-arrow nav. Typing/modal suppression lives in useShortcuts.
  //
  // Declining protocol: a handler returns `false` to DECLINE (no preventDefault,
  // event keeps propagating for the Esc-ladder etc.); returning nothing/undefined
  // means "I handled it" (preventDefault fires).
  //
  // NO `representative` binding on Corpus — that verb is Focus-page only.
  useShortcuts({
    // ── Navigation: sample-axis ──────────────────────────────────────────────
    prevSample: () => clampSample(-1),
    nextSample: () => clampSample(1),
    // ── Navigation: detail-axis (frame within active sample) ─────────────────
    prevDetail: () => clampFrame(-1),
    nextDetail: () => clampFrame(1),
    // ── Navigation: open Focus / Loupe ───────────────────────────────────────
    openFocus: () => {
      if (activeSample == null) return false;
      navigate(`/sample/${activeSample.id}`);
      return undefined;
    },
    openLoupe: () => {
      if (activeSample == null) return false;
      navigate(`/sample/${activeSample.id}/loupe`);
      return undefined;
    },
    // ── Selection: multi-select verbs ────────────────────────────────────────
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
      // Select all frames of the active sample (or all frames in scope).
      const allFrameIds: number[] = [];
      for (const sam of sortedSamples) {
        for (const f of corpusExposures.byId.get(sam.id) ?? []) {
          allFrameIds.push(f.id);
        }
      }
      if (allFrameIds.length === 0) return false;
      setSelected(new Set(allFrameIds));
      return undefined;
    },
    // ── Screen verbs: drop/keep/restore ──────────────────────────────────────
    // When a selection is live → act on the whole selection (via batchSet).
    // When no selection → act on the active frame (the single exposure at the
    // cursor). If there is genuinely no active frame, return false so the key
    // stays inert (don't crash on an empty list).
    drop: () => {
      if (selected.size > 0) { batchSet("rejected"); return undefined; }
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      batch.mutate({ sampleId: s.id, exposureId: frame.id, status: "rejected" });
      showToast("1 frame dropped", "success");
      return undefined;
    },
    keep: () => {
      if (selected.size > 0) { batchSet("accepted"); return undefined; }
      const s = activeSample;
      const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
      const frame = frames[cursor.frameIndex];
      if (s == null || frame == null) return false;
      batch.mutate({ sampleId: s.id, exposureId: frame.id, status: "accepted" });
      showToast("1 frame kept", "success");
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
    // ── Dismiss: clear selection → Esc ladder continues if nothing to clear ──
    // Returns `false` (declines) when there is nothing to clear, so the Esc
    // ladder can continue to the next handler.
    dismiss: () => {
      if (selected.size === 0) return false;
      setSelected(new Set());
      anchorRef.current = null;
      return undefined;
    },
    // NO representative binding on Corpus.
  });

  // SA-F2: a double-clicked frame carries its exposure id into the loupe via
  // ?exposure= so the loupe opens AT that frame. A search param (not router
  // state) keeps the frame on refresh/permalink, the house slug-permalink
  // direction. Frameless openings (name click) omit it → loupe default frame.
  function loupeHref(id: number, exposureId?: number): string {
    const params = new URLSearchParams();
    if (experimentId !== undefined) params.set("experiment", String(experimentId));
    if (exposureId !== undefined) params.set("exposure", String(exposureId));
    const qs = params.toString();
    return `/sample/${id}/loupe${qs ? `?${qs}` : ""}`;
  }

  return (
    <PageFrame width="sheet" className="px-10 py-8">
      <div data-testid="samples-page">
        {/* ── Head ──────────────────────────────────────────────────────────── */}
        <div className="flex items-end justify-between gap-8 mb-5">
          <div>
            <Kicker tone="accent">Contact sheet</Kicker>
            <h1 className="text-display text-ink">{title}</h1>
            <p className="text-body text-ink-soft mt-2 max-w-[62ch]">
              Flip each frame and drop the ones with flares or artifacts. Tags
              are optional; you name the ordering variable later, when you scope
              a series.
            </p>
          </div>
          <div className="shrink-0 text-right">
            <div className="text-headline-lg text-ink leading-none">
              {screened}
              <b className="font-medium text-ink-soft"> / {total}</b>
            </div>
            <Kicker tone="soft" className="mt-0.5">samples screened</Kicker>
            <ProgressBar
              value={screened}
              total={total}
              label="samples screened"
              className="w-40 mt-1"
            />
          </div>
        </div>

        {/* SA-CULLHINT: a contextual, registry-driven teaching hint for the
            X/K cull gesture — shown only while nothing is selected. Once frames
            are selected the floating CullBar shows the real Drop/Keep buttons,
            so the hint retires (no duplicate affordance). */}
        {selected.size === 0 && sortedSamples.length > 0 && (
          <div className="mb-4 flex items-center gap-2">
            <HintText>Select frames, then</HintText>
            <KbdLegend ids={["drop", "keep"]} testId="samples-cull-hint" />
          </div>
        )}

        {/* ── Table ─────────────────────────────────────────────────────────── */}
        {corpusQuery.isError ? (
          <div data-testid="samples-error">
            <EmptyState
              title="Could not load the corpus"
              body="The sample list request failed."
              action={
                <Button
                  variant="outline"
                  disabled={corpusQuery.isFetching}
                  onClick={() => void corpusQuery.refetch()}
                >
                  Reload the corpus
                </Button>
              }
            />
          </div>
        ) : filter.unknown ? (
          // SA-F5: an address filtering by an experiment that names nothing gets
          // an honest empty state, never the bare "No samples in this experiment"
          // placeholder (which would imply the experiment exists). The single
          // action embodies the way forward: clear the filter.
          <div data-testid="samples-unknown-experiment">
            <EmptyState
              title={UNKNOWN_EXPERIMENT_LABEL}
              body="The address filters by an experiment that is not in this corpus."
              action={
                <Button variant="outline" onClick={clearExperimentFilter}>
                  Show all experiments
                </Button>
              }
            />
          </div>
        ) : (
          <Skeleton
            name="contact-sheet"
            className="block"
            loading={corpusQuery.isLoading}
            stagger={50}
            transition={200}
            fixture={CONTACT_SHEET_FIXTURE}
            fallback={
              <div className="p-8 text-sm text-ink-soft">Loading samples…</div>
            }
          >
            <SheetTable
              checkboxColumn
              sort={sort}
              onSort={applySort}
              empty={
                // ON-EMPTY: a real door, consistent with the error /
                // unknown-experiment states. The corpus is populated by ingesting
                // and analyzing an experiment (CLI — there is no in-app upload),
                // so the body names that path and the action reloads once it has
                // run, rather than dead-ending on a lone sentence.
                <EmptyState
                  title={
                    experimentId === undefined
                      ? "No samples yet"
                      : "No samples in this experiment yet"
                  }
                  body="Samples appear after an experiment is ingested and analyzed. If you just added data, reload to pull it in."
                  action={
                    <Button
                      variant="outline"
                      disabled={corpusQuery.isFetching}
                      onClick={() => void corpusQuery.refetch()}
                    >
                      Reload the corpus
                    </Button>
                  }
                />
              }
            >
              {sortedSamples.map((s, rowIndex) => {
                const loadedExposures = corpusExposures.byId.get(s.id);
                const m = toSampleRowModel(s, loadedExposures);
                // SA-ZEROEXP: distinguish a CONFIRMED-empty sample (query resolved
                // with []) from one still loading (undefined). Only the former
                // gets the terminal "No exposures" status + suppressed door.
                const noExposures = loadedExposures !== undefined && loadedExposures.length === 0;
                // T2.4: pointer click on the row sets the page cursor so keyboard
                // nav and pointer nav share a single active-item concept.
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
                    checked={checkedSamples.has(s.id)}
                    onCheck={() => toggleSampleCheck(s.id)}
                    selectedExposureIds={selected}
                    onSelectExposure={(id) => {
                      setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                      toggleSelect(s.id, id);
                    }}
                    onActivateExposure={(id) => {
                      setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                      navigate(loupeHref(s.id, id), { state: { sampleOrder } });
                    }}
                    onOpenLoupe={() => {
                      setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                      navigate(loupeHref(s.id), { state: { sampleOrder } });
                    }}
                    onOpenFocus={() => {
                      setCursor((c) => ({ ...c, sampleIndex: rowIndex }));
                      navigate(`/sample/${s.id}`);
                    }}
                    {...(hasDrop ? { onRestore: () => {
                      // T2.4: restore verb applied to this specific sample's
                      // dropped exposures. Mirrors the CullBar restore (null
                      // status batchSet) but scoped to one row's drops.
                      const drops = (corpusExposures.byId.get(s.id) ?? [])
                        .filter((e) => e.status === "rejected")
                        .map((e) => e.id);
                      if (drops.length === 0) return;
                      for (const id of drops) {
                        batch.mutate({ sampleId: s.id, exposureId: id, status: null });
                      }
                      showToast(
                        `${drops.length} frame${drops.length === 1 ? "" : "s"} restored`,
                        "success",
                      );
                    }} : {})}
                  />
                );
              })}
            </SheetTable>
          </Skeleton>
        )}

        {/* ── Footer keyboard legend ────────────────────────────────────────── */}
        {/* SA-RESCORE3 F11: the legend documents cell navigation / selection —
            meaningless under an empty, errored, or unknown-experiment table (all
            of which leave sortedSamples empty). Render it only when there are
            rows to drive. */}
        {sortedSamples.length > 0 && (
        <KbLegend
          className="mt-4"
          shortcuts={[
            { keyLabel: "↑ ↓ ← →", description: "move between cells" },
            { keyLabel: "Enter", description: "open, sort a header, or enter a cell" },
            { keyLabel: "click", description: "select a frame" },
            { keyLabel: "⇧ click", description: "extend the range" },
            // SA-KBDDUP: X/K (drop/keep) live in the contextual cull hint above
            // the table (shown while nothing is selected) and on the CullBar once
            // frames are selected; not repeated here to keep the legend shorter.
            { keyLabel: "double-click", description: "open the loupe" },
            { keyLabel: "⇧ Enter", description: "open the loupe at a frame" },
            { keyLabel: "Esc", description: "exit a cell or clear" },
            { keyLabel: "⌘K", description: "find a sample" },
          ]}
        />
        )}
      </div>

      {/* ── Floating cull bar (page root) ─────────────────────────────────────── */}
      <CullBar
        count={selected.size}
        sampleCount={selectionSpread(selected, ownerOf)}
        show={selected.size > 0}
        onReject={() => batchSet("rejected")}
        onKeep={() => batchSet("accepted")}
        onRestore={() => batchSet(null)}
        onClear={() => {
          setSelected(new Set());
          anchorRef.current = null;
        }}
      />

      {/* ── Floating sample-grain compose bar (sits above CullBar) ─────────────── */}
      <ComposeBar
        count={checkedSamples.size}
        show={checkedSamples.size > 0}
        onNewSeries={() => navigateToNewSeries(checkedSamples, navigate)}
        onClear={() => setCheckedSamples(new Set())}
      />

      {/* ── Contextual bottom dock (Corpus grammar §3.3) ─────────────────────────
          ‹ Experiments · Sample↑↓ · Frame‹› · Drop · Keep · Restore · Loupe · Focus
          NO Set-representative on Corpus (keyboard comment above confirms this).
          Each verb calls the SAME callback the keyboard shortcut uses — no divergence. */}
      <Dock>
        {/* Up-link */}
        <a
          href="/experiments"
          onClick={(e) => { e.preventDefault(); navigate("/experiments"); }}
          className="text-meta font-semibold text-print-accent hover:underline mr-1"
          data-testid="dock-up-link"
        >
          ‹ Experiments
        </a>

        <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

        {/* Sample↑↓ stepper */}
        <IconButton
          label="Previous sample"
          tone="ghost"
          disabled={cursor.sampleIndex === 0}
          onClick={() => clampSample(-1)}
          data-testid="dock-prev-sample"
        >
          ↑
        </IconButton>
        <IconButton
          label="Next sample"
          tone="ghost"
          disabled={cursor.sampleIndex >= sortedSamples.length - 1}
          onClick={() => clampSample(1)}
          data-testid="dock-next-sample"
        >
          ↓
        </IconButton>

        <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

        {/* Frame‹› stepper */}
        <IconButton
          label="Previous frame"
          tone="ghost"
          disabled={cursor.frameIndex === 0}
          onClick={() => clampFrame(-1)}
          data-testid="dock-prev-frame"
        >
          ‹
        </IconButton>
        <IconButton
          label="Next frame"
          tone="ghost"
          disabled={(() => {
            const frames = corpusExposures.byId.get(sortedSamples[cursor.sampleIndex]?.id ?? -1) ?? [];
            return cursor.frameIndex >= frames.length - 1;
          })()}
          onClick={() => clampFrame(1)}
          data-testid="dock-next-frame"
        >
          ›
        </IconButton>

        <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

        {/* Cull verbs: Drop / Keep / Restore — coloured outlines (spec §3.3) */}
        <Button
          variant="outlineDanger"
          onClick={() => {
            const s = activeSample;
            const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
            const frame = frames[cursor.frameIndex];
            if (s == null || frame == null) return;
            batch.mutate({ sampleId: s.id, exposureId: frame.id, status: "rejected" });
          }}
          data-testid="dock-drop"
        >
          Drop
        </Button>
        <Button
          variant="outlineSuccess"
          onClick={() => {
            const s = activeSample;
            const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
            const frame = frames[cursor.frameIndex];
            if (s == null || frame == null) return;
            batch.mutate({ sampleId: s.id, exposureId: frame.id, status: "accepted" });
          }}
          data-testid="dock-keep"
        >
          Keep
        </Button>
        <Button
          variant="ghost"
          onClick={() => {
            const s = activeSample;
            const frames = s != null ? (corpusExposures.byId.get(s.id) ?? []) : [];
            const frame = frames[cursor.frameIndex];
            if (s == null || frame == null) return;
            batch.mutate({ sampleId: s.id, exposureId: frame.id, status: null });
          }}
          data-testid="dock-restore"
        >
          Restore
        </Button>

        <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

        {/* Destination buttons */}
        <Button
          variant="ghost"
          onClick={() => {
            if (activeSample == null) return;
            navigate(`/sample/${activeSample.id}/loupe`);
          }}
          data-testid="dock-loupe"
        >
          Loupe
        </Button>
        <Button
          variant="accent"
          onClick={() => {
            if (activeSample == null) return;
            navigate(`/sample/${activeSample.id}`);
          }}
          data-testid="dock-focus"
        >
          Focus
        </Button>
      </Dock>
    </PageFrame>
  );
}
