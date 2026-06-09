import { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { SheetTable } from "../components/SheetTable";
import { SampleTableRow } from "../components/SampleTableRow";
import { CullBar } from "../components/CullBar";
import { Kicker, KbLegend, ProgressBar, ComposeBar } from "../ui";
import {
  useCorpusSamples,
  useCorpusExposures,
  useScreenedProgress,
  useExperiments,
  useSetExposureStatusBatch,
} from "../../queries";
import { navigateToNewSeries } from "../../lib/series/newSeriesNav";
import { toSampleRowModel } from "./samplesAdapters";

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
 * SamplesPage (greenfield) — the contact-sheet table at /samples.
 *
 * Reimplemented from src/print composites + the sample-table mockup: a
 * `PageFrame width="sheet"` body with the contact-sheet head (kicker + serif
 * beamtime title + screened progress), a slotted `SheetTable` of
 * `SampleTableRow`s, a footer keyboard legend, and a floating `CullBar`.
 *
 * Carried logic only (queries / mutators / the row adapter); no legacy
 * presentation. Selection is page-owned (a Set of exposure ids); a Shift-anchor
 * ref drives the contiguous-range cull within a single sample's exposure order.
 */
export function SamplesPage(): JSX.Element {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();

  // ── filter + title (carried logic) ──────────────────────────────────────────
  const raw = searchParams.get("beamtime");
  const beamtime = raw !== null && /^\d+$/.test(raw) ? Number(raw) : undefined;
  const beamtimeQuery = beamtime !== undefined ? `?beamtime=${beamtime}` : "";

  const corpusQuery = useCorpusSamples();
  const experimentsQuery = useExperiments();
  const samples = useMemo(
    () => corpusQuery.data ?? [],
    [corpusQuery.data],
  );
  const filtered = useMemo(
    () =>
      beamtime === undefined
        ? samples
        : samples.filter((s) => s.experiment_id === beamtime),
    [samples, beamtime],
  );

  const title =
    beamtime === undefined
      ? "The corpus"
      : (experimentsQuery.data?.find((e) => e.id === beamtime)?.name ??
        `experiment ${beamtime}`);

  const corpusExposures = useCorpusExposures(filtered);
  const { screened, total } = useScreenedProgress(filtered);
  const batch = useSetExposureStatusBatch();

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

  function batchSet(status: "rejected" | null): void {
    for (const id of selected) {
      const sampleId = ownerOf.get(id);
      if (sampleId !== undefined) batch.mutate({ sampleId, exposureId: id, status });
    }
    setSelected(new Set());
    anchorRef.current = null;
  }

  // ── keyboard: X drops, Esc clears (only meaningful with a selection) ─────────
  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const tag = (e.target as HTMLElement | null)?.tagName;
      if (tag === "INPUT" || tag === "TEXTAREA") return;
      if (selected.size === 0) return;
      if (e.key === "x" || e.key === "X") batchSet("rejected");
      else if (e.key === "Escape") {
        setSelected(new Set());
        anchorRef.current = null;
      }
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selected, ownerOf]);

  function loupeHref(id: number): string {
    return `/samples/loupe/${id}${beamtimeQuery}`;
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
              Flip the frames and drop the ones with flares or artifacts. Tags
              are a light, optional note on what each sample is — the ordering
              variable is named later, when you scope a series.
            </p>
          </div>
          <div className="shrink-0 text-right">
            <div className="text-headline-lg text-ink leading-none">
              {screened}
              <b className="font-medium text-ink-faint"> / {total}</b>
            </div>
            <Kicker tone="faint" className="mt-0.5">samples screened</Kicker>
            <ProgressBar
              value={screened}
              total={total}
              label="samples screened"
              className="w-40 mt-1"
            />
          </div>
        </div>

        {/* ── Table ─────────────────────────────────────────────────────────── */}
        {corpusQuery.isError ? (
          <div
            data-testid="samples-error"
            className="rounded-md border border-hair-strong p-8 text-center text-ink-faint"
          >
            Could not load the corpus. Try reloading.
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
              <div className="p-8 text-sm text-ink-faint">Loading samples…</div>
            }
          >
            <SheetTable
              checkboxColumn
              empty={
                <div className="p-10 text-center text-ink-faint">
                  No samples{" "}
                  {beamtime === undefined ? "in the corpus" : "in this beamtime"}{" "}
                  yet.
                </div>
              }
            >
              {filtered.map((s) => {
                const m = toSampleRowModel(s, corpusExposures.byId.get(s.id));
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
                    tags={m.tags}
                    {...(m.phase !== undefined ? { phase: m.phase } : {})}
                    checked={checkedSamples.has(s.id)}
                    onCheck={() => toggleSampleCheck(s.id)}
                    selectedExposureIds={selected}
                    onSelectExposure={(id) => toggleSelect(s.id, id)}
                    onActivateExposure={() => navigate(loupeHref(s.id))}
                    onOpenLoupe={() => navigate(loupeHref(s.id))}
                    onOpenFocus={() => navigate(`/sample/${s.id}`)}
                  />
                );
              })}
            </SheetTable>
          </Skeleton>
        )}

        {/* ── Footer keyboard legend ────────────────────────────────────────── */}
        <KbLegend
          className="mt-4"
          shortcuts={[
            { keyLabel: "click", description: "select a frame" },
            { keyLabel: "⇧ click", description: "extend the range" },
            { keyLabel: "X", description: "drop the selected frames" },
            { keyLabel: "double-click", description: "open the loupe" },
            { keyLabel: "Esc", description: "clear" },
          ]}
        />
      </div>

      {/* ── Floating cull bar (page root) ─────────────────────────────────────── */}
      <CullBar
        count={selected.size}
        show={selected.size > 0}
        onReject={() => batchSet("rejected")}
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
    </PageFrame>
  );
}
