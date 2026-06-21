import { useMemo } from "react";
import { Link, useNavigate, useParams } from "react-router-dom";
import { useAppState } from "../../state";
import { useExperiment, useLoads } from "../../queries";
import { Badge } from "../ui/Badge";
import { ProgressBar } from "../ui/ProgressBar";
import { ScanFailedPage } from "./ScanFailedPage";
import { GroupingReviewPage } from "../components/GroupingReviewPage";
import { SheetTable } from "../components/SheetTable";
import { SampleTableRow } from "../components/SampleTableRow";

/**
 * ExperimentCorpusPage — the Corpus tab body.
 *
 * Implements the §6.2 state machine:
 *   - scanning   (inFlight.status==="scanning")      → GroupingReviewPage (live-unfold)
 *   - rescanning (inFlight.status==="analyzing")     → inline ProgressBar
 *   - failed     (!processing && status==="failed")  → ScanFailedPage  [T4.2, preserved]
 *   - has-flags  (reviewCount > 0)                   → banner + SheetTable
 *   - clean      (reviewCount === 0)                 → SheetTable only
 *
 * The SheetTable is scoped to this experiment's samples by flattening
 * useLoads(expId) → LoadSample[]. No per-sample exposure fan-out is needed for
 * the corpus view: exposures inside each LoadSample provide the frame count.
 */
export function ExperimentCorpusPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const navigate = useNavigate();

  const inFlight = useAppState((s) => s.ingestInFlight?.[expId]);
  const loads = useLoads(expId);
  const exp = useExperiment(expId);

  // --- State machine derivation ---
  // scanning: initial scan in progress (live-unfold surface)
  const scanning = inFlight?.status === "scanning";
  // rescanning: re-analysis pass (progress bar, table data not ready yet)
  const rescanning = !scanning && inFlight?.status === "analyzing";
  // processing: either active pass (scanning or rescanning)
  const processing = scanning || rescanning;
  // failed: terminal failure from the server (or ephemeral inFlight failure)
  const failed =
    !processing &&
    (inFlight?.status === "failed" || exp.data?.ingest_status === "failed");

  // Flatten Load▸Sample tree for review-count + table rows.
  const loadSamples = useMemo(
    () => (loads.data ?? []).flatMap((l) => l.samples),
    [loads.data],
  );

  // Review count: LoadSamples whose `flag` is non-null (merge/split discrepancy).
  const reviewCount = useMemo(
    () => loadSamples.filter((s) => s.flag !== null).length,
    [loadSamples],
  );

  // scanning → live GroupingReviewPage (the live-unfold surface for the initial scan)
  if (scanning) {
    return (
      <GroupingReviewPage
        experimentId={expId}
        onBack={() => navigate(`/experiments/${expId}/corpus`)}
      />
    );
  }

  // rescanning → inline progress driven by ingestInFlight
  if (rescanning) {
    return (
      <div data-testid="live-ingest-slot" className="flex flex-col gap-3">
        <p className="text-sm text-ink-soft">Analyzing exposures…</p>
        <ProgressBar
          value={inFlight ? inFlight.processed : 0}
          total={inFlight ? Math.max(inFlight.total, 1) : 1}
          label="Analysis progress"
        />
      </div>
    );
  }

  // failed → ScanFailedPage (T4.2; preserved exactly)
  if (failed) {
    return (
      <ScanFailedPage
        experimentId={expId}
        unmatched={[]}
        parsedCount={0}
      />
    );
  }

  // clean / has-flags → SheetTable scoped to this experiment's samples
  return (
    <div className="flex flex-col gap-4">
      {/* grouping-review banner (has-flags branch) */}
      {reviewCount > 0 && (
        <div
          data-testid="grouping-review-banner"
          className="sticky top-0 z-10 flex items-center gap-3 rounded-sm border border-hair-strong bg-paper-sunk px-4 py-2.5"
        >
          <span className="text-sm text-ink">
            {reviewCount} {reviewCount === 1 ? "sample needs" : "samples need"} grouping review
            <Badge>{reviewCount}</Badge>
          </span>
          <Link
            to={`/experiments/${expId}/grouping`}
            data-testid="grouping-review-link"
            className="ml-auto text-sm font-semibold text-accent hover:underline"
          >
            Review grouping →
          </Link>
        </div>
      )}

      {/* Corpus sheet — one SampleTableRow per LoadSample, scoped to expId.
          Exposures are taken from the nested LoadSample.exposures (the roll-up
          already fetched by useLoads), so no per-sample fan-out is needed here. */}
      <SheetTable>
        {loadSamples.map((s) => (
          <SampleTableRow
            key={s.sample_id}
            name={s.name}
            sampleId={`#${s.sample_id}`}
            exposures={[]}
            kept={0}
            total={s.exposures.length}
            tags={[]}
            onOpenFocus={() => navigate(`/sample/${s.sample_id}`)}
            onOpenLoupe={() => navigate(`/sample/${s.sample_id}/loupe`)}
          />
        ))}
      </SheetTable>
    </div>
  );
}
