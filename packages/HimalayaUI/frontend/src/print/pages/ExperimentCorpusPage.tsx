import { Link, useParams } from "react-router-dom";
import { useAppState } from "../../state";
import { useLoads } from "../../queries";
import { Badge } from "../ui/Badge";

/**
 * ExperimentCorpusPage — the Corpus tab body. Reuses the shipped SheetTable
 * contact sheet scoped to the experiment (E2 wires the table), with a sticky
 * grouping-review banner above it linking to E2's GroupingReviewPage
 * (/experiments/:id/grouping). The live-ingest unfold (E2 LiveIngestUnfold)
 * replaces the table while ingestInFlight[id] is active.
 */
export function ExperimentCorpusPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const inFlight = useAppState((s) => s.ingestInFlight?.[expId]);
  const loads = useLoads(expId);

  // Review count: LoadSamples across all loads whose `flag` is non-null
  // (a flagged merge/split discrepancy). Derived from useLoads — tests mock
  // api.listLoads via vi.spyOn, same pattern as Tasks 11–14.
  const reviewCount = (loads.data ?? []).reduce(
    (n, l) => n + l.samples.filter((s) => s.flag !== null).length,
    0,
  );
  const processing = inFlight?.status === "scanning" || inFlight?.status === "analyzing";

  return (
    <div className="flex flex-col gap-4">
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

      {processing ? (
        // E2 LiveIngestUnfold slots here (ProgressBar + StatBar counters +
        // skeleton rows driven by ingestInFlight). E1 renders the labelled slot.
        <div data-testid="live-ingest-slot" className="text-sm text-ink-soft">
          Processing exposures…
        </div>
      ) : (
        // E2 mounts the scoped SheetTable here. E1 renders the labelled slot so
        // the page is assemblable without the corpus query wiring.
        <div data-testid="corpus-sheet-slot" className="text-sm text-ink-soft">
          {/* SheetTable scoped to experiment {expId} — wired in E2 */}
        </div>
      )}
    </div>
  );
}
