import { useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useSeries } from "../queries";

/** Static skeleton for boneyard's headless capture: the plate area + a rail. */
const BUILDER_FIXTURE = (
  <div className="grid grid-cols-[1fr_336px] gap-0">
    <div className="m-4 rounded border border-hair bg-plate" style={{ aspectRatio: "10 / 3" }} />
    <div className="border-l border-hair p-4">
      <div className="h-4 w-1/2 rounded bg-paper-sunk" />
    </div>
  </div>
);

/**
 * SeriesBuilderPage — the series builder visual surface at /series/:id
 * (#175 / I3.5a). Read-only: reads one series via useSeries(id) and composes
 * the existing MultiTracePlot render core (see Task 3+). Mutations (recipe
 * edits, plate commit, permalink) are I3.5b — NOT here. Mounted under the
 * CorpusShell layout route, the destination the I3.3 folio card links to.
 *
 * URL-owned: the series id comes from the route param, never Zustand.
 */
export function SeriesBuilderPage(): JSX.Element {
  const { id: idParam } = useParams<{ id: string }>();
  const seriesId = Number(idParam);
  const query = useSeries(Number.isFinite(seriesId) ? seriesId : undefined);

  if (query.isError) {
    return (
      <div data-testid="series-builder-error" className="p-8 text-sm text-ink-soft">
        Could not load this series. It may have been deleted.
      </div>
    );
  }

  const s = query.data;
  const title = s && s.title.trim() !== "" ? s.title : "Untitled series";

  return (
    <div data-testid="series-builder-page" className="flex h-full min-h-0 flex-col">
      <Skeleton
        name="series-builder"
        className="flex-1 min-h-0 flex flex-col"
        loading={query.isLoading}
        fixture={BUILDER_FIXTURE}
      >
        {s && (
          <header data-testid="series-builder-header" className="shrink-0 px-6 pt-5">
            <div className="text-xs font-semibold uppercase tracking-wide text-print-accent">
              Series
            </div>
            <h1 className="font-medium text-ink">{title}</h1>
          </header>
        )}
        {/* Plot + rail composed in Task 3+. */}
      </Skeleton>
    </div>
  );
}
