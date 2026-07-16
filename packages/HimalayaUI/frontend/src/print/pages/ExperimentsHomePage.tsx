import { useNavigate } from "react-router-dom";
import { useExperiments } from "../../queries";
import { useAppState } from "../../state";
import { PageFrame } from "../components/PageFrame";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import { Dot } from "../ui/Dot";
import { EmptyState } from "../ui/EmptyState";
import { NoticePill } from "../ui/NoticePill";
import type { NoticePillTone } from "../ui/NoticePill";
import type { Experiment } from "../../api";

// ─── State chip ─────────────────────────────────────────────────────────────

/** Derive the status chip tone + label for a gallery card. Priority:
 *  1. scanning/analyzing → "scanning…" (accent tone)
 *  2. review_count > 0   → "N to review" (warning/amber tone)
 *  3. otherwise          → "up to date" (success/sage tone) */
function stateChip(e: Experiment): { tone: NoticePillTone; label: string } {
  if (e.ingest_status === "scanning" || e.ingest_status === "analyzing") {
    return { tone: "scanning", label: "scanning…" };
  }
  const rc = e.review_count ?? 0;
  if (rc > 0) return { tone: "warning", label: `${rc} to review` };
  return { tone: "success", label: "up to date" };
}

// ─── Year grouping ───────────────────────────────────────────────────────────

/** Extract the year from an ISO timestamp string, or null on failure. */
function yearOf(ts: string | null | undefined): number | null {
  if (!ts) return null;
  const y = parseInt(ts.slice(0, 4), 10);
  return isNaN(y) ? null : y;
}

/** Group experiments by year (from stats.started_at → fallback created_at),
 *  sorted DESC by year; within a year sorted DESC by started_at. */
function groupByYear(
  experiments: Experiment[],
): { year: number; items: Experiment[] }[] {
  const yearMap = new Map<number, Experiment[]>();

  for (const e of experiments) {
    const ts = e.stats?.started_at ?? e.created_at;
    const year = yearOf(ts) ?? new Date().getFullYear();
    const bucket = yearMap.get(year) ?? [];
    bucket.push(e);
    yearMap.set(year, bucket);
  }

  // Sort within each bucket DESC by started_at
  for (const bucket of yearMap.values()) {
    bucket.sort((a, b) => {
      const ta = a.stats?.started_at ?? a.created_at ?? "";
      const tb = b.stats?.started_at ?? b.created_at ?? "";
      return tb.localeCompare(ta);
    });
  }

  // Sort years DESC
  return Array.from(yearMap.entries())
    .sort(([ya], [yb]) => yb - ya)
    .map(([year, items]) => ({ year, items }));
}

// ─── Page ────────────────────────────────────────────────────────────────────

/**
 * ExperimentsHomePage — /experiments gallery (M2). Year-grouped 2-column card
 * grid with a slim left timeline rail. Each card shows the experiment name,
 * state chip (scanning / N to review / up to date), data_dir path, and counts.
 * Empty → EmptyState + New CTA.
 */
export function ExperimentsHomePage(): JSX.Element {
  const navigate = useNavigate();
  const setActiveExperiment = useAppState((s) => s.setActiveExperiment);
  const experiments = useExperiments();

  const open = (id: number): void => {
    setActiveExperiment(id);
    navigate(`/experiments/${id}/corpus`);
  };

  const list = experiments.data ?? [];
  const groups = groupByYear(list);
  // Empty gallery: the prominent banner CTA carries the call to action, so the
  // header CTA is redundant — disable it and let the banner drive (note 1).
  const isEmpty = experiments.isSuccess && list.length === 0;

  return (
    <PageFrame width="home" className="px-6 py-8">
      {/* Page header */}
      <div className="flex items-start justify-between gap-6 mb-8">
        <div>
          {/* +4px below the accent kicker so it doesn't crowd the title (note 2). */}
          <Kicker tone="accent" className="mb-1">Experiments</Kicker>
          <h1 className="text-display text-ink">All experiments</h1>
        </div>
        <Button
          variant="accent"
          data-testid="new-experiment-cta"
          disabled={isEmpty}
          onClick={() => navigate("/experiments/new")}
        >
          + New experiment
        </Button>
      </div>

      {experiments.isSuccess && list.length === 0 ? (
        <EmptyState
          title="No experiments yet"
          body={
            // Narrower measure: the body wrapped tight so the line length suits
            // the font size instead of stretching the full page width (note 3).
            <span className="block max-w-[34ch] mx-auto">
              Point Himalaya to your experiment directory. It will discover your
              setup, scan your exposures, group them into samples. Easy.
            </span>
          }
          action={
            <Button size="lg" variant="accent" onClick={() => navigate("/experiments/new")}>
              + New experiment
            </Button>
          }
        />
      ) : (
        /* Two-column layout: slim timeline rail + year-grouped card content */
        <div className="flex gap-8 items-start">
          {/* Left timeline rail */}
          <aside className="shrink-0 w-36 pt-1">
            <Kicker tone="soft" className="mb-3">Timeline</Kicker>
            <div className="flex flex-col gap-4">
              {groups.map(({ year, items }) => (
                <div key={year} className="flex items-start gap-2">
                  <Dot tone="neutral" className="mt-1 shrink-0" />
                  <div>
                    <div className="text-body text-ink font-medium">{year}</div>
                    <div className="text-xs text-ink-soft">
                      {items.length} {items.length === 1 ? "session" : "sessions"}
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </aside>

          {/* Year-grouped card sections */}
          <div className="flex-1 min-w-0 flex flex-col gap-8">
            {groups.map(({ year, items }) => (
              <section key={year}>
                {/* Year heading + hairline */}
                <div className="flex items-center gap-4 mb-4">
                  <h2 className="text-headline text-ink shrink-0">{year}</h2>
                  <div className="flex-1 border-t border-hair" />
                </div>

                {/* 2-column card grid */}
                <div className="grid grid-cols-2 gap-4">
                  {items.map((e) => {
                    const { tone, label } = stateChip(e);
                    const stats = e.stats;
                    const isScanning =
                      e.ingest_status === "scanning" ||
                      e.ingest_status === "analyzing";

                    return (
                      <Card
                        key={e.id}
                        as="button"
                        interactive
                        padding="md"
                        className="w-full text-left"
                        data-testid={`experiment-card-${e.id}`}
                        onClick={() => open(e.id)}
                      >
                        {/* Card header: name + state chip */}
                        <div className="flex items-start justify-between gap-3 mb-2">
                          <span className="text-body text-ink font-medium leading-snug">
                            {e.name ?? `Experiment ${e.id}`}
                          </span>
                          <NoticePill
                            tone={tone}
                            className="shrink-0 mt-0.5"
                          >
                            {label}
                          </NoticePill>
                        </div>

                        {/* data_dir path */}
                        <div className="font-mono text-xs text-ink-soft mb-2 truncate">
                          {e.data_dir}
                        </div>

                        {/* Counts or indexing placeholder */}
                        {isScanning || !stats ? (
                          <div className="text-xs text-ink-soft">indexing…</div>
                        ) : (
                          <div className="text-xs text-ink-soft">
                            <span className="font-mono text-ink">{stats.samples}</span>
                            {" samples · "}
                            <span className="font-mono text-ink">{stats.exposures}</span>
                            {" exp · "}
                            <span className="font-mono text-ink">{stats.loads}</span>
                            {" loads"}
                          </div>
                        )}
                      </Card>
                    );
                  })}
                </div>
              </section>
            ))}
          </div>
        </div>
      )}
    </PageFrame>
  );
}
