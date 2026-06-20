import { useNavigate } from "react-router-dom";
import { useExperiments } from "../../queries";
import { useAppState } from "../../state";
import { PageFrame } from "../components/PageFrame";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import { Dot } from "../ui/Dot";
import { EmptyState } from "../ui/EmptyState";
import type { DotTone } from "../ui/Dot";
import type { Experiment } from "../../api";

// Dot tones are accent | success | muted | neutral (Dot.tsx — there is NO
// `ok`/`danger` tone). failed → muted (a quiet hollow ring, NOT a red alarm —
// the row's "Scan failed" text carries the severity), scanning/analyzing →
// accent, complete → success.
function statusTone(s: Experiment["ingest_status"]): DotTone {
  if (s === "failed") return "muted";
  if (s === "scanning" || s === "analyzing") return "accent";
  if (s === "complete") return "success";
  return "neutral";
}

/**
 * ExperimentsHomePage — /experiments gallery (spec §7/§8.7). Cards carry a
 * status dot, serif name, `date range · directory` meta, counts, and a
 * "N need grouping review" hint; selecting one persists activeExperimentId and
 * routes to its corpus. Empty → EmptyState + New CTA.
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

  return (
    <PageFrame width="home" className="px-6 py-8">
      <div className="flex items-start justify-between gap-6 mb-6">
        <div>
          <Kicker>Experiments</Kicker>
          <h1 className="text-display text-ink">Your beamtimes</h1>
        </div>
        <Button
          variant="accent"
          data-testid="new-experiment-cta"
          onClick={() => navigate("/experiments/new")}
        >
          + New experiment
        </Button>
      </div>

      {experiments.isSuccess && list.length === 0 ? (
        <EmptyState
          title="No experiments yet"
          body="Point Himalaya at a directory of exposures and it scans, groups them into samples, and derives the geometry. No manifest needed."
          action={
            <Button variant="accent" onClick={() => navigate("/experiments/new")}>
              + New experiment
            </Button>
          }
        />
      ) : (
        <ul className="flex flex-col gap-3">
          {list.map((e) => (
            <Card
              as="li"
              key={e.id}
              interactive
              padding="md"
              data-testid={`experiment-card-${e.id}`}
              onClick={() => open(e.id)}
            >
              <div className="flex items-center gap-3">
                <Dot tone={statusTone(e.ingest_status)} />
                <span className="text-headline text-ink">{e.name ?? `Experiment ${e.id}`}</span>
                <span className="ml-auto font-mono text-sm text-ink-soft">{e.data_dir}</span>
              </div>
            </Card>
          ))}
        </ul>
      )}
    </PageFrame>
  );
}
