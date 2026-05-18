import { useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useCorpusSamples, useExperiments } from "../queries";
import { HintText } from "../components/ui/HintText";
import {
  ContactSheetRow,
  CONTACT_SHEET_COLS,
} from "../components/ContactSheetRow";

/** Static skeleton shape for the boneyard headless capture (boneyard.md). */
const CONTACT_SHEET_FIXTURE = (
  <div>
    {[0, 1, 2, 3].map((i) => (
      <div
        key={i}
        className={`${CONTACT_SHEET_COLS} border-b border-hair px-4 py-3`}
      >
        <div className="h-8 rounded bg-paper-sunk" />
        <div className="h-16 rounded bg-paper-sunk" />
        <div className="h-8 rounded bg-paper-sunk" />
        <div className="h-8 rounded bg-paper-sunk" />
        <div className="h-8 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

/**
 * SamplesPage — the corpus contact-sheet table at /samples (#160 / I1.4).
 *
 * Owns the corpus query, the ?beamtime= URL filter, and the boneyard
 * skeleton. ?beamtime=<experiment_id> is read here (the topbar chip is the
 * writer — CorpusTopbar); the corpus is filtered client-side because
 * useCorpusSamples() fetches the whole corpus and exposes no filter.
 */
export function SamplesPage(): JSX.Element {
  const [searchParams] = useSearchParams();
  const beamtimeRaw = searchParams.get("beamtime");
  const beamtime =
    beamtimeRaw !== null && /^\d+$/.test(beamtimeRaw)
      ? Number(beamtimeRaw)
      : undefined;

  const corpusQuery = useCorpusSamples();
  const experimentsQuery = useExperiments();

  const samples = corpusQuery.data ?? [];
  const filtered =
    beamtime === undefined
      ? samples
      : samples.filter((s) => s.experiment_id === beamtime);

  const scopeName =
    beamtime === undefined
      ? "the whole corpus"
      : (experimentsQuery.data?.find((e) => e.id === beamtime)?.name ??
        `experiment ${beamtime}`);

  return (
    <div data-testid="samples-page" className="flex flex-col gap-4 p-6">
      <header className="flex flex-col gap-1">
        <div className="text-xs font-semibold uppercase tracking-wide text-print-accent">
          Contact sheet
        </div>
        <p data-testid="samples-scope" className="text-sm text-ink-faint">
          Showing {scopeName}.
        </p>
      </header>

      <div
        className={`${CONTACT_SHEET_COLS} border-b border-hair-strong px-4 pb-2
                    text-xs font-semibold uppercase tracking-wide text-ink-faint`}
      >
        <div>Sample</div>
        <div>Exposures</div>
        <div>Kept</div>
        <div>Tags</div>
        <div>Status</div>
      </div>

      {corpusQuery.isError ? (
        <div data-testid="samples-error" className="px-4 py-8 text-sm text-ink-soft">
          Could not load samples. Try reloading the page.
        </div>
      ) : (
        <Skeleton
          name="contact-sheet"
          className="w-full"
          loading={corpusQuery.isLoading}
          stagger={50}
          transition={200}
          fixture={CONTACT_SHEET_FIXTURE}
          fallback={<HintText>Loading samples…</HintText>}
        >
          <div data-testid="contact-sheet-rows">
            {filtered.map((s) => (
              <ContactSheetRow key={s.id} sample={s} />
            ))}
          </div>
        </Skeleton>
      )}

      {/* Provisional preview copy from the sample-table.html mockup. These
          keybindings are wired by culling (#162) and the loupe (#161), not
          #160 — the strip is static legend text here. */}
      <div className="flex gap-4 px-4 text-xs text-ink-faint">
        <span>click — select a frame</span>
        <span>X — drop the selected frames</span>
        <span>double-click — open the loupe</span>
        <span>Esc — clear</span>
      </div>
    </div>
  );
}
