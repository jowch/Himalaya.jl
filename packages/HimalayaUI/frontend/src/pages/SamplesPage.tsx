import { useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSamples,
  useCorpusExposures,
  useExperiments,
  useScreenedProgress,
} from "../queries";
import { HintText } from "../components/ui/HintText";
import {
  ContactSheetRow,
  CONTACT_SHEET_COLS,
} from "../components/ContactSheetRow";

/** Footer keyboard legend (L-5/L-8/L-11) — five keycap-chip hints. */
const KB_HINTS: ReadonlyArray<{ key: string; label: string }> = [
  { key: "click", label: "select a frame" },
  { key: "⇧ click", label: "extend the range" },
  { key: "X", label: "drop the selected frames" },
  { key: "double-click", label: "open the loupe" },
  { key: "Esc", label: "clear" },
];

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

  // L-3: the title is the beamtime (or the whole corpus when unfiltered).
  const title =
    beamtime === undefined ? "The corpus" : scopeName;

  // M-1: screened-progress aggregate over the visible samples.
  const { screened, total } = useScreenedProgress(filtered);
  const pct = total === 0 ? 0 : Math.round((screened / total) * 100);

  // R2-M14: one bulk hook drives every row's exposures, so the page-level
  // observer count is O(1) instead of O(samples) — no `ERR_INSUFFICIENT_
  // RESOURCES` from the 139× JSON fan-out. The cache rows are the same
  // `queryKeys.exposures(id)` slots `useExposures` reads from, so loupe +
  // mutator paths inherit these rows for free (no double-fetch).
  const corpusExposures = useCorpusExposures(filtered);

  return (
    <div data-testid="samples-page" className="flex flex-col gap-4 p-6">
      {/* L-3: kicker + beamtime serif h1 + descriptive sub, with the M-1
          progress block on the right. */}
      <header className="flex items-end justify-between gap-8">
        <div className="flex flex-col gap-1">
          <div className="text-xs font-semibold uppercase tracking-[0.14em] text-print-accent">
            Contact sheet
          </div>
          <h1 data-testid="samples-title" className="text-display text-ink">
            {title}
          </h1>
          <p
            data-testid="samples-sub"
            className="mt-1 max-w-[62ch] text-sm text-ink-soft"
          >
            Flip the frames and drop the ones with flares or artifacts. Tags
            are a light, optional note on what each sample is — the ordering
            variable is named later, when you scope a series.
          </p>
          {/* Kept for back-compat with the scope assertion; doubles as a
              screen-reader-friendly statement of the active filter. */}
          <p data-testid="samples-scope" className="sr-only">
            Showing {scopeName}.
          </p>
        </div>

        <div
          data-testid="screened-progress"
          className="flex shrink-0 flex-col items-end text-right"
        >
          <div className="text-display !text-[25px] text-ink">
            {screened}
            <b className="font-medium text-ink-faint"> / {total}</b>
          </div>
          <div className="mt-0.5 text-[10.5px] font-bold uppercase tracking-[0.08em] text-ink-faint">
            samples screened
          </div>
          <div
            data-testid="screened-progress-bar"
            className="mt-2 h-1 w-[152px] overflow-hidden rounded-full bg-hair"
          >
            <div
              className="h-full bg-print-accent"
              style={{ width: `${pct}%` }}
            />
          </div>
        </div>
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
          {filtered.length === 0 ? (
            <div
              data-testid="samples-empty"
              className="px-4 py-8 text-sm text-ink-faint"
            >
              {beamtime === undefined
                ? "No samples in the corpus yet."
                : "No samples in this beamtime."}
            </div>
          ) : (
            <div data-testid="contact-sheet-rows">
              {filtered.map((s) => (
                <ContactSheetRow
                  key={s.id}
                  sample={s}
                  exposures={corpusExposures.byId.get(s.id)}
                  exposuresLoading={
                    !corpusExposures.byId.has(s.id) && corpusExposures.isLoading
                  }
                />
              ))}
            </div>
          )}
        </Skeleton>
      )}

      {/* L-5/L-8/L-11: footer keycap-chip legend. Every hint is a live
          control on each ContactSheetRow: click / shift-click select &
          extend a contiguous range over a sample's frames, double-click
          opens that sample's loupe, and while a selection exists X
          batch-rejects it and Esc clears it.

          R2-M15: sticky-bottom so the legend stays in view on long scrolls
          (the corpus has ~139 samples). The CullBar already carries the X /
          Esc keycaps for live selections, so this legend is purely the
          rest-state "what you can do here" anchor. */}
      <div
        data-testid="kb-legend"
        className="sticky bottom-0 z-10 -mx-6 -mb-6 mt-2 flex flex-wrap gap-5
                   border-t border-hair bg-paper/90 px-10 py-2 text-[11.5px]
                   text-ink-faint backdrop-blur-sm"
      >
        {KB_HINTS.map((h) => (
          <span key={h.label}>
            <span
              data-kb-key
              className="mr-1.5 rounded border border-hair-strong border-b-2
                         bg-plate px-1.5 py-px font-mono text-[10.5px] text-ink-soft"
            >
              {h.key}
            </span>
            {h.label}
          </span>
        ))}
      </div>
    </div>
  );
}
