import { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../../state";
import { useCorpusSamples, useExperiments, useSamples } from "../../queries";
import type { CorpusSample, Experiment, Sample } from "../../api";
import { IconButton, ModalShell } from "../ui";

function navSkeletonRows(n: number): JSX.Element {
  return (
    <>
      {Array.from({ length: n }, (_, idx) => (
        <div
          key={idx}
          className={
            "w-full text-left px-3 py-2 flex flex-col gap-0.5 text-base " +
            (idx === 0 ? "bg-paper-sunk text-ink" : "text-ink")
          }
        >
          <span className="font-medium">&nbsp;</span>
          <span className="text-ink-soft text-sm font-sans">&nbsp;</span>
        </div>
      ))}
    </>
  );
}

// Cap for the direct-sample group at the experiment step (SA-F4). The corpus is
// ~139 samples; anything past the cap is disclosed honestly as "+N more".
const SAMPLE_HIT_CAP = 8;

/**
 * NavModal — cascading experiment → sample picker.
 *
 * Behavior:
 * - Opens with chips for whatever is already committed in the store.
 * - Step "experiment": filters the experiment list; Enter/Tab commits + advances to "sample".
 *   A query that also matches sample names across the whole corpus surfaces a
 *   capped "Samples" group under the experiment matches (SA-F4) — selecting one
 *   skips the cascade and lands straight on /sample/:id. One flat highlight
 *   order spans both groups, experiments first.
 * - Step "sample": filters samples in the chosen experiment; Enter/Tab commits + closes modal.
 * - Backspace on empty input rewinds one step (removes sample chip, then experiment chip).
 * - Clicking a chip × is equivalent to Backspace at that chip's position.
 * - Esc closes without committing further changes.
 *
 * All state reads/writes go through `useAppState`, so the modal is self-contained.
 */
export function NavModal(): JSX.Element | null {
  const open          = useAppState((s) => s.navModalOpen);
  const step          = useAppState((s) => s.navModalStep);
  const committedExp  = useAppState((s) => s.activeExperimentId);
  const committedSamp = useAppState((s) => s.activeSampleId);
  const closeModal    = useAppState((s) => s.closeNavModal);
  const setStep       = useAppState((s) => s.setNavModalStep);
  const setExperiment = useAppState((s) => s.setActiveExperiment);
  const setSample     = useAppState((s) => s.setActiveSample);
  const navigate      = useNavigate();

  // Local "pending" selection inside the modal — lets us rewind without nuking committed state
  // until the user explicitly commits.
  const [pendingExp, setPendingExp]   = useState<number | undefined>(committedExp);
  const [pendingSamp, setPendingSamp] = useState<number | undefined>(committedSamp);
  const [query, setQuery] = useState("");
  const [selIdx, setSelIdx] = useState(0);

  const inputRef  = useRef<HTMLInputElement>(null);

  // When the modal opens, sync pending values with committed, reset query.
  useEffect(() => {
    if (open) {
      setPendingExp(committedExp);
      setPendingSamp(committedSamp);
      setQuery("");
      setSelIdx(0);
      // Focus the input. Synchronous call works in jsdom; rAF is unreliable there.
      inputRef.current?.focus();
    }
  }, [open, committedExp, committedSamp]);

  const experimentsQ = useExperiments();
  const samplesQ     = useSamples(pendingExp ?? 0);
  const corpusQ      = useCorpusSamples();

  const filteredExperiments: Experiment[] = useMemo(() => {
    const list = experimentsQ.data ?? [];
    if (!query) return list;
    const needle = query.toLowerCase();
    return list.filter((e) =>
      (e.name ?? "").toLowerCase().includes(needle) ||
      e.path.toLowerCase().includes(needle),
    );
  }, [experimentsQ.data, query]);

  const filteredSamples: Sample[] = useMemo(() => {
    const list = samplesQ.data ?? [];
    if (!query) return list;
    const needle = query.toLowerCase();
    return list.filter((s) =>
      s.name.toLowerCase().includes(needle),
    );
  }, [samplesQ.data, query]);

  // SA-F4: while at the experiment step, a non-empty query also matches sample
  // names across the WHOLE corpus (same matcher as the sample step), so a known
  // sample name is one keystroke sequence away — no experiment commit needed.
  const corpusSampleHits: CorpusSample[] = useMemo(() => {
    if (step !== "experiment" || !query) return [];
    const needle = query.toLowerCase();
    return (corpusQ.data ?? []).filter((s) =>
      s.name.toLowerCase().includes(needle),
    );
  }, [step, query, corpusQ.data]);

  const visibleSampleHits = corpusSampleHits.slice(0, SAMPLE_HIT_CAP);
  const sampleHitOverflow = corpusSampleHits.length - visibleSampleHits.length;

  // Reset selection cursor on query/step change
  useEffect(() => { setSelIdx(0); }, [query, step]);

  const activeList: readonly { id: number; primary: string; secondary: string }[] =
    step === "experiment"
      ? filteredExperiments.map((e) => ({
          id: e.id,
          primary: e.name ?? `Experiment ${e.id}`,
          secondary: e.path,
        }))
      : filteredSamples.map((s) => ({
          id: s.id,
          primary: s.name || `Sample ${s.id}`,
          secondary: "",
        }));

  // One flat highlight order across both groups, experiments first (SA-F4).
  const totalRows = activeList.length + visibleSampleHits.length;

  const experimentName = (id: number): string => {
    const exp = experimentsQ.data?.find((e) => e.id === id);
    return exp?.name ?? `Experiment ${id}`;
  };

  const commitExperiment = (id: number): void => {
    setPendingExp(id);
    setPendingSamp(undefined);
    setStep("sample");
    setQuery("");
  };

  const commitSample = (id: number): void => {
    // write both to the store in the right order
    if (pendingExp !== undefined && pendingExp !== committedExp) {
      setExperiment(pendingExp);
    }
    setSample(id);
    closeModal();
    // M1 on-ramp: committing a sample is a door into the indexing workspace.
    // Before this the picker only mutated the store and left the URL where it
    // was; now it actually lands you on /sample/:id (the third door, beside the
    // contact-sheet status cell and the loupe "Open in the Index stage" link).
    navigate(`/sample/${id}`);
  };

  // SA-F4: a direct sample hit from the experiment step skips the cascade and
  // lands on exactly the destination the cascading flow ends at — /sample/:id —
  // with the sample's own experiment committed to the store.
  const commitDirectSample = (s: CorpusSample): void => {
    if (s.experiment_id !== committedExp) {
      setExperiment(s.experiment_id);
    }
    setSample(s.id);
    closeModal();
    navigate(`/sample/${s.id}`);
  };

  const popSampleChip = (): void => {
    setPendingSamp(undefined);
    setStep("sample");
    setQuery("");
  };

  const popExperimentChip = (): void => {
    setPendingExp(undefined);
    setPendingSamp(undefined);
    setStep("experiment");
    setQuery("");
  };

  const onInputKeyDown = (e: React.KeyboardEvent<HTMLInputElement>): void => {
    if (e.key === "Escape") {
      e.preventDefault();
      closeModal();
      return;
    }
    if (e.key === "Backspace" && query === "") {
      e.preventDefault();
      // Pop chips one at a time, right-to-left: sample first, then experiment.
      if (pendingSamp !== undefined) {
        popSampleChip();
      } else if (pendingExp !== undefined) {
        popExperimentChip();
      }
      return;
    }
    if (e.key === "ArrowDown") {
      e.preventDefault();
      setSelIdx((i) => Math.min(totalRows - 1, i + 1));
      return;
    }
    if (e.key === "ArrowUp") {
      e.preventDefault();
      setSelIdx((i) => Math.max(0, i - 1));
      return;
    }
    if (e.key === "Enter" || e.key === "Tab") {
      if (step === "experiment" && selIdx >= activeList.length) {
        // Cursor sits in the Samples group (SA-F4) — direct sample commit.
        const hit = visibleSampleHits[selIdx - activeList.length];
        if (!hit) return;
        e.preventDefault();
        commitDirectSample(hit);
        return;
      }
      const picked = activeList[selIdx];
      if (!picked) return;
      e.preventDefault();
      if (step === "experiment") commitExperiment(picked.id);
      else                        commitSample(picked.id);
    }
  };

  const expChipLabel = (() => {
    if (pendingExp === undefined) return null;
    const exp = experimentsQ.data?.find((e) => e.id === pendingExp);
    return `Experiment ${exp?.name ?? pendingExp}`;
  })();

  const sampChipLabel = (() => {
    if (pendingSamp === undefined) return null;
    const samp = samplesQ.data?.find((s) => s.id === pendingSamp);
    return `Sample ${samp?.name || pendingSamp}`;
  })();

  return (
    <ModalShell
      open={open}
      onClose={closeModal}
      size="md"
      align="top"
      closeOnEsc={false}
      testId="nav-modal"
      aria-label="Navigate to experiment or sample"
      className="max-h-[72vh]"
    >
        <div className="flex items-center gap-2 flex-wrap px-3 py-2.5 border-b border-hair-strong">
          {expChipLabel && (
            <Chip
              label={expChipLabel}
              onRemove={popExperimentChip}
              testId="nav-chip-experiment"
            />
          )}
          {sampChipLabel && (
            <Chip
              label={sampChipLabel}
              onRemove={popSampleChip}
              testId="nav-chip-sample"
            />
          )}
          <input
            ref={inputRef}
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={onInputKeyDown}
            placeholder={step === "experiment" ? "find experiment…" : "find sample…"}
            className="flex-1 min-w-[120px] bg-transparent border-0 outline-0 text-ink
                       placeholder:text-ink-soft font-sans text-base py-1"
            data-testid="nav-modal-input"
            spellCheck={false}
          />
          <span className="text-xs text-ink-soft px-1.5 py-0.5
                           border border-hair-strong rounded">esc</span>
        </div>

        <Skeleton
          name={step === "experiment" ? "nav-experiments" : "nav-samples"}
          className="flex-1 min-h-0 flex flex-col"
          loading={step === "experiment" ? experimentsQ.isLoading : samplesQ.isLoading}
          stagger={50}
          transition={200}
          fixture={navSkeletonRows(4)}
          fallback={<div className="px-4 py-6 text-center text-ink-soft italic text-base">{step === "experiment" ? "loading experiments…" : "loading samples…"}</div>}
        >
          <div className="flex-1 overflow-y-auto py-1" data-testid="nav-modal-results">
            {activeList.length === 0 && visibleSampleHits.length === 0 ? (
              <div className="px-4 py-6 text-center text-ink-soft italic text-base">
                {step === "experiment"
                  ? "no experiments"
                  : pendingExp === undefined
                    ? "pick an experiment first"
                    : "no samples"}
              </div>
            ) : (
              <>
                {activeList.map((item, idx) => (
                  <button
                    key={`${step}-${item.id}`}
                    type="button"
                    data-testid={`nav-item-${step}-${item.id}`}
                    data-selected={idx === selIdx || undefined}
                    // UI-RINGCLIP: full-bleed row in an overflow-y-auto list —
                    // inset the focus ring so its left/right edges aren't clipped.
                    data-focus-ring="inset"
                    onMouseEnter={() => setSelIdx(idx)}
                    onClick={() => {
                      if (step === "experiment") commitExperiment(item.id);
                      else                        commitSample(item.id);
                    }}
                    className={
                      "w-full text-left px-3 py-2 flex flex-col gap-0.5 text-base " +
                      (idx === selIdx ? "bg-paper-sunk text-ink" : "text-ink hover:bg-paper-sunk")
                    }
                  >
                    <span className="font-medium">{item.primary}</span>
                    {item.secondary && (
                      <span className="text-ink-soft text-sm font-sans">{item.secondary}</span>
                    )}
                  </button>
                ))}
                {visibleSampleHits.length > 0 && (
                  <>
                    <div
                      data-testid="nav-samples-group-label"
                      className="px-3 pt-2 pb-1 text-xs text-ink-soft"
                    >
                      Samples
                    </div>
                    {visibleSampleHits.map((s, i) => {
                      const idx = activeList.length + i;
                      return (
                        <button
                          key={`corpus-sample-${s.id}`}
                          type="button"
                          data-testid={`nav-item-corpus-sample-${s.id}`}
                          data-selected={idx === selIdx || undefined}
                          data-focus-ring="inset"
                          onMouseEnter={() => setSelIdx(idx)}
                          onClick={() => commitDirectSample(s)}
                          className={
                            "w-full text-left px-3 py-2 flex flex-col gap-0.5 text-base " +
                            (idx === selIdx ? "bg-paper-sunk text-ink" : "text-ink hover:bg-paper-sunk")
                          }
                        >
                          <span className="font-medium">
                            {s.name || `Sample ${s.id}`}
                          </span>
                          <span className="text-ink-soft text-sm font-sans">
                            {experimentName(s.experiment_id)}
                          </span>
                        </button>
                      );
                    })}
                    {sampleHitOverflow > 0 && (
                      <div
                        data-testid="nav-samples-overflow"
                        className="px-3 py-2 text-ink-soft italic text-sm"
                      >
                        +{sampleHitOverflow} more {sampleHitOverflow === 1 ? "match" : "matches"}
                      </div>
                    )}
                  </>
                )}
              </>
            )}
          </div>
        </Skeleton>

        <div className="flex items-center gap-3 px-3 py-2 border-t border-hair-strong
                        text-xs text-ink-soft">
          <span><kbd className="border border-hair-strong rounded px-1">↑↓</kbd> navigate</span>
          <span><kbd className="border border-hair-strong rounded px-1">⏎</kbd> select</span>
          <span><kbd className="border border-hair-strong rounded px-1">⌫</kbd> back</span>
          <span className="flex-1" />
          <span>{step === "experiment" ? "experiment" : "sample"}</span>
        </div>
    </ModalShell>
  );
}

interface ChipProps {
  label: string;
  onRemove: () => void;
  testId?: string;
}

function Chip({ label, onRemove, testId }: ChipProps): JSX.Element {
  return (
    <span
      data-testid={testId}
      className="inline-flex items-center gap-1 px-2 py-1 rounded-md
                 bg-paper-sunk border border-hair-strong text-sm text-ink"
    >
      {label}
      <IconButton
        label={`Remove ${label}`}
        dismiss
        tone="accent"
        onClick={onRemove}
        data-testid={testId ? `${testId}-remove` : undefined}
      />
    </span>
  );
}
