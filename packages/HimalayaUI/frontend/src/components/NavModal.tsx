import { useEffect, useMemo, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import { useExperiments, useSamples } from "../queries";
import type { Experiment, Sample } from "../api";
import { ModalShell } from "./ui";

const NAV_FIXTURE_EXPERIMENTS: { id: number; primary: string; secondary: string }[] = [
  { id: 1, primary: "Experiment A", secondary: "/data/lipids/expA" },
  { id: 2, primary: "Experiment B", secondary: "/data/lipids/expB" },
  { id: 3, primary: "Experiment C", secondary: "/data/lipids/expC" },
  { id: 4, primary: "Experiment D", secondary: "/data/lipids/expD" },
];

const NAV_FIXTURE_SAMPLES: { id: number; primary: string; secondary: string }[] = [
  { id: 1, primary: "DOPE 70%",        secondary: "JC001" },
  { id: 2, primary: "DOPE 80%",        secondary: "JC002" },
  { id: 3, primary: "DPPC 100%",       secondary: "JC003" },
  { id: 4, primary: "DPPC/DOPE 50/50", secondary: "JC004" },
];

function navFixtureItems(items: { id: number; primary: string; secondary: string }[]): JSX.Element {
  return (
    <>
      {items.map((item, idx) => (
        <div
          key={item.id}
          className={
            "w-full text-left px-3 py-2 flex flex-col gap-0.5 text-base " +
            (idx === 0 ? "bg-paper-sunk text-ink" : "text-ink")
          }
        >
          <span className="font-medium">{item.primary}</span>
          <span className="text-ink-soft text-sm font-sans">{item.secondary}</span>
        </div>
      ))}
    </>
  );
}

const NAV_EXPERIMENTS_FIXTURE = navFixtureItems(NAV_FIXTURE_EXPERIMENTS);
const NAV_SAMPLES_FIXTURE     = navFixtureItems(NAV_FIXTURE_SAMPLES);

/**
 * NavModal — cascading experiment → sample picker.
 *
 * Behavior:
 * - Opens with chips for whatever is already committed in the store.
 * - Step "experiment": filters the experiment list; Enter/Tab commits + advances to "sample".
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
      [s.name ?? "", s.display_name ?? ""].some(h => h.toLowerCase().includes(needle)),
    );
  }, [samplesQ.data, query]);

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
          primary: s.display_name || s.name || `Sample ${s.id}`,
          secondary: s.name && s.display_name && s.name !== s.display_name ? s.name : "",
        }));

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
      setSelIdx((i) => Math.min(activeList.length - 1, i + 1));
      return;
    }
    if (e.key === "ArrowUp") {
      e.preventDefault();
      setSelIdx((i) => Math.max(0, i - 1));
      return;
    }
    if (e.key === "Enter" || e.key === "Tab") {
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
    return `Sample ${samp?.display_name || samp?.name || pendingSamp}`;
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
                       placeholder:text-ink-faint font-sans text-base py-1"
            data-testid="nav-modal-input"
            spellCheck={false}
          />
          <span className="text-xs text-ink-faint px-1.5 py-0.5
                           border border-hair-strong rounded">esc</span>
        </div>

        <Skeleton
          name={step === "experiment" ? "nav-experiments" : "nav-samples"}
          className="flex-1 min-h-0 flex flex-col"
          loading={step === "experiment" ? experimentsQ.isLoading : samplesQ.isLoading}
          stagger={50}
          transition={200}
          fixture={step === "experiment" ? NAV_EXPERIMENTS_FIXTURE : NAV_SAMPLES_FIXTURE}
          fallback={<div className="px-4 py-6 text-center text-ink-soft italic text-base">{step === "experiment" ? "loading experiments…" : "loading samples…"}</div>}
        >
          <div className="flex-1 overflow-y-auto py-1" data-testid="nav-modal-results">
            {activeList.length === 0 ? (
              <div className="px-4 py-6 text-center text-ink-soft italic text-base">
                {step === "experiment"
                  ? "no experiments"
                  : pendingExp === undefined
                    ? "pick an experiment first"
                    : "no samples"}
              </div>
            ) : (
              activeList.map((item, idx) => (
                <button
                  key={`${step}-${item.id}`}
                  type="button"
                  data-testid={`nav-item-${step}-${item.id}`}
                  data-selected={idx === selIdx || undefined}
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
              ))
            )}
          </div>
        </Skeleton>

        <div className="flex items-center gap-3 px-3 py-2 border-t border-hair-strong
                        text-xs text-ink-faint">
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
      <button
        type="button"
        aria-label={`Remove ${label}`}
        onClick={onRemove}
        className="text-ink-soft hover:text-error px-0.5 leading-none"
        data-testid={testId ? `${testId}-remove` : undefined}
      >
        ×
      </button>
    </span>
  );
}
