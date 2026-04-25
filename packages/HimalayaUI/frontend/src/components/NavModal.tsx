import { useEffect, useMemo, useRef, useState } from "react";
import { useAppState } from "../state";
import { useExperiments, useSamples } from "../queries";
import type { Experiment, Sample } from "../api";
import { useFocusTrap } from "../hooks/useFocusTrap";

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

  // Local "pending" selection inside the modal — lets us rewind without nuking committed state
  // until the user explicitly commits.
  const [pendingExp, setPendingExp]   = useState<number | undefined>(committedExp);
  const [pendingSamp, setPendingSamp] = useState<number | undefined>(committedSamp);
  const [query, setQuery] = useState("");
  const [selIdx, setSelIdx] = useState(0);

  const inputRef  = useRef<HTMLInputElement>(null);
  const dialogRef = useRef<HTMLDivElement>(null);
  useFocusTrap(dialogRef, open);

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
      (s.name  ?? "").toLowerCase().includes(needle) ||
      (s.label ?? "").toLowerCase().includes(needle),
    );
  }, [samplesQ.data, query]);

  // Reset selection cursor on query/step change
  useEffect(() => { setSelIdx(0); }, [query, step]);

  if (!open) return null;

  const activeList: readonly { id: number; primary: string; secondary: string }[] =
    step === "experiment"
      ? filteredExperiments.map((e) => ({
          id: e.id,
          primary: e.name ?? `Experiment ${e.id}`,
          secondary: e.path,
        }))
      : filteredSamples.map((s) => ({
          id: s.id,
          primary: s.name ?? s.label ?? `Sample ${s.id}`,
          secondary: s.label && s.name && s.label !== s.name ? s.label : "",
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
    return `Sample ${samp?.name ?? samp?.label ?? pendingSamp}`;
  })();

  return (
    <div
      data-testid="nav-modal"
      className="fixed inset-0 z-50 flex items-start justify-center pt-[12vh]
                 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm
                 anim-pal-in"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) closeModal(); }}
    >
      <div
        ref={dialogRef}
        role="dialog"
        aria-modal="true"
        className="w-[min(640px,calc(100vw-48px))] max-h-[72vh]
                   bg-bg-elevated border border-border rounded-xl shadow-2xl
                   flex flex-col overflow-hidden
                   anim-pal-scale"
      >
        <div className="flex items-center gap-2 flex-wrap px-3 py-2.5 border-b border-border">
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
            className="flex-1 min-w-[120px] bg-transparent border-0 outline-0 text-fg
                       placeholder:text-fg-dim font-sans text-base py-1"
            data-testid="nav-modal-input"
            spellCheck={false}
          />
          <span className="text-[10px] text-fg-dim px-1.5 py-0.5
                           border border-border rounded">esc</span>
        </div>

        <div className="flex-1 overflow-y-auto py-1" data-testid="nav-modal-results">
          {activeList.length === 0 ? (
            <div className="px-4 py-6 text-center text-fg-muted italic text-[13px]">
              {step === "experiment"
                ? (experimentsQ.isPending ? "loading experiments…" : "no experiments")
                : pendingExp === undefined
                  ? "pick an experiment first"
                  : (samplesQ.isPending ? "loading samples…" : "no samples")}
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
                  "w-full text-left px-3 py-2 flex flex-col gap-0.5 text-[13px] " +
                  (idx === selIdx ? "bg-bg-hover text-fg" : "text-fg hover:bg-bg-hover")
                }
              >
                <span className="font-medium">{item.primary}</span>
                {item.secondary && (
                  <span className="text-fg-muted text-[11px] font-sans">{item.secondary}</span>
                )}
              </button>
            ))
          )}
        </div>

        <div className="flex items-center gap-3 px-3 py-2 border-t border-border
                        text-[10px] text-fg-dim">
          <span><kbd className="border border-border rounded px-1">↑↓</kbd> navigate</span>
          <span><kbd className="border border-border rounded px-1">⏎</kbd> select</span>
          <span><kbd className="border border-border rounded px-1">⌫</kbd> back</span>
          <span className="flex-1" />
          <span>{step === "experiment" ? "experiment" : "sample"}</span>
        </div>
      </div>
    </div>
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
                 bg-bg-hover border border-border text-[11.5px] text-fg"
    >
      {label}
      <button
        type="button"
        aria-label={`Remove ${label}`}
        onClick={onRemove}
        className="text-fg-muted hover:text-error px-0.5 leading-none"
        data-testid={testId ? `${testId}-remove` : undefined}
      >
        ×
      </button>
    </span>
  );
}
