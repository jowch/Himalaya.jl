import { useEffect, useId, useMemo, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { Button, Input, Popover } from "../ui";
import type { CorpusSample, Experiment } from "../../api";
import { safeScrollIntoView } from "../../lib/safeScrollIntoView";

export interface AddSamplePickerProps {
  /** Addable corpus samples (those not already in the recipe). */
  options: CorpusSample[];
  /** For experiment group headers + search-by-experiment. */
  experiments: Experiment[];
  onAdd: (sampleId: number) => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

function sampleLabel(s: CorpusSample): string {
  return s.name || `Sample ${s.id}`;
}

/**
 * AddSamplePicker (BU-PICKER) — replaces the Series Builder's ~130-option flat
 * native add-sample <select>. A search-first combobox in a Popover: type to
 * filter by sample name, `smp_{id}`, or experiment; results group under
 * experiment headers; every row carries its mono `smp_{id}` so duplicate
 * display names are distinguishable. Arrow keys move the active option (focus
 * stays in the input — the ARIA combobox pattern, aria-activedescendant), Enter
 * adds it, click adds it; the popover stays open for multi-add (the added
 * sample leaves the list as the recipe grows). Escape / outside-click close it
 * (Popover) and return focus to the trigger.
 */
export function AddSamplePicker({
  options,
  experiments,
  onAdd,
  className,
}: AddSamplePickerProps): JSX.Element {
  const [query, setQuery] = useState("");
  const [activeId, setActiveId] = useState<number | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);
  const listboxId = useId();

  const expName = useMemo(() => {
    const m = new Map<number, string>();
    for (const e of experiments) m.set(e.id, e.name ?? `Experiment ${e.id}`);
    return m;
  }, [experiments]);

  const q = query.trim().toLowerCase();
  const filtered = useMemo(
    () =>
      options.filter((s) => {
        if (q === "") return true;
        const name = sampleLabel(s).toLowerCase();
        const exp = (expName.get(s.experiment_id) ?? "").toLowerCase();
        return (
          name.includes(q) || `smp_${s.id}`.includes(q) || exp.includes(q)
        );
      }),
    [options, q, expName],
  );

  // Group by experiment, ordered by the experiments list (then any leftovers).
  const groups = useMemo(() => {
    const byExp = new Map<number, CorpusSample[]>();
    for (const s of filtered) {
      const arr = byExp.get(s.experiment_id);
      if (arr) arr.push(s);
      else byExp.set(s.experiment_id, [s]);
    }
    const ordered = [
      ...experiments.map((e) => e.id).filter((id) => byExp.has(id)),
      ...[...byExp.keys()].filter((id) => !experiments.some((e) => e.id === id)),
    ];
    return ordered.map((id) => ({
      expId: id,
      name: expName.get(id) ?? `Experiment ${id}`,
      samples: byExp.get(id)!,
    }));
  }, [filtered, experiments, expName]);

  // Flat order (grouped) for arrow-key traversal.
  const flat = useMemo(() => groups.flatMap((g) => g.samples), [groups]);

  // Keep the active option valid as the list narrows/grows (e.g. after a query
  // change or an add). Default to the first option; clear when empty.
  useEffect(() => {
    setActiveId((cur) =>
      flat.length === 0
        ? null
        : cur !== null && flat.some((s) => s.id === cur)
          ? cur
          : flat[0]!.id,
    );
  }, [flat]);

  function moveActive(delta: 1 | -1): void {
    if (flat.length === 0) return;
    const i = flat.findIndex((s) => s.id === activeId);
    const next = i < 0 ? 0 : (i + delta + flat.length) % flat.length;
    const id = flat[next]!.id;
    setActiveId(id);
    // Keep the active option in view (no-op under jsdom — see safeScrollIntoView).
    safeScrollIntoView(document.getElementById(`add-opt-${id}`));
  }

  function add(id: number): void {
    onAdd(id);
    // Stay open for multi-add; the added sample drops out as `options` shrinks.
    inputRef.current?.focus();
  }

  function onInputKeyDown(e: KeyboardEvent<HTMLInputElement>): void {
    if (e.key === "ArrowDown") {
      e.preventDefault();
      moveActive(1);
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      moveActive(-1);
    } else if (e.key === "Enter") {
      if (activeId !== null) {
        e.preventDefault();
        add(activeId);
      }
    }
  }

  return (
    <Popover
      label="Add a sample to the series"
      initialFocusRef={inputRef}
      className={`w-72${className ? ` ${className}` : ""}`}
      trigger={
        <Button
          variant="outline"
          data-testid="builder-add-sample"
          className="w-full"
        >
          + Add sample…
        </Button>
      }
    >
      <div className="flex flex-col gap-2">
        <Input
          testId="add-sample-search"
          inputRef={inputRef}
          value={query}
          onValueChange={setQuery}
          placeholder="Search samples…"
          aria-label="Search samples to add"
          role="combobox"
          aria-expanded
          aria-controls={listboxId}
          aria-autocomplete="list"
          {...(activeId !== null
            ? { "aria-activedescendant": `add-opt-${activeId}` }
            : {})}
          onKeyDown={onInputKeyDown}
        />
        <div
          id={listboxId}
          role="listbox"
          aria-label="Addable samples"
          data-testid="add-sample-listbox"
          className="max-h-64 overflow-y-auto -mx-1 px-1"
        >
          {flat.length === 0 ? (
            <div className="px-2 py-3 text-caption text-ink-soft">
              No samples match.
            </div>
          ) : (
            groups.map((g) => (
              <div key={g.expId} role="group" aria-label={g.name}>
                <div className="px-2 pt-2 pb-1 text-caption text-ink-soft truncate">
                  {g.name}
                </div>
                {g.samples.map((s) => (
                  <button
                    key={s.id}
                    type="button"
                    id={`add-opt-${s.id}`}
                    role="option"
                    aria-selected={s.id === activeId}
                    data-testid={`add-opt-${s.id}`}
                    onMouseEnter={() => setActiveId(s.id)}
                    onClick={() => add(s.id)}
                    className={`w-full text-left rounded-sm px-2 py-1.5 flex items-center gap-2 text-meta ${
                      s.id === activeId ? "bg-paper-sunk text-ink" : "text-ink"
                    }`}
                  >
                    <span className="flex-1 min-w-0 truncate">
                      {sampleLabel(s)}
                    </span>
                    <span className="flex-shrink-0 font-mono text-caption text-ink-soft">
                      smp_{s.id}
                    </span>
                  </button>
                ))}
              </div>
            ))
          )}
        </div>
      </div>
    </Popover>
  );
}
