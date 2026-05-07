/**
 * GroupingModeToggle — segmented control for the Compare-page grouping mode
 * (Plan §Phase 9, Task 9.2; spec §Trace coloring).
 *
 * Three options: By sample / By phase / Distinct. Reads/writes Zustand
 * `groupingMode` (per-tab; not persisted on the comparison or in storage).
 *
 * Spec selectors: container `data-testid="grouping-mode"` with
 * `data-mode={"bySample" | "byPhase" | "distinct"}` reflecting the
 * active option (used by E2E tests).
 */
import { useAppState } from "../state";
import type { GroupingMode } from "../lib/comparison/coloring";

const OPTIONS: Array<{ value: GroupingMode; label: string }> = [
  { value: "bySample", label: "By sample" },
  { value: "byPhase",  label: "By phase"  },
  { value: "distinct", label: "Distinct"  },
];

export function GroupingModeToggle(): JSX.Element {
  const mode = useAppState((s) => s.groupingMode);
  const setMode = useAppState((s) => s.setGroupingMode);

  return (
    <div
      data-testid="grouping-mode"
      data-mode={mode}
      role="radiogroup"
      aria-label="Trace grouping mode"
      className="inline-flex items-center gap-0 rounded border border-border overflow-hidden"
    >
      <span className="px-2 py-1 text-xs text-fg-muted bg-bg-elevated/40 border-r border-border">
        Color
      </span>
      {OPTIONS.map((opt) => {
        const active = opt.value === mode;
        return (
          <button
            key={opt.value}
            type="button"
            role="radio"
            aria-checked={active}
            data-active={active ? "true" : "false"}
            data-value={opt.value}
            onClick={() => setMode(opt.value)}
            className={
              "px-2 py-1 text-xs " +
              (active
                ? "bg-accent text-bg"
                : "text-fg-muted hover:text-fg hover:bg-bg-elevated/40")
            }
          >
            {opt.label}
          </button>
        );
      })}
    </div>
  );
}
