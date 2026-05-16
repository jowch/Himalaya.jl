/**
 * GroupingModeToggle — text-link toggle for the Compare-page grouping mode
 * (Plan §Phase 9, Task 9.2; spec §Trace coloring).
 *
 * Three options: By sample / By phase / Distinct. Accepts `mode` and
 * `onChange` as props; the parent (ComparePage) resolves the effective mode
 * via effectiveGroupingMode(draft, comparison) and routes the write to
 * setDraftViewGroupingMode (Compare UX C-4).
 *
 * Spec selectors: container `data-testid="grouping-mode"` with
 * `data-mode={"bySample" | "byPhase" | "distinct"}` reflecting the
 * active option (used by E2E tests). Each option carries `data-value` and
 * `data-active` for assertions.
 *
 * **Styling — text-link parity with PlotCard.** Mirrors the q-range-reset /
 * XScaleToggle vocabulary used by the Index page: `text-fg-dim` at rest,
 * `text-fg + bg-bg-hover + border-border` on hover, `bg-bg-subtle text-fg`
 * when active. No outer bordered wrapper, no leading "Color" cell — the
 * `aria-label="Trace grouping mode"` carries that semantic.
 */
import type { GroupingMode } from "../lib/comparison/coloring";

const OPTIONS: Array<{ value: GroupingMode; label: string }> = [
  { value: "bySample", label: "By sample" },
  { value: "byPhase",  label: "By phase"  },
  { value: "distinct", label: "Distinct"  },
];

export interface GroupingModeToggleProps {
  mode: GroupingMode;
  onChange: (mode: GroupingMode) => void;
}

export function GroupingModeToggle({
  mode,
  onChange,
}: GroupingModeToggleProps): JSX.Element {
  return (
    <div
      data-testid="grouping-mode"
      data-mode={mode}
      role="radiogroup"
      aria-label="Trace grouping mode"
      className="inline-flex items-center gap-1"
    >
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
            onClick={() => onChange(opt.value)}
            className={[
              "px-1.5 py-0.5 rounded text-xs transition-colors",
              "border border-transparent hover:border-border",
              active
                ? "bg-bg-subtle text-fg"
                : "text-fg-dim hover:text-fg hover:bg-bg-hover",
            ].join(" ")}
          >
            {opt.label}
          </button>
        );
      })}
    </div>
  );
}
