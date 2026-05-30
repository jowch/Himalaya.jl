/**
 * GroupingModeToggle — text-link toggle for the Compare-page grouping mode
 * (Plan §Phase 9, Task 9.2; spec §Trace coloring).
 *
 * Three options: By sample / By phase / Distinct. Accepts `value` and
 * `onChange` as props (unified with the sibling segmented controls
 * `ScaleToggle` and `RepresentationToggle` after PR #251 round 1); the parent
 * (ComparePage) resolves the effective mode via
 * effectiveGroupingMode(draft, comparison) and routes the write to
 * setDraftViewGroupingMode (Compare UX C-4).
 *
 * Spec selectors: container `data-testid="grouping-mode"` with
 * `data-mode={"bySample" | "byPhase" | "distinct"}` reflecting the
 * active option (used by E2E tests). Each option carries `data-value` and
 * `data-active` for assertions.
 *
 * Thin wrapper over SegmentedControl<GroupingMode> with `role="radiogroup"` +
 * `variant="plain"` + `stateAttr="data-mode"` (aliases `data-state` →
 * `data-mode` on the container). Gains roving arrow-key nav (WAI-ARIA
 * radiogroup pattern). Active fill is canonical ink-on-paper (`bg-ink text-paper`).
 */
import { SegmentedControl, type SegmentOption } from "./ui/SegmentedControl";
import type { GroupingMode } from "../lib/comparison/coloring";

const OPTIONS: ReadonlyArray<SegmentOption<GroupingMode>> = [
  { value: "bySample", label: "By sample" },
  { value: "byPhase",  label: "By phase"  },
  { value: "distinct", label: "Distinct"  },
];

export interface GroupingModeToggleProps {
  value: GroupingMode;
  onChange: (next: GroupingMode) => void;
}

export function GroupingModeToggle({
  value,
  onChange,
}: GroupingModeToggleProps): JSX.Element {
  return (
    <SegmentedControl<GroupingMode>
      aria-label="Trace grouping mode"
      role="radiogroup"
      variant="plain"
      testId="grouping-mode"
      stateAttr="data-mode"
      options={OPTIONS}
      value={value}
      onChange={onChange}
    />
  );
}
