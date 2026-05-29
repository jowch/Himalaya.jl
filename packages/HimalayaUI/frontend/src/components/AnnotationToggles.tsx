/**
 * AnnotationToggles — review-mode header toggles (Plan §Phase 9, Task 9.3;
 * spec §Plot rendering / "Annotation toggles").
 *
 * Two flags, both default `true`, both per-tab Zustand:
 *   - `showPeakTicks`  — render peak triangle marks above each trace.
 *   - `showPeakLabels` — render the `q` value label above each labeled peak.
 *
 * Both are AND-ed with the per-member `peak_display` state in
 * `MemberTraceLayer.buildMemberMarks` — turning the global flag off hides
 * everything; turning it on still respects the per-peak hide / label flags.
 *
 * **Hidden in edit mode.** Edit mode has its own peak-click cycle for
 * per-peak control, so the global flags would conflict with the
 * per-peak state the user is trying to manipulate. Caller (ComparePage vs.
 * ComparePageEdit) decides whether to mount this component.
 *
 * **Styling — text-link parity with PlotCard.** Mirrors the q-range-reset /
 * XScaleToggle vocabulary used by the Index page: `text-ink-faint` at rest,
 * `text-ink + bg-paper-sunk + border-hair-strong` on hover, `bg-paper-sunk text-ink`
 * when active. No native checkbox — `aria-pressed` carries the toggle
 * semantics and `data-active` lets E2E selectors assert state.
 *
 * No predicted-phase-ratio toggle. Per spec §Annotation toggles, v1
 * doesn't render predicted-q ticks at all — the figure is the result of
 * curation, not helpful suggestions. If users later need predicted-q
 * overlays, that ships behind a separate considered design.
 */
import { useAppState } from "../state";

interface ToggleButtonProps {
  testId: string;
  label: string;
  active: boolean;
  onToggle: () => void;
}

function ToggleButton({
  testId, label, active, onToggle,
}: ToggleButtonProps): JSX.Element {
  return (
    <button
      type="button"
      data-testid={testId}
      data-active={active ? "true" : "false"}
      aria-pressed={active}
      onClick={onToggle}
      className={[
        "px-1.5 py-0.5 rounded text-xs transition-colors",
        "border border-transparent hover:border-hair-strong",
        active
          ? "bg-paper-sunk text-ink"
          : "text-ink-faint hover:text-ink hover:bg-paper-sunk",
      ].join(" ")}
    >
      {label}
    </button>
  );
}

export function AnnotationToggles(): JSX.Element {
  const showPeakTicks  = useAppState((s) => s.showPeakTicks);
  const showPeakLabels = useAppState((s) => s.showPeakLabels);
  const setShowPeakTicks  = useAppState((s) => s.setShowPeakTicks);
  const setShowPeakLabels = useAppState((s) => s.setShowPeakLabels);

  return (
    <div
      data-testid="annotation-toggles"
      role="group"
      aria-label="Annotation toggles"
      className="inline-flex items-center gap-1 text-xs"
    >
      <ToggleButton
        testId="annotation-toggle-peaks"
        label="Peak ticks"
        active={showPeakTicks}
        onToggle={() => setShowPeakTicks(!showPeakTicks)}
      />
      <ToggleButton
        testId="annotation-toggle-labels"
        label="Peak labels"
        active={showPeakLabels}
        onToggle={() => setShowPeakLabels(!showPeakLabels)}
      />
    </div>
  );
}
