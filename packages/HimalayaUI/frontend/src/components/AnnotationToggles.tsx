/**
 * AnnotationToggles — review-mode header checkboxes (Plan §Phase 9, Task 9.3;
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
 * No predicted-phase-ratio toggle. Per spec §Annotation toggles, v1
 * doesn't render predicted-q ticks at all — the figure is the result of
 * curation, not helpful suggestions. If users later need predicted-q
 * overlays, that ships behind a separate considered design.
 */
import { useAppState } from "../state";

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
      className="inline-flex items-center gap-3 text-xs text-fg-muted"
    >
      <label className="inline-flex items-center gap-1 cursor-pointer">
        <input
          type="checkbox"
          data-testid="annotation-toggle-peaks"
          checked={showPeakTicks}
          onChange={(e) => setShowPeakTicks(e.currentTarget.checked)}
          className="cursor-pointer"
        />
        <span>Peak ticks</span>
      </label>
      <label className="inline-flex items-center gap-1 cursor-pointer">
        <input
          type="checkbox"
          data-testid="annotation-toggle-labels"
          checked={showPeakLabels}
          onChange={(e) => setShowPeakLabels(e.currentTarget.checked)}
          className="cursor-pointer"
        />
        <span>Peak labels</span>
      </label>
    </div>
  );
}
