import { ModalShell } from "../ui/ModalShell";
import { KbdLegend } from "./KbdLegend";
import { useAppState } from "../../state";
import type { ShortcutId } from "./shortcuts";

/**
 * KbdOverlay — the `?` keyboard shortcut legend modal. Opened by the
 * `helpOverlay` binding in `useGlobalShortcuts`; Esc closes via ModalShell's
 * built-in `closeOnEsc`. Lists the global/navigation keymap drawn directly from
 * the shortcut registry so the displayed keys can never drift from the live
 * bindings.
 *
 * Mounted app-wide in `PrintApp` (like NavModal) so it is reachable from any
 * surface.
 */

// The ids rendered in the overlay: the global + navigation + screen vocabulary.
// Surface-local in-surface verb keys (Focus P, Series A/⌘Enter) are NOT
// included per spec §4.1 — the overlay lists only the global/navigation keymap.
const OVERLAY_IDS: ShortcutId[] = [
  "prevSample",
  "nextSample",
  "prevDetail",
  "nextDetail",
  "openFocus",
  "openLoupe",
  "drop",
  "keep",
  "representative",
  "restore",
  "dismiss",
  "find",
  "helpOverlay",
];

export function KbdOverlay(): JSX.Element | null {
  const open = useAppState((s) => s.helpOverlayOpen);
  const close = useAppState((s) => s.closeHelpOverlay);

  return (
    <ModalShell
      open={open}
      onClose={close}
      size="sm"
      aria-label="Keyboard shortcuts"
      testId="kbd-overlay"
    >
      <div className="px-6 py-5">
        <KbdLegend ids={OVERLAY_IDS} testId="kbd-overlay-legend" />
      </div>
    </ModalShell>
  );
}
