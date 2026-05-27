import { useRef } from "react";
import { useFocusTrap } from "../hooks/useFocusTrap";

interface Props {
  open: boolean;
  orderingKey: string | undefined;
  count: number;
  onConfirm: () => void;
  onClose: () => void;
}

/**
 * ScopingConfirmModal — the confirm-and-build gate for series scoping (#174).
 * Focus-trapped (useFocusTrap) like NavModal/OnboardingFlow. Summarizes the
 * (key, N) sample_tags write; the build button is the single durable action.
 */
export function ScopingConfirmModal({
  open, orderingKey, count, onConfirm, onClose,
}: Props): JSX.Element | null {
  const dialogRef = useRef<HTMLDivElement>(null);
  useFocusTrap(dialogRef, open);
  if (!open) return null;
  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/40"
      onClick={onClose}
    >
      <div
        ref={dialogRef}
        data-testid="scoping-confirm-modal"
        role="dialog"
        aria-modal="true"
        aria-label="Confirm series scoping"
        className="rounded-lg bg-paper p-6 shadow-xl"
        onClick={(e) => e.stopPropagation()}
        onKeyDown={(e) => { if (e.key === "Escape") onClose(); }}
      >
        <h2 className="text-lg font-semibold text-ink">Confirm &amp; build</h2>
        <p className="mt-2 text-sm text-ink-soft">
          Records the ordering variable{" "}
          <span className="font-mono text-ink">{orderingKey ?? "—"}</span> on{" "}
          <span className="font-semibold text-ink">{count}</span>{" "}
          sample{count === 1 ? "" : "s"} as scoping tags.
        </p>
        <div className="mt-5 flex justify-end gap-2">
          <button type="button" data-testid="scoping-confirm-cancel"
            onClick={onClose}
            className="rounded px-3 py-1.5 text-sm text-ink-soft hover:underline">
            Cancel
          </button>
          <button type="button" data-testid="scoping-confirm-build"
            onClick={onConfirm}
            className="rounded bg-accent px-3 py-1.5 text-sm font-semibold text-paper">
            Confirm &amp; build →
          </button>
        </div>
      </div>
    </div>
  );
}
