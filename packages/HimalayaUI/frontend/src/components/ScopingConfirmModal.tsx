import { ModalShell } from "./ui";
import { Button } from "./ui";

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
  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="sm"
      testId="scoping-confirm-modal"
      aria-label="Confirm series scoping"
      className="p-6"
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
        <Button type="button" data-testid="scoping-confirm-build" variant="solid"
          onClick={onConfirm} className="px-3 py-1.5 text-sm font-semibold">
          Confirm &amp; build →
        </Button>
      </div>
    </ModalShell>
  );
}
