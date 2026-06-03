import type { ReactNode } from "react";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface ModalFooterProps {
  note?: ReactNode;
  actions: ReactNode;
  className?: string;
}

export function ModalFooter({ note, actions, className }: ModalFooterProps): JSX.Element {
  return (
    <div
      data-testid="modal-foot"
      className={cx("flex items-center justify-between px-5 pt-3.5 pb-4", className)}
    >
      {note != null ? (
        <span className="text-caption text-ink-faint">{note}</span>
      ) : (
        <span />
      )}
      <div className="flex gap-2">{actions}</div>
    </div>
  );
}
