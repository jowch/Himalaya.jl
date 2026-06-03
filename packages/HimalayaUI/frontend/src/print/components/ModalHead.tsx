import type { ReactNode } from "react";
import { Kicker, IconButton } from "../ui";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface ModalHeadProps {
  kicker: ReactNode;
  title: ReactNode;
  titleId?: string;
  onClose: () => void;
  className?: string;
}

export function ModalHead({
  kicker,
  title,
  titleId,
  onClose,
  className,
}: ModalHeadProps): JSX.Element {
  return (
    <div
      data-testid="modal-head"
      className={cx(
        "flex items-start justify-between px-5 py-4 border-b border-hair",
        className,
      )}
    >
      <div>
        <Kicker tone="accent">{kicker}</Kicker>
        <h2
          {...(titleId ? { id: titleId } : {})}
          className="text-headline mt-0.5"
        >
          {title}
        </h2>
      </div>
      <IconButton
        dismiss
        label="Close"
        tone="ghost"
        onClick={onClose}
        className="shrink-0 -mr-1"
      />
    </div>
  );
}
