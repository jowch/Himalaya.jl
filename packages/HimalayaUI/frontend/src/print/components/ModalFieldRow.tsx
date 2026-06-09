import type { ReactNode } from "react";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface ModalFieldRowProps {
  label: ReactNode;
  labelSuffix?: ReactNode;
  children: ReactNode;
  className?: string;
}

export function ModalFieldRow({
  label,
  labelSuffix,
  children,
  className,
}: ModalFieldRowProps): JSX.Element {
  return (
    <div
      data-testid="modal-field-row"
      className={cx("flex items-center gap-3.5", className)}
    >
      <span
        data-testid="modal-field-label"
        className="text-label shrink-0"
        style={{ width: 100 }}
      >
        {label}
        {labelSuffix != null && (
          <span className="font-mono normal-case text-ink-soft"> {labelSuffix}</span>
        )}
      </span>
      <div className="flex-1 min-w-0 flex items-center">{children}</div>
    </div>
  );
}
