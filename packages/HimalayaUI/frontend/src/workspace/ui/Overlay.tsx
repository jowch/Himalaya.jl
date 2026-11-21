import { useRef, type ReactNode } from "react";
import { useFocusTrap } from "../../hooks/useFocusTrap";

interface Props {
  label: string;
  onClose: () => void;
  children: ReactNode;
  testId?: string;
  closeKeys?: readonly string[];
  onKey?: (key: string, repeat: boolean) => boolean;
}

export function Overlay({ label, onClose, children, testId, closeKeys = ["Escape"], onKey }: Props): JSX.Element {
  const panel = useRef<HTMLDivElement | null>(null);
  useFocusTrap(panel, true);
  return (
    <div
      className="fixed inset-0 z-40 flex items-start justify-center bg-scrim px-4 pt-[12vh]"
      onMouseDown={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div
        ref={panel}
        role="dialog"
        aria-modal="true"
        aria-label={label}
        data-keys-trap
        data-testid={testId}
        tabIndex={-1}
        className="w-full max-w-xl overflow-hidden rounded-sm border border-hair bg-plate text-ink shadow-lg outline-none"
        onKeyDown={(e) => {
          if (closeKeys.includes(e.key)) {
            e.preventDefault();
            if (!e.repeat) onClose();
          } else if (!e.metaKey && !e.ctrlKey && !e.altKey && onKey?.(e.key, e.repeat)) {
            e.preventDefault();
          }
        }}
      >
        {children}
      </div>
    </div>
  );
}
