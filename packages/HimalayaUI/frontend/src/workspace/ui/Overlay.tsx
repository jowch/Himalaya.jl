import { useEffect, useRef, type ReactNode } from "react";

interface Props {
  label: string;
  onClose: () => void;
  children: ReactNode;
  testId?: string;
}

export function Overlay({ label, onClose, children, testId }: Props): JSX.Element {
  const panel = useRef<HTMLDivElement | null>(null);
  useEffect(() => {
    const prev = document.activeElement instanceof HTMLElement ? document.activeElement : null;
    const target = panel.current?.querySelector<HTMLElement>("input,button,[tabindex]") ?? panel.current;
    target?.focus();
    return () => prev?.focus();
  }, []);
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
          if (e.key === "Escape") {
            e.preventDefault();
            onClose();
          } else if (e.key === "Tab") {
            const items = Array.from(panel.current?.querySelectorAll<HTMLElement>("input,button:not([tabindex='-1'])") ?? []);
            if (items.length === 0) {
              e.preventDefault();
              return;
            }
            const i = items.indexOf(document.activeElement as HTMLElement);
            const next = items[(i + (e.shiftKey ? -1 : 1) + items.length) % items.length]!;
            e.preventDefault();
            next.focus();
          }
        }}
      >
        {children}
      </div>
    </div>
  );
}
