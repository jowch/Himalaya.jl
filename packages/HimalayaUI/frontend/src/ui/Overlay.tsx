import { useEffect, useId, useLayoutEffect, useRef, useState, type ReactNode } from "react";
import { FOCUSABLE, useFocusTrap } from "../hooks/useFocusTrap";
import { IconButton } from "./IconButton";

interface Props {
  label: string;
  onClose: () => void;
  children: ReactNode;
  testId?: string;
  closeKeys?: readonly string[];
  onKey?: (key: string, repeat: boolean) => boolean;
  heading?: boolean;
}

function useBodyScroll(body: React.RefObject<HTMLDivElement | null>): { more: boolean; bare: boolean } {
  const [state, setState] = useState({ more: false, bare: false });
  const measure = (): void => {
    const el = body.current;
    if (!el) return;
    const more = el.scrollTop + el.clientHeight < el.scrollHeight - 1;
    const bare = el.querySelector(FOCUSABLE) === null;
    setState((s) => (s.more === more && s.bare === bare ? s : { more, bare }));
  };
  useLayoutEffect(measure);
  useEffect(() => {
    const el = body.current;
    window.addEventListener("resize", measure);
    el?.addEventListener("scroll", measure, { passive: true });
    return () => {
      window.removeEventListener("resize", measure);
      el?.removeEventListener("scroll", measure);
    };
  });
  return state;
}

export function Overlay({ label, onClose, children, testId, closeKeys = ["Escape"], onKey, heading = false }: Props): JSX.Element {
  const panel = useRef<HTMLDivElement | null>(null);
  const body = useRef<HTMLDivElement | null>(null);
  const headingId = useId();
  useFocusTrap(panel, true);
  const { more, bare } = useBodyScroll(body);
  useEffect(() => {
    if (body.current && body.current.querySelector(FOCUSABLE) === null) body.current.focus();
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
        aria-label={heading ? undefined : label}
        aria-labelledby={heading ? headingId : undefined}
        data-keys-trap
        data-testid={testId}
        tabIndex={-1}
        className="flex max-h-[calc(88vh-1.5rem)] w-full max-w-xl flex-col overflow-hidden rounded-sm border border-hair bg-plate text-ink shadow-lg outline-none"
        onKeyDown={(e) => {
          if (closeKeys.includes(e.key)) {
            e.preventDefault();
            if (!e.repeat) onClose();
          } else if (!e.metaKey && !e.ctrlKey && !e.altKey && onKey?.(e.key, e.repeat)) {
            e.preventDefault();
          }
        }}
      >
        {heading && (
          <header data-testid="overlay-header" className="flex shrink-0 items-center justify-between gap-4 border-b border-hair py-2 pl-5 pr-3">
            <h2 id={headingId} className="text-label text-ink-soft">{label}</h2>
            <IconButton label="Close" dismiss onClick={onClose} />
          </header>
        )}
        <div
          ref={body}
          data-testid="overlay-body"
          data-more={more || undefined}
          tabIndex={bare ? 0 : -1}
          className="min-h-0 overflow-y-auto overscroll-contain outline-none"
        >
          {children}
          {more && <div aria-hidden="true" data-testid="overlay-more" className="pointer-events-none sticky bottom-0 -mt-8 h-8 bg-linear-to-t from-plate to-transparent" />}
        </div>
      </div>
    </div>
  );
}
