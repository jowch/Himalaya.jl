import { useEffect, useRef, useState, useCallback } from "react";
import { setToastImpl, type ToastKind } from "../../lib/toast";
import { IconButton } from "./IconButton";

interface ToastItem {
  id: number;
  msg: string;
  kind: ToastKind;
}

const DURATIONS: Record<ToastKind, number> = {
  error: 5000,
  info: 3000,
  success: 3000,
  warning: 3000,
};

// Severity is conveyed by a leading status icon + word, NOT an edge hue.
// The icon's color is decorative; the accessible severity is the aria-label
// (on the icon) plus the visible word.
const KIND_WORD: Record<ToastKind, string> = {
  info: "Info",
  success: "Success",
  warning: "Warning",
  error: "Error",
};

const KIND_GLYPH: Record<ToastKind, string> = {
  info: "i",
  success: "✓", // ✓
  warning: "!",
  error: "✕", // ✕
};

const KIND_ICON_COLOR: Record<ToastKind, string> = {
  info: "text-accent",
  success: "text-success",
  warning: "text-warning",
  error: "text-error",
};

/**
 * ToastContainer — mounts the global toast surface, wires showToast into a
 * local queue on mount, and restores the console.warn fallback on unmount.
 *
 * Render once at the root of the app (App.tsx).
 */
export function ToastContainer(): JSX.Element {
  const [items, setItems] = useState<ToastItem[]>([]);
  const idRef = useRef(0);
  // Tracks the auto-dismiss timer per toast id so we can clear it on manual
  // dismissal (close button) and on unmount, preventing stale setItems calls
  // and harmless-but-noisy "filter an already-absent id" follow-ups.
  const timersRef = useRef<Map<number, ReturnType<typeof setTimeout>>>(new Map());

  const dismiss = useCallback((id: number): void => {
    const handle = timersRef.current.get(id);
    if (handle !== undefined) {
      clearTimeout(handle);
      timersRef.current.delete(id);
    }
    setItems((curr) => curr.filter((t) => t.id !== id));
  }, []);

  useEffect(() => {
    setToastImpl((msg, kind) => {
      idRef.current += 1;
      const id = idRef.current;
      const ttl = DURATIONS[kind] ?? 3000;
      setItems((curr) => [...curr, { id, msg, kind }]);
      const handle = setTimeout(() => {
        timersRef.current.delete(id);
        setItems((curr) => curr.filter((t) => t.id !== id));
      }, ttl);
      timersRef.current.set(id, handle);
    });
    return () => setToastImpl(null);
  }, []);

  // Clear any in-flight auto-dismiss timers when the container unmounts so
  // they don't fire setItems against stale state.
  useEffect(() => {
    const timers = timersRef.current;
    return () => {
      for (const handle of timers.values()) clearTimeout(handle);
      timers.clear();
    };
  }, []);

  return (
    <div
      data-testid="toast-container"
      className="fixed top-4 right-4 z-50 flex flex-col gap-2 pointer-events-none"
    >
      {items.map((t) => (
        <div
          key={t.id}
          data-testid="toast"
          data-toast-kind={t.kind}
          // WCAG 4.1.3: error/warning are interruptive (assertive); info/success
          // are advisory (polite, the default). The severity word stays the
          // human-readable channel; the role/aria-live carry it to assistive tech.
          role={t.kind === "error" || t.kind === "warning" ? "alert" : "status"}
          aria-live={t.kind === "error" || t.kind === "warning" ? "assertive" : "polite"}
          className={
            "anim-toast-in pointer-events-auto flex items-start gap-2 rounded-md border border-hair " +
            "bg-plate text-ink px-3 py-2 shadow-lg text-body min-w-[220px] max-w-[360px]"
          }
        >
          <span
            aria-label={KIND_WORD[t.kind]}
            className={`flex-shrink-0 font-bold leading-tight ${KIND_ICON_COLOR[t.kind]}`}
          >
            {KIND_GLYPH[t.kind]}
          </span>
          <span className="flex-1 break-words">
            <span className="font-semibold">{KIND_WORD[t.kind]}</span>{" "}
            {t.msg}
          </span>
          <IconButton label="Dismiss" dismiss tone="ghost" onClick={() => dismiss(t.id)} />
        </div>
      ))}
    </div>
  );
}
