import { useEffect, useRef, useState, useCallback } from "react";
import { setToastImpl, type ToastKind } from "../../lib/toast";

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

const KIND_CLASS: Record<ToastKind, string> = {
  // Subtle elevated surface with a left-edge accent stripe matching the kind.
  info:    "bg-bg-elevated border-accent text-fg",
  success: "bg-bg-elevated border-success text-fg",
  warning: "bg-bg-elevated border-warning text-fg",
  error:   "bg-bg-elevated border-error text-fg",
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

  const dismiss = useCallback((id: number): void => {
    setItems((curr) => curr.filter((t) => t.id !== id));
  }, []);

  useEffect(() => {
    setToastImpl((msg, kind) => {
      idRef.current += 1;
      const id = idRef.current;
      const ttl = DURATIONS[kind] ?? 3000;
      setItems((curr) => [...curr, { id, msg, kind }]);
      window.setTimeout(() => {
        setItems((curr) => curr.filter((t) => t.id !== id));
      }, ttl);
    });
    return () => setToastImpl(null);
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
          role="status"
          className={
            "pointer-events-auto flex items-start gap-2 rounded-md border-l-4 " +
            "px-3 py-2 shadow-lg text-body min-w-[220px] max-w-[360px] " +
            KIND_CLASS[t.kind]
          }
        >
          <span className="flex-1 break-words">{t.msg}</span>
          <button
            type="button"
            aria-label="Dismiss"
            onClick={() => dismiss(t.id)}
            className="text-fg-muted hover:text-fg leading-none px-1"
          >
            ×
          </button>
        </div>
      ))}
    </div>
  );
}
