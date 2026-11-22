import { useEffect, useRef, useState, useCallback } from "react";
import { setToastImpl, type ToastAction, type ToastKind } from "../lib/toast";
import { IconButton } from "./IconButton";

interface ToastItem {
  id: number;
  msg: string;
  kind: ToastKind;
  action?: ToastAction;
}

const DURATIONS: Record<ToastKind, number> = {
  error: 5000,
  info: 3000,
  success: 3000,
  warning: 3000,
};

const KIND_WORD: Record<ToastKind, string> = {
  info: "Info",
  success: "Success",
  warning: "Warning",
  error: "Error",
};

const KIND_GLYPH: Record<ToastKind, string> = {
  info: "i",
  success: "✓",
  warning: "!",
  error: "✕",
};

const KIND_ICON_COLOR: Record<ToastKind, string> = {
  info: "text-accent",
  success: "text-success",
  warning: "text-warning",
  error: "text-error",
};

export function ToastContainer(): JSX.Element {
  const [items, setItems] = useState<ToastItem[]>([]);
  const idRef = useRef(0);
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
    setToastImpl((msg, kind, action) => {
      idRef.current += 1;
      const id = idRef.current;
      const ttl = DURATIONS[kind] ?? 3000;
      setItems((curr) => [...curr, { id, msg, kind, ...(action !== undefined ? { action } : {}) }]);
      const handle = setTimeout(() => {
        timersRef.current.delete(id);
        setItems((curr) => curr.filter((t) => t.id !== id));
      }, ttl);
      timersRef.current.set(id, handle);
    });
    return () => setToastImpl(null);
  }, []);
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
          {t.action ? (
            <button
              type="button"
              className="ml-3 text-xs font-bold text-accent"
              onClick={t.action.onClick}
            >
              {t.action.label}
            </button>
          ) : null}
          <IconButton label="Dismiss" dismiss tone="ghost" onClick={() => dismiss(t.id)} />
        </div>
      ))}
    </div>
  );
}