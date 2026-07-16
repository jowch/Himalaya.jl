// Minimal toast stub — M2.2 replaces with a real Toast UI mounted at App.tsx.
// Until then this surfaces validation errors via console.warn so they are
// visible in dev. Production behavior changes when M2.2 ships.

export type ToastKind = "info" | "success" | "warning" | "error";

/** Optional inline action (e.g. Undo) rendered inside a toast. */
export interface ToastAction {
  label: string;
  onClick: () => void;
}

let activeImpl: (msg: string, kind: ToastKind, action?: ToastAction) => void = (msg, kind) => {
  // eslint-disable-next-line no-console
  console.warn(`[toast:${kind}] ${msg}`);
};

export function showToast(msg: string, kind: ToastKind = "info", action?: ToastAction): void {
  // Only forward `action` when present so the impl is still called with the
  // 2-arg (msg, kind) shape for the common no-action case — keeps existing
  // toast-impl spies (toHaveBeenCalledWith(msg, kind)) exact, no trailing undefined.
  if (action !== undefined) activeImpl(msg, kind, action);
  else activeImpl(msg, kind);
}

/**
 * Replace the toast implementation. M2.2's ToastProvider calls this on mount.
 * Subsequent calls override; passing `null` restores the console.warn fallback.
 */
export function setToastImpl(
  impl: ((msg: string, kind: ToastKind, action?: ToastAction) => void) | null,
): void {
  if (impl === null) {
    activeImpl = (msg, kind) => {
      // eslint-disable-next-line no-console
      console.warn(`[toast:${kind}] ${msg}`);
    };
  } else {
    activeImpl = impl;
  }
}
