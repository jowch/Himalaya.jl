export type ToastKind = "info" | "success" | "warning" | "error";

export interface ToastAction {
  label: string;
  onClick: () => void;
}

type ToastImpl = (msg: string, kind: ToastKind, action?: ToastAction) => void;

let activeImpl: ToastImpl | null = null;
let pending: [string, ToastKind, ToastAction | undefined][] = [];

export function showToast(msg: string, kind: ToastKind = "info", action?: ToastAction): void {
  if (activeImpl === null) {
    if (pending.length === 0) queueMicrotask(dropPending);
    pending.push([msg, kind, action]);
    return;
  }
  if (action !== undefined) activeImpl(msg, kind, action);
  else activeImpl(msg, kind);
}

function dropPending(): void {
  for (const [msg, kind] of pending) console.warn(`[toast:${kind}] ${msg}`);
  pending = [];
}

export function setToastImpl(impl: ToastImpl | null): void {
  activeImpl = impl;
  if (impl === null) return;
  const queued = pending;
  pending = [];
  for (const [msg, kind, action] of queued) showToast(msg, kind, action);
}
