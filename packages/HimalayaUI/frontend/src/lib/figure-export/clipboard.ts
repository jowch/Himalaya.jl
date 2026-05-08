// clipboard.ts — write a PNG blob to the system clipboard.

/**
 * Pre-flight feature check. The Copy button uses this to decide whether to
 * disable itself (insecure origin / unsupported browser) before the user
 * even clicks. Spec §Browser compatibility / Disabled states.
 */
export function canCopyPngToClipboard(): boolean {
  if (typeof navigator === "undefined") return false;
  if (!navigator.clipboard) return false;
  if (typeof (globalThis as { ClipboardItem?: unknown }).ClipboardItem === "undefined") {
    return false;
  }
  return true;
}

/**
 * Write a PNG blob to the clipboard. Caller is expected to have gated
 * on `canCopyPngToClipboard()`; we don't re-check here so a runtime
 * rejection (user denial, third failure mode) propagates as a real error
 * the caller can toast.
 */
export async function copyPngToClipboard(blob: Blob): Promise<void> {
  // ClipboardItem is global at runtime when canCopyPngToClipboard() is true.
  // TS doesn't pick up DOM-lib types in the test env; cast at call time.
  const Ctor = (globalThis as { ClipboardItem: new (data: Record<string, Blob>) => unknown })
    .ClipboardItem;
  const item = new Ctor({ "image/png": blob }) as unknown;
  await navigator.clipboard.write([item as ClipboardItem]);
}
