import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import {
  canCopyPngToClipboard,
  copyPngToClipboard,
} from "../../src/lib/figure-export/clipboard";

const origNav = globalThis.navigator;
const origClipboardItem = (globalThis as { ClipboardItem?: unknown }).ClipboardItem;

afterEach(() => {
  vi.unstubAllGlobals();
  // Restore in case stubs leak.
  Object.defineProperty(globalThis, "navigator", {
    value: origNav,
    configurable: true,
  });
  if (origClipboardItem !== undefined) {
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem = origClipboardItem;
  }
});

describe("canCopyPngToClipboard", () => {
  it("returns false when navigator.clipboard is missing", () => {
    vi.stubGlobal("navigator", {});
    expect(canCopyPngToClipboard()).toBe(false);
  });
  it("returns false when ClipboardItem is missing", () => {
    vi.stubGlobal("navigator", { clipboard: { write: vi.fn() } });
    delete (globalThis as { ClipboardItem?: unknown }).ClipboardItem;
    expect(canCopyPngToClipboard()).toBe(false);
  });
  it("returns true when both are present", () => {
    vi.stubGlobal("navigator", { clipboard: { write: vi.fn() } });
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem = function () { /* stub */ };
    expect(canCopyPngToClipboard()).toBe(true);
  });
});

describe("copyPngToClipboard", () => {
  it("writes a ClipboardItem with image/png MIME", async () => {
    const writeSpy = vi.fn().mockResolvedValue(undefined);
    vi.stubGlobal("navigator", { clipboard: { write: writeSpy } });
    const ciCalls: Array<Record<string, Blob>> = [];
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem =
      function (this: unknown, items: Record<string, Blob>) {
        ciCalls.push(items);
      };
    const blob = new Blob([new Uint8Array([1, 2, 3])], { type: "image/png" });

    await copyPngToClipboard(blob);

    expect(writeSpy).toHaveBeenCalledTimes(1);
    expect(ciCalls.length).toBe(1);
    expect(ciCalls[0]!["image/png"]).toBe(blob);
  });

  it("rejects when navigator.clipboard.write rejects", async () => {
    const writeSpy = vi.fn().mockRejectedValue(new Error("denied"));
    vi.stubGlobal("navigator", { clipboard: { write: writeSpy } });
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem =
      function (this: unknown, _items: Record<string, Blob>) { /* stub */ };
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });

    await expect(copyPngToClipboard(blob)).rejects.toThrow("denied");
  });
});
