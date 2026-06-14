import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { downloadBlob } from "../../src/lib/figure-export/download";

// The download anchor's `download` attribute is the figure's filename. If the
// object URL is revoked SYNCHRONOUSLY after click() it can race the browser's
// async download and make it fall back to the blob UUID (and skip the downloads
// folder). These tests pin: the anchor carries the right name, click() fires,
// and revocation is DEFERRED, not synchronous.
describe("downloadBlob", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    // jsdom lacks createObjectURL/revokeObjectURL — stub them.
    Object.defineProperty(URL, "createObjectURL", {
      configurable: true,
      value: vi.fn(() => "blob:mock-uuid"),
    });
    Object.defineProperty(URL, "revokeObjectURL", {
      configurable: true,
      value: vi.fn(),
    });
  });
  afterEach(() => {
    vi.useRealTimers();
    vi.restoreAllMocks();
  });

  it("sets the download filename on the anchor and clicks it", () => {
    let captured: HTMLAnchorElement | undefined;
    const clickSpy = vi
      .spyOn(HTMLAnchorElement.prototype, "click")
      .mockImplementation(function (this: HTMLAnchorElement) {
        captured = this;
      });
    downloadBlob(new Blob(["x"], { type: "image/svg+xml" }), "himalaya-series-fig-2026-05-08.svg");
    expect(clickSpy).toHaveBeenCalledOnce();
    expect(captured?.getAttribute("download")).toBe("himalaya-series-fig-2026-05-08.svg");
    expect(captured?.getAttribute("href")).toBe("blob:mock-uuid");
  });

  it("does NOT revoke the object URL synchronously (avoids the UUID-name race)", () => {
    vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(() => {});
    downloadBlob(new Blob(["x"]), "fig.png");
    // Synchronously after the call, the URL is still alive.
    expect(URL.revokeObjectURL).not.toHaveBeenCalled();
    // It is cleaned up on a later tick.
    vi.runAllTimers();
    expect(URL.revokeObjectURL).toHaveBeenCalledWith("blob:mock-uuid");
  });
});
