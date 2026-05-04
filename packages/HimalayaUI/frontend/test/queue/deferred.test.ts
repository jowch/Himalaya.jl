import { describe, it, expect, beforeEach } from "vitest";
import {
  makeDeferred,
  clearDeferred,
  getDeferred,
  pendingDeferreds,
} from "../../src/lib/queue/deferred";

describe("makeDeferred / pendingDeferreds", () => {
  beforeEach(() => {
    // Each test starts with a clean registry to avoid cross-test pollution.
    pendingDeferreds.clear();
  });

  it("registers a deferred under its clientOpId and resolves to the supplied value", async () => {
    const d = makeDeferred<{ id: number }>("op-1");
    setTimeout(() => d.resolve({ id: 42 }), 0);
    const result = await d.promise;
    expect(result).toEqual({ id: 42 });
  });

  it("getDeferred finds the deferred by clientOpId", () => {
    const d = makeDeferred<unknown>("op-2");
    expect(getDeferred("op-2")).toBe(d);
  });

  it("getDeferred returns undefined when no deferred is registered", () => {
    expect(getDeferred("missing")).toBeUndefined();
  });

  it("clearDeferred removes the deferred from the registry", () => {
    makeDeferred<unknown>("op-3");
    clearDeferred("op-3");
    expect(getDeferred("op-3")).toBeUndefined();
  });

  it("clearDeferred is idempotent (no throw on missing key)", () => {
    expect(() => clearDeferred("never-registered")).not.toThrow();
  });

  it("AbortSignal-driven rejection propagates", async () => {
    const ctrl = new AbortController();
    const d = makeDeferred<unknown>("op-4");
    ctrl.signal.addEventListener("abort", () =>
      d.reject(new DOMException("aborted", "AbortError")),
    );
    ctrl.abort();
    await expect(d.promise).rejects.toThrow("aborted");
  });

  it("two deferreds with the same key — second overwrites first (caller's responsibility to mint unique IDs)", () => {
    const d1 = makeDeferred<number>("dup");
    const d2 = makeDeferred<number>("dup");
    expect(getDeferred("dup")).toBe(d2);
    expect(getDeferred("dup")).not.toBe(d1);
  });
});
