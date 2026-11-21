import { describe, it, expect, vi, afterEach } from "vitest";
import { request, NetworkError } from "../src/api";

describe("api.request failure classes", () => {
  afterEach(() => vi.unstubAllGlobals());

  it("turns a refused connection into a NetworkError", async () => {
    vi.stubGlobal("fetch", vi.fn().mockRejectedValue(new TypeError("Failed to fetch")));
    await expect(request("PUT", "/api/x", {})).rejects.toBeInstanceOf(NetworkError);
  });

  it("leaves a header the browser refuses as a plain TypeError, since retrying cannot fix it", async () => {
    const fetchSpy = vi.fn();
    vi.stubGlobal("fetch", fetchSpy);
    const err = await request("PUT", "/api/x", {}, { username: "名前" }).catch((e: unknown) => e);
    expect(err).toBeInstanceOf(TypeError);
    expect(err).not.toBeInstanceOf(NetworkError);
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("passes an abort through unchanged", async () => {
    vi.stubGlobal("fetch", vi.fn().mockRejectedValue(new DOMException("aborted", "AbortError")));
    const err = await request("GET", "/api/x").catch((e: unknown) => e);
    expect(err).toBeInstanceOf(DOMException);
  });
});
