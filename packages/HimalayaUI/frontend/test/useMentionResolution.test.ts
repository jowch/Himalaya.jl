import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import { createElement } from "react";
import { makeClient } from "./test-utils";
import * as api from "../src/api";
import { useMentionResolution } from "../src/hooks/useMentionResolution";
import type { MentionToken } from "../src/lib/renderMentions";

const PEAK: api.Peak = {
  id: 42, exposure_id: 1, q: 1.223, intensity: 841.2,
  prominence: 4.2, sharpness: 0.3, source: "auto", excluded: false,
};

function wrapper() {
  const client = makeClient();
  return ({ children }: { children: React.ReactNode }) =>
    createElement(QueryClientProvider, { client }, children);
}

describe("useMentionResolution", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("returns 'loading' initially when entity is not in cache", async () => {
    vi.spyOn(api, "getPeak").mockResolvedValue(PEAK);
    const tokens: MentionToken[] = [{ kind: "mention", type: "peak", id: 42 }];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    expect(result.current.get("peak:42")).toBe("loading");
  });

  it("resolves to entity data after fetch", async () => {
    vi.spyOn(api, "getPeak").mockResolvedValue(PEAK);
    const tokens: MentionToken[] = [{ kind: "mention", type: "peak", id: 42 }];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    await waitFor(() => {
      const entry = result.current.get("peak:42");
      expect(entry).not.toBe("loading");
      expect(entry).not.toBe("dead");
    });
    const entry = result.current.get("peak:42");
    expect(entry).toMatchObject({ type: "peak", data: PEAK });
  });

  it("marks entity as 'dead' on 404", async () => {
    vi.spyOn(api, "getPeak").mockRejectedValue(
      new api.ApiError(404, "not found", null)
    );
    const tokens: MentionToken[] = [{ kind: "mention", type: "peak", id: 999 }];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.get("peak:999")).toBe("dead"));
  });

  it("returns empty map for empty token list", () => {
    const { result } = renderHook(() => useMentionResolution([]), { wrapper: wrapper() });
    expect(result.current.size).toBe(0);
  });

  // ─── Comparison resolution (Phase 10) ────────────────────────────────────
  it("resolves a comparison mention to entity data", async () => {
    const COMP: api.Comparison = {
      id: 7, title: "DOPE titration",
      description: null,
      content_hash: "a1b2c3d4e5f60718",
      created_by: 1,
      created_at: "2026-05-02 10:00:00",
      updated_at: "2026-05-02 10:00:00",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      members: [],
    };
    vi.spyOn(api, "getComparison").mockResolvedValue(COMP);
    const tokens: MentionToken[] = [
      { kind: "mention", type: "comparison", id: 7, hash: "a1b2c3d4" },
    ];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    await waitFor(() => {
      const entry = result.current.get("comparison:7");
      expect(entry).not.toBe("loading");
      expect(entry).not.toBe("dead");
    });
    expect(result.current.get("comparison:7")).toMatchObject({
      type: "comparison", data: COMP,
    });
  });

  it("marks a missing comparison as 'dead' on 404", async () => {
    vi.spyOn(api, "getComparison").mockRejectedValue(
      new api.ApiError(404, "not found", null),
    );
    const tokens: MentionToken[] = [
      { kind: "mention", type: "comparison", id: 999 },
    ];
    const { result } = renderHook(() => useMentionResolution(tokens), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.get("comparison:999")).toBe("dead"));
  });
});
