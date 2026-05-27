import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";
import type { Series } from "../src/api";

// I3.5b — commitSeriesPlate maps a 409 (content_hash drift) to the typed
// `ConflictError`, mirroring `saveComparison`. This is the ONLY series fetcher
// that throws ConflictError: recipe-save (PATCH /api/series/:id) never reads
// `expected_content_hash` and never 409s, so `saveSeries` is left untouched.

function serverSeries(id: number): Series {
  return {
    id, title: "Server title", description: null, content_hash: "sha256:server",
    created_by: 1, created_at: null, updated_at: "2026-05-02",
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, ordering_variable: null, order_rule: "manual",
    state: "committed", members: [], samples: [],
  };
}

describe("commitSeriesPlate conflict mapping", () => {
  beforeEach(() => {
    vi.stubGlobal("fetch", vi.fn());
  });
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  function mock409(body: unknown) {
    (globalThis.fetch as ReturnType<typeof vi.fn>).mockResolvedValue({
      ok: false,
      status: 409,
      headers: { get: () => "application/json" },
      json: async () => body,
      text: async () => JSON.stringify(body),
    } as unknown as Response);
  }

  function mockOk(body: unknown) {
    (globalThis.fetch as ReturnType<typeof vi.fn>).mockResolvedValue({
      ok: true,
      status: 200,
      headers: { get: () => "application/json" },
      json: async () => body,
      text: async () => JSON.stringify(body),
    } as unknown as Response);
  }

  it("throws ConflictError carrying current_hash + current_state (a Series) on 409", async () => {
    const current = serverSeries(5);
    mock409({ current_hash: "sha256:server", current_state: current });
    await expect(
      api.commitSeriesPlate(5, { members: [] }),
    ).rejects.toBeInstanceOf(api.ConflictError);
    try {
      await api.commitSeriesPlate(5, { members: [] });
    } catch (err) {
      const ce = err as api.ConflictError;
      expect(ce.status).toBe(409);
      expect(ce.current_hash).toBe("sha256:server");
      // current_state is now widened to Comparison | Series | null; the
      // commit conflict populates it with a Series.
      expect((ce.current_state as Series).id).toBe(5);
      expect((ce.current_state as Series).state).toBe("committed");
    }
  });

  it("POSTs /api/series/:id/commit and returns the Series on 200", async () => {
    const committed = serverSeries(5);
    mockOk(committed);
    const out = await api.commitSeriesPlate(5, { members: [] });
    const url = (globalThis.fetch as ReturnType<typeof vi.fn>).mock.calls[0][0];
    const init = (globalThis.fetch as ReturnType<typeof vi.fn>).mock.calls[0][1];
    expect(String(url)).toContain("/api/series/5/commit");
    expect(init?.method).toBe("POST");
    expect(out.id).toBe(5);
  });
});
