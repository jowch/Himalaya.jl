import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";
import type { SeriesSummary, Series } from "../src/api";

// A fixture whose keys are exactly what _series_listing_rows emits
// (packages/HimalayaUI/src/series.jl). If SeriesSummary drifts from the
// backend projection, this assignment stops type-checking.
const LISTING_ROW: SeriesSummary = {
  id: 1,
  title: "Lipid dose response",
  description: null,
  content_hash: "abc123",
  created_by: 2,
  created_at: "2026-05-01T10:00:00",
  updated_at: "2026-05-02T10:00:00",
  forked_from_id: null,
  forked_at_hash: null,
  view_grouping_mode: null,
  view_show_peak_ticks: null,
  view_show_peak_labels: null,
  last_event_at: "2026-05-02 10:00:00",
  author_username: "jc",
  member_count: 3,
  member_phases: ["Pn3m", "Lamellar"],
  member_phase_count: 2,
  has_stale_members: false,
};

describe("series read fetchers", () => {
  beforeEach(() => {
    vi.stubGlobal("fetch", vi.fn());
  });
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  function mockJson(body: unknown) {
    (globalThis.fetch as ReturnType<typeof vi.fn>).mockResolvedValue({
      ok: true,
      status: 200,
      headers: { get: () => "application/json" },
      json: async () => body,
      text: async () => JSON.stringify(body),
    } as unknown as Response);
  }

  it("listSeries GETs /api/series and returns SeriesSummary[]", async () => {
    mockJson([LISTING_ROW]);
    const rows = await api.listSeries();
    const init = (globalThis.fetch as ReturnType<typeof vi.fn>).mock.calls[0][1];
    expect((globalThis.fetch as ReturnType<typeof vi.fn>).mock.calls[0][0]).toContain(
      "/api/series",
    );
    expect(init?.method ?? "GET").toBe("GET");
    expect(rows[0].member_count).toBe(3);
  });

  it("getSeries GETs /api/series/:id and returns the full Series shape", async () => {
    const full: Series = {
      ...LISTING_ROW,
      forked_from_title: null,
      ordering_variable: null,
      order_rule: "ascending",
      state: "draft",
      members: [],
      samples: [],
    };
    mockJson(full);
    const s = await api.getSeries(7);
    expect((globalThis.fetch as ReturnType<typeof vi.fn>).mock.calls[0][0]).toContain(
      "/api/series/7",
    );
    expect(s.members).toEqual([]);
  });

  it("forksOfSeries GETs /api/series/:id/forks", async () => {
    mockJson([]);
    await api.forksOfSeries(7);
    expect((globalThis.fetch as ReturnType<typeof vi.fn>).mock.calls[0][0]).toContain(
      "/api/series/7/forks",
    );
  });
});
