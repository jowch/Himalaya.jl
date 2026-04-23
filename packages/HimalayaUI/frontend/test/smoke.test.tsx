import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { App } from "../src/App";
import { renderWithProviders } from "./test-utils";

function mockFetch(map: Record<string, unknown>): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    const key = Object.keys(map).find((k) => url.endsWith(k));
    if (!key) return new Response("not found", { status: 404 });
    return new Response(JSON.stringify(map[key]), {
      status: 200, headers: { "Content-Type": "application/json" },
    });
  });
}

describe("App smoke", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    localStorage.clear();
    mockFetch({
      "/api/users": [],
      "/api/experiments/1": {
        id: 1, name: "demo", path: "/x", data_dir: "/x/data",
        analysis_dir: "/x/analysis", manifest_path: null,
        created_at: "2026-04-22T00:00:00Z",
      },
      "/api/experiments/1/samples": [
        { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
      ],
    });
  });

  it("renders breadcrumb and sample list after boot", async () => {
    renderWithProviders(<App />);
    await waitFor(() => expect(screen.getByText(/demo/)).toBeInTheDocument());
    expect(screen.getByText("A1")).toBeInTheDocument();
  });
});
