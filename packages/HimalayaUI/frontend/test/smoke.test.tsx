import { render, screen } from "@testing-library/react";
import { describe, it, expect, beforeEach, vi } from "vitest";
import { App } from "../src/App";

describe("App smoke", () => {
  beforeEach(() => {
    // Stub fetch so boot doesn't throw unhandled rejections in test logs
    global.fetch = vi.fn(async () =>
      new Response(JSON.stringify([]), { status: 200, headers: { "Content-Type": "application/json" } })
    ) as unknown as typeof fetch;
    // Clean persisted state between tests
    localStorage.clear();
  });

  it("renders the Himalaya navbar logo", () => {
    render(<App />);
    expect(screen.getByText("Himalaya")).toBeInTheDocument();
  });
});
