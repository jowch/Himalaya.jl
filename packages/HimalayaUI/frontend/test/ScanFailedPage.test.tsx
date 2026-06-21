// test/ScanFailedPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { ScanFailedPage } from "../src/print/pages/ScanFailedPage";
import type { ManifestUnmatched } from "../src/api";

// Module-level navigate mock — canonical pattern.
const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderScanFailed({
  unmatched,
  parsedCount,
}: {
  unmatched: ManifestUnmatched[];
  parsedCount: number;
}) {
  render(
    <MemoryRouter initialEntries={["/experiments/7/corpus"]}>
      <Routes>
        <Route
          path="/experiments/:id/corpus"
          element={
            <ScanFailedPage
              experimentId={7}
              unmatched={unmatched}
              parsedCount={parsedCount}
            />
          }
        />
      </Routes>
    </MemoryRouter>,
  );
}

describe("ScanFailedPage", () => {
  beforeEach(() => {
    navigate.mockClear();
  });
  afterEach(() => vi.restoreAllMocks());

  it("groups misses, offers per-type test + ingest-N confirm", () => {
    renderScanFailed({
      unmatched: [
        { file: "a", miss: "metadata" },
        { file: "b", miss: "integration" },
      ],
      parsedCount: 5,
    });

    // Open Configuration primary action
    expect(screen.getByRole("button", { name: /open configuration/i })).toBeInTheDocument();

    // Unmatched stem shown
    expect(screen.getByText("a")).toBeInTheDocument();

    // Miss-type label shown
    expect(screen.getByText(/metadata/i)).toBeInTheDocument();

    // One textbox per affected type (2 types: metadata + integration)
    expect(screen.getAllByRole("textbox").length).toBe(2);

    // Two-stage confirm: first click shows Confirm
    fireEvent.click(screen.getByRole("button", { name: /ingest 5/i }));
    expect(screen.getByText(/confirm/i)).toBeInTheDocument();
  });

  it("renders only one textbox when only one miss type is present", () => {
    renderScanFailed({
      unmatched: [
        { file: "x", miss: "metadata" },
        { file: "y", miss: "metadata" },
      ],
      parsedCount: 3,
    });
    expect(screen.getAllByRole("textbox").length).toBe(1);
  });

  it("navigates to config tab when Open Configuration is clicked", () => {
    renderScanFailed({ unmatched: [], parsedCount: 0 });
    fireEvent.click(screen.getByRole("button", { name: /open configuration/i }));
    expect(navigate).toHaveBeenCalledWith("/experiments/7/config");
  });
});
