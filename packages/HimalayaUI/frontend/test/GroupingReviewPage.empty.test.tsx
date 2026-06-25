import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

// Mutable so each test can pose a different empty case. Regression for the
// real-data sweep finding #5: all three empty cases used to collapse to the
// misleading `No samples match ""`.
let LOADS: Load[] = [];
vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return { ...actual, useLoads: () => ({ data: LOADS, isLoading: false }) };
});

function wrap(node: ReactNode) {
  const qc = new QueryClient();
  return render(<QueryClientProvider client={qc}>{node}</QueryClientProvider>);
}

const unflaggedLoad: Load = {
  load_id: 1, load_index: 1, session_id: null, start_time: "10:02", end_time: "10:38",
  frame_count: 0, note: null,
  samples: [{ sample_id: 10, name: "HA85 (S01P15)", slot_index: 15, grouping_source: "auto_position",
    name_source: "auto", merged_into_id: null, flag: null, exposures: [] }],
};

describe("GroupingReviewPage empty states (sweep finding #5)", () => {
  beforeEach(() => { LOADS = []; });

  it("no loads at all → 'No loads yet' + rescan hint, NOT 'No samples match'", () => {
    LOADS = [];
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    expect(screen.getByText("No loads yet")).toBeInTheDocument();
    expect(screen.getByText(/rescan/i)).toBeInTheDocument();
    expect(screen.queryByText(/No samples match/)).toBeNull();
  });

  it("loads present but nothing flagged → 'Nothing needs review', NOT 'No samples match'", () => {
    LOADS = [unflaggedLoad];
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    expect(screen.getByText("Nothing needs review")).toBeInTheDocument();
    expect(screen.queryByText(/No samples match/)).toBeNull();
  });

  it("active search with no match → keeps the literal 'No samples match \"q\"'", () => {
    LOADS = [unflaggedLoad];
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "ZZZ-nomatch" } });
    expect(screen.getByText('No samples match "ZZZ-nomatch"')).toBeInTheDocument();
  });

  it("never renders the empty-quote bug `No samples match \"\"`", () => {
    LOADS = [];
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    expect(screen.queryByText('No samples match ""')).toBeNull();
  });
});
