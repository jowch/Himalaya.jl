import { describe, expect, it, vi } from "vitest";
import { render, screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { IndexEntry } from "../../src/api";
import { CandidateList } from "../../src/workspace/ui/CandidateList";

function candidate(id: number): IndexEntry {
  return {
    id, exposure_id: 1, phase: "Lamellar", basis: 0.1, score: 1 - id / 100, r_squared: null, lattice_d: 40 + id,
    ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
  };
}

describe("CandidateList", () => {
  it("lists every candidate, gives the first nine a number key, and leaves the rest clickable without one", async () => {
    const onToggle = vi.fn();
    const candidates = Array.from({ length: 12 }, (_, i) => candidate(i + 1));
    render(
      <CandidateList candidates={candidates} assigned={new Set()} disabled={false} formFactor={false}
        onToggle={onToggle} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows).toHaveLength(12);
    const key = (row: HTMLElement) => {
      const copy = row.cloneNode(true) as HTMLElement;
      copy.querySelectorAll("[aria-hidden=true]").forEach((n) => n.remove());
      return copy.textContent?.match(/^\d+/)?.[0] ?? null;
    };
    expect(rows.map(key)).toEqual(["1", "2", "3", "4", "5", "6", "7", "8", "9", null, null, null]);
    await userEvent.click(rows[11]!);
    expect(onToggle).toHaveBeenCalledWith(candidates[11]);
  });
});
