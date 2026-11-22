import { describe, expect, it, vi } from "vitest";
import { render, screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { IndexEntry } from "../../src/api";
import { CandidateList } from "../../src/workspace/ui/CandidateList";

function candidate(id: number, phase = "Lamellar", lattice: number | null = 40 + id): IndexEntry {
  return {
    id, exposure_id: 1, phase, basis: 0.1, score: 1 - id / 100, r_squared: null, lattice_d: lattice,
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

  it("shows the mean Gaussian curvature in nm⁻² on bicontinuous cubic rows only, on one line", () => {
    render(
      <CandidateList candidates={[candidate(1, "Pn3m", 108), candidate(2, "Im3m", 138.5), candidate(3, "Hexagonal", 60), candidate(4, "Ia3d", null)]}
        assigned={new Set([1])} disabled={false} formFactor={false} onToggle={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent(/108\.0 Å\s*κ −0\.0561 nm⁻²/);
    expect(rows[1]).toHaveTextContent(/138\.5 Å\s*κ −0\.0559 nm⁻²/);
    expect(rows[2]).not.toHaveTextContent("κ");
    expect(rows[3]).not.toHaveTextContent("κ");
    expect(rows[0]!.querySelector("[data-testid=workspace-candidate-curvature]")?.parentElement).toHaveClass("whitespace-nowrap");
  });
});
