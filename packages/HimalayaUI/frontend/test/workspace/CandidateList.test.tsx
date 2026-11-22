import { describe, expect, it, vi } from "vitest";
import { render, screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { IndexEntry } from "../../src/api";
import { CandidateList } from "../../src/workspace/ui/CandidateList";
import type { CallRow } from "../../src/workspace/bonnetGroups";

function candidate(id: number, phase = "Lamellar", lattice: number | null = 40 + id): IndexEntry {
  return {
    id, exposure_id: 1, phase, basis: 0.1, score: 1 - id / 100, r_squared: null, lattice_d: lattice,
    status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
  };
}

const one = (candidate: IndexEntry): CallRow => ({ kind: "one", candidate });

describe("CandidateList", () => {
  it("lists every candidate, gives the first nine a number key, and leaves the rest clickable without one", async () => {
    const onToggle = vi.fn();
    const candidates = Array.from({ length: 12 }, (_, i) => candidate(i + 1));
    render(
      <CandidateList rows={candidates.map(one)} assigned={new Set()} disabled={false} formFactor={false}
        onToggle={onToggle} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
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

  it("shows a Bonnet group as one numbered header over numbered members", async () => {
    const onToggle = vi.fn();
    const onToggleGroup = vi.fn();
    const onHover = vi.fn();
    const [a, b, c] = [candidate(1), candidate(2), candidate(3)];
    render(
      <CandidateList rows={[{ kind: "group", members: [a!, b!], peaks: 8 }, one(c!)]} assigned={new Set([a!.id])}
        disabled={false} formFactor={false} onToggle={onToggle} onToggleGroup={onToggleGroup} onFormFactor={() => {}} onHover={onHover} />,
    );
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows.map((r) => r.textContent?.match(/^\d+/)?.[0])).toEqual(["1", "2", "3", "4"]);
    const header = screen.getByTestId("workspace-bonnet-group");
    expect(header).toHaveTextContent("Bonnet group");
    expect(header).toHaveTextContent("8 peaks");
    expect(header).toHaveAttribute("aria-pressed", "mixed");
    await userEvent.hover(header);
    expect(onHover).toHaveBeenLastCalledWith([a!.id, b!.id]);
    await userEvent.click(header);
    expect(onToggleGroup).toHaveBeenCalledWith([a, b]);
    await userEvent.click(rows[2]!);
    expect(onToggle).toHaveBeenCalledWith(b);
  });

  it("marks a group pressed only when every member is in the call, and dashes its border while only some are", () => {
    const [a, b] = [candidate(1), candidate(2)];
    const { rerender } = render(
      <CandidateList rows={[{ kind: "group", members: [a!, b!], peaks: 8 }]} assigned={new Set()} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    expect(screen.getByTestId("workspace-bonnet-group")).toHaveAttribute("aria-pressed", "false");
    expect(screen.getByTestId("workspace-bonnet-group").className).not.toContain("border-dashed");
    rerender(
      <CandidateList rows={[{ kind: "group", members: [a!, b!], peaks: 8 }]} assigned={new Set([a!.id])} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    expect(screen.getByTestId("workspace-bonnet-group")).toHaveAttribute("aria-pressed", "mixed");
    expect(screen.getByTestId("workspace-bonnet-group").className).toContain("border-dashed");
    rerender(
      <CandidateList rows={[{ kind: "group", members: [a!, b!], peaks: 8 }]} assigned={new Set([a!.id, b!.id])} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    expect(screen.getByTestId("workspace-bonnet-group")).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByTestId("workspace-bonnet-group").className).not.toContain("border-dashed");
  });

  it("shows the mean Gaussian curvature in nm⁻² on bicontinuous cubic rows only, on one line", () => {
    render(
      <CandidateList rows={[candidate(1, "Pn3m", 108), candidate(2, "Im3m", 138.5), candidate(3, "Hexagonal", 60), candidate(4, "Ia3d", null)].map(one)}
        assigned={new Set([1])} disabled={false} formFactor={false} onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent(/108\.0 Å\s*κ −0\.0561 nm⁻²/);
    expect(rows[1]).toHaveTextContent(/138\.5 Å\s*κ −0\.0559 nm⁻²/);
    expect(rows[2]).not.toHaveTextContent("κ");
    expect(rows[3]).not.toHaveTextContent("κ");
    expect(rows[0]!.querySelector("[data-testid=workspace-candidate-curvature]")?.parentElement).toHaveClass("whitespace-nowrap");
  });

  it("shows the group's shared curvature on its header", () => {
    render(
      <CandidateList rows={[{ kind: "group", members: [candidate(1, "Im3m", 138.5), candidate(2, "Pn3m", 108)], peaks: 8 }]}
        assigned={new Set()} disabled={false} formFactor={false} onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    expect(screen.getByTestId("workspace-bonnet-curvature")).toHaveTextContent("κ −0.0560 nm⁻²");
  });
});
