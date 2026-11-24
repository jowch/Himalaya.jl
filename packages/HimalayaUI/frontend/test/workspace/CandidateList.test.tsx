import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { IndexEntry } from "../../src/api";
import { CandidateList } from "../../src/ui/CandidateList";
import { phaseColor } from "../../src/phases";
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

  it("keeps a blank keycap on a choice past the ninth, so its tick still shows", () => {
    const candidates = Array.from({ length: 11 }, (_, i) => candidate(i + 1));
    render(
      <CandidateList rows={candidates.map(one)} assigned={new Set([11])} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[10]!.querySelector("[data-blank-key] kbd")).toHaveAttribute("data-pressed", "true");
    expect(rows[9]!.querySelector("[data-blank-key] kbd")).not.toHaveAttribute("data-pressed");
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

  it("marks a group pressed only when every member is in the call, and outlines its keycap while only some are", () => {
    const [a, b] = [candidate(1), candidate(2)];
    const props = { rows: [{ kind: "group" as const, members: [a!, b!], peaks: 8 }], disabled: false, formFactor: false,
      onToggle: () => {}, onToggleGroup: () => {}, onFormFactor: () => {}, onHover: () => {} };
    const header = () => screen.getByTestId("workspace-bonnet-group");
    const keycap = () => header().querySelector("kbd")!;
    const { rerender } = render(<CandidateList {...props} assigned={new Set()} />);
    expect(header()).toHaveAttribute("aria-pressed", "false");
    expect(keycap()).not.toHaveAttribute("data-pressed");
    rerender(<CandidateList {...props} assigned={new Set([a!.id])} />);
    expect(header()).toHaveAttribute("aria-pressed", "mixed");
    expect(keycap()).toHaveAttribute("data-pressed", "mixed");
    expect(keycap().className).toMatch(/(^| )border-accent( |$)/);
    expect(keycap().className).not.toContain("bg-accent-fill");
    rerender(<CandidateList {...props} assigned={new Set([a!.id, b!.id])} />);
    expect(header()).toHaveAttribute("aria-pressed", "true");
    expect(keycap()).toHaveAttribute("data-pressed", "true");
  });

  it("shows a ticked choice by filling its keycap, leaving the row untinted so selection never reads as keyboard focus", () => {
    const [a, b] = [candidate(1), candidate(2)];
    render(
      <CandidateList rows={[a!, b!].map(one)} assigned={new Set([a!.id])} disabled={false} formFactor
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const [on, off] = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    const ff = screen.getByTestId("workspace-form-factor");
    for (const el of [on!, ff, off!]) {
      expect(el.className).not.toContain("bg-row-selected");
      expect(el.className).not.toContain("shadow-");
      expect(el.className).not.toMatch(/border-ink/);
    }
    for (const el of [on!, ff]) {
      expect(el.querySelector("kbd")).toHaveAttribute("data-pressed", "true");
      expect(el.querySelector("kbd")!.className).toMatch(/(^| )bg-accent-fill( |$)/);
      expect(el.querySelector("kbd")!.className).toMatch(/(^| )text-on-accent( |$)/);
    }
    expect(off!.querySelector("kbd")).not.toHaveAttribute("data-pressed");
  });

  it("names its phase as plain text in the phase colour, so the chosen keycap is the only fill in the row", () => {
    render(
      <CandidateList rows={[one(candidate(1, "Pn3m"))]} assigned={new Set([1])} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const name = within(screen.getByTestId("workspace-candidates")).getByTestId("phase-chip");
    expect(name).toHaveAttribute("data-variant", "text");
    expect(name.style.background).toBe("");
    expect(name.style.color).toBe(phaseColor("Pn3m"));
  });

  it("keeps a ticked, a partly ticked and a hovered choice apart by the keycap alone, since a row tint for ticks reads like hover", () => {
    const [a, b, c] = [candidate(1), candidate(2), candidate(3)];
    render(
      <CandidateList rows={[{ kind: "group", members: [a!, b!], peaks: 8 }, one(c!)]} assigned={new Set([a!.id])} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const [group, ticked, , unticked] = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    const rowLook = (el: HTMLElement) => el.className.split(" ").filter((c) => c.startsWith("bg-") || c.startsWith("hover:bg-")).sort();
    expect(rowLook(ticked!)).toEqual(rowLook(unticked!));
    expect(rowLook(group!)).toEqual(rowLook(unticked!));
    const keys = [group, ticked, unticked].map((el) => el!.querySelector("kbd")!.getAttribute("data-pressed"));
    expect(keys).toEqual(["mixed", "true", null]);
  });

  it("reports hover on pointer movement, so a row that appears under a resting pointer previews as soon as the pointer moves", () => {
    const onHover = vi.fn();
    const [a, b] = [candidate(1), candidate(2)];
    render(
      <CandidateList rows={[{ kind: "group", members: [a!, b!], peaks: 8 }]} assigned={new Set()} disabled={false} formFactor={false}
        onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={onHover} />,
    );
    const [group, first] = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    fireEvent.mouseMove(first!);
    expect(onHover).toHaveBeenLastCalledWith([a!.id]);
    fireEvent.mouseMove(group!);
    expect(onHover).toHaveBeenLastCalledWith([a!.id, b!.id]);
  });

  it("leaves the curvature off single candidate rows, cubic or not", () => {
    render(
      <CandidateList rows={[candidate(1, "Pn3m", 108), candidate(2, "Im3m", 138.5), candidate(3, "Hexagonal", 60)].map(one)}
        assigned={new Set([1])} disabled={false} formFactor={false} onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent("108.0 Å");
    for (const row of rows) expect(row).not.toHaveTextContent("κ");
  });

  it("shows the group's shared curvature on its header", () => {
    render(
      <CandidateList rows={[{ kind: "group", members: [candidate(1, "Im3m", 138.5), candidate(2, "Pn3m", 108)], peaks: 8 }]}
        assigned={new Set()} disabled={false} formFactor={false} onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />,
    );
    expect(screen.getByTestId("workspace-bonnet-curvature")).toHaveTextContent("κ −0.0560 nm⁻²");
    for (const member of within(screen.getByTestId("workspace-candidates")).getAllByRole("button").slice(1)) expect(member).not.toHaveTextContent("κ");
  });
});
