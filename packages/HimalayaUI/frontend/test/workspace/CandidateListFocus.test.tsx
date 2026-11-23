import { describe, expect, it } from "vitest";
import { act, fireEvent, render, screen } from "@testing-library/react";
import type { IndexEntry } from "../../src/api";
import { CandidateList } from "../../src/ui/CandidateList";

function candidate(id: number): IndexEntry {
  return {
    id, exposure_id: 1, phase: "Pn3m", basis: 0.1, score: 0.9, r_squared: null, lattice_d: 100,
    status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
  };
}

function setup() {
  render(
    <>
      <CandidateList rows={[candidate(1), candidate(2)].map((c) => ({ kind: "one" as const, candidate: c }))} assigned={new Set()} disabled={false}
        formFactor={false} onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />
      <button type="button">outside</button>
    </>,
  );
  const group = screen.getByRole("group", { name: "Phase choices" });
  const [first, second] = screen.getAllByRole("button", { pressed: false });
  return { group, first: first!, second: second! };
}

describe("CandidateList focus ring", () => {
  it("hides the ring on a choice focused by the mouse", () => {
    const { group, first } = setup();
    expect(group).not.toHaveAttribute("data-pointer-focus");
    fireEvent.pointerDown(first);
    first.focus();
    expect(group).toHaveAttribute("data-pointer-focus");
    expect(first).toHaveFocus();
  });

  it("shows the ring again once Tab or Shift+Tab moves through the choices", () => {
    const { group, first } = setup();
    fireEvent.pointerDown(first);
    first.focus();
    fireEvent.keyDown(first, { key: "ArrowRight" });
    expect(group).toHaveAttribute("data-pointer-focus");
    fireEvent.keyDown(first, { key: "Tab", shiftKey: true });
    expect(group).not.toHaveAttribute("data-pointer-focus");
    fireEvent.pointerDown(first);
    fireEvent.keyDown(first, { key: "Tab" });
    expect(group).not.toHaveAttribute("data-pointer-focus");
  });

  it("forgets the mouse once focus leaves the choices", () => {
    const { group, first, second } = setup();
    fireEvent.pointerDown(first);
    act(() => first.focus());
    act(() => second.focus());
    expect(group).toHaveAttribute("data-pointer-focus");
    act(() => screen.getByRole("button", { name: "outside" }).focus());
    expect(group).not.toHaveAttribute("data-pointer-focus");
  });
});
