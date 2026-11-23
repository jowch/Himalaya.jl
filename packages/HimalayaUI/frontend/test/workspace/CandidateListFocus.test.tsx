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

const list = (ids: readonly number[]) => (
  <>
    <CandidateList rows={ids.map((id) => ({ kind: "one" as const, candidate: candidate(id) }))} assigned={new Set()} disabled={false}
      formFactor={false} onToggle={() => {}} onToggleGroup={() => {}} onFormFactor={() => {}} onHover={() => {}} />
    <button type="button">outside</button>
  </>
);

function setup() {
  const view = render(list([1, 2]));
  const group = screen.getByRole("group", { name: "Phase choices" });
  const [first, second] = screen.getAllByRole("button", { pressed: false });
  return { group, first: first!, second: second!, rerender: view.rerender };
}

describe("CandidateList focus ring", () => {
  it("hides the ring on a choice focused by the mouse", () => {
    const { group, first } = setup();
    expect(group).not.toHaveAttribute("data-pointer-focus");
    fireEvent.pointerDown(first);
    act(() => first.focus());
    fireEvent.pointerUp(first);
    expect(group).toHaveAttribute("data-pointer-focus");
    expect(first).toHaveFocus();
  });

  it("keeps the ring hidden through other keys, and shows it on the next choice focused without the mouse", () => {
    const { group, first, second } = setup();
    fireEvent.pointerDown(first);
    act(() => first.focus());
    fireEvent.pointerUp(first);
    fireEvent.keyDown(first, { key: "ArrowRight" });
    expect(group).toHaveAttribute("data-pointer-focus");
    act(() => second.focus());
    expect(group).not.toHaveAttribute("data-pointer-focus");
  });

  it("hides the ring when the mouse presses a choice that already has focus", () => {
    const { group, first } = setup();
    act(() => first.focus());
    expect(group).not.toHaveAttribute("data-pointer-focus");
    fireEvent.pointerDown(first);
    expect(group).toHaveAttribute("data-pointer-focus");
  });

  it("shows the ring on keyboard focus after the mouse-focused choice was replaced", () => {
    const { first, rerender } = setup();
    fireEvent.pointerDown(first);
    act(() => first.focus());
    fireEvent.pointerUp(first);
    rerender(list([3, 4]));
    expect(document.body).toHaveFocus();
    const group = screen.getByRole("group", { name: "Phase choices" });
    act(() => screen.getAllByRole("button", { pressed: false })[0]!.focus());
    expect(group).not.toHaveAttribute("data-pointer-focus");
  });

  it("forgets the mouse once focus leaves the choices", () => {
    const { group, first, second } = setup();
    fireEvent.pointerDown(first);
    act(() => first.focus());
    fireEvent.pointerUp(first);
    fireEvent.pointerDown(second);
    act(() => second.focus());
    fireEvent.pointerUp(second);
    expect(group).toHaveAttribute("data-pointer-focus");
    act(() => screen.getByRole("button", { name: "outside" }).focus());
    expect(group).not.toHaveAttribute("data-pointer-focus");
  });
});
