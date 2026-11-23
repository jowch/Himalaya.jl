import { describe, it, expect } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ExposureStrip } from "../../src/ui/ExposureStrip";
import type { Exposure } from "../../src/api";

const exp = (id: number): Exposure => ({ id, status: "accepted", image_path: null } as unknown as Exposure);
const strip = (currentId: number) => (
  <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1), exp(2), exp(3)] }]} currentId={currentId} onPick={() => {}} onToggleDrop={() => {}} />
);
const image = (n: number) => screen.getByRole("button", { name: `Image ${n}` });
const drops = () => screen.getAllByRole("button", { name: "Drop" });

describe("ExposureStrip focus", () => {
  it("makes only the current image and its Drop button Tab stops", () => {
    render(strip(2));
    expect(image(1)).toHaveAttribute("tabindex", "-1");
    expect(image(2)).toHaveAttribute("tabindex", "0");
    expect(image(3)).toHaveAttribute("tabindex", "-1");
    expect(drops().map((b) => b.getAttribute("tabindex"))).toEqual(["-1", "0", "-1"]);
  });

  it("does not leave focus on an image or Drop button clicked with the mouse", () => {
    render(strip(1));
    image(2).focus();
    fireEvent.click(image(2), { detail: 1 });
    expect(image(2)).not.toHaveFocus();
    const drop = drops()[0]!;
    drop.focus();
    fireEvent.click(drop, { detail: 1 });
    expect(drop).not.toHaveFocus();
  });

  it("keeps focus on an image picked from the keyboard", () => {
    render(strip(1));
    image(2).focus();
    fireEvent.click(image(2), { detail: 0 });
    expect(image(2)).toHaveFocus();
  });

  it("moves focus to the new current image when the selection moves while an image has focus", () => {
    const { rerender } = render(strip(1));
    image(1).focus();
    rerender(strip(3));
    expect(image(3)).toHaveFocus();
  });

  it("moves focus to the new current image's Drop button when a Drop button had focus", () => {
    const { rerender } = render(strip(1));
    drops()[0]!.focus();
    rerender(strip(2));
    expect(drops()[1]).toHaveFocus();
  });
});
