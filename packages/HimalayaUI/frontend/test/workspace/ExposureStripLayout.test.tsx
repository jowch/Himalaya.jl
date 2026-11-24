import { describe, it, expect } from "vitest";
import { render, screen, within } from "@testing-library/react";
import { ExposureStrip } from "../../src/ui/ExposureStrip";
import { beamOffWarning } from "../../src/workspace/model";
import { makeExposure } from "../helpers/exposure";

describe("ExposureStrip tile layout", () => {
  it("gives every tile one footer line, so a warned or dropped image lines up with the rest", () => {
    const plain = makeExposure({ id: 1 });
    const dark = makeExposure({ id: 2, beam_off: true });
    const dropped = makeExposure({ id: 3, status: "rejected" });
    render(
      <ExposureStrip
        groups={[{ key: 1, title: undefined, exposures: [plain, dark, dropped] }]}
        currentId={1}
        onPick={() => {}}
        onToggleDrop={() => {}}
        warningOf={beamOffWarning}
      />,
    );
    const tiles = screen.getAllByRole("listitem");
    for (const tile of tiles) {
      expect(tile.querySelectorAll("[data-tile-footer]")).toHaveLength(1);
      expect(tile.children).toHaveLength(2);
    }
    const footer = (i: number) => tiles[i]!.querySelector<HTMLElement>("[data-tile-footer]")!;
    expect(within(footer(1)).getByText("Beam likely off")).toBeTruthy();
    expect(within(footer(2)).getByText("Dropped")).toBeTruthy();
    expect(within(footer(2)).getByRole("button", { name: "Restore" })).toBeTruthy();
    expect(within(footer(0)).getByRole("button", { name: "Drop" })).toBeTruthy();
    expect(within(footer(0)).getByText("Image 1")).toBeTruthy();
  });
});
