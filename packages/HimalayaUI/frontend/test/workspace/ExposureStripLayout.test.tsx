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

  it("puts the Images heading, the selection caption and Select all/none on one line", () => {
    const a = makeExposure({ id: 1 });
    const b = makeExposure({ id: 2 });
    const selection = { selectedIds: new Set([1, 2]), selectableIds: new Set([1, 2]), onToggle: () => {}, onRange: () => {}, onAll: () => {}, onNone: () => {} };
    render(
      <ExposureStrip
        title="Images"
        hint="2 of 3 kept"
        groups={[{ key: 1, title: undefined, exposures: [a, b], selection }]}
        currentId={1}
        onPick={() => {}}
        onToggleDrop={() => {}}
      />,
    );
    const heading = screen.getByRole("heading", { level: 2, name: "Images" });
    const row = heading.closest("[data-strip-heading]")!;
    expect(row).not.toBeNull();
    expect(within(row as HTMLElement).getByText("2 of 3 kept")).toBeTruthy();
    expect(within(row as HTMLElement).getByTestId("exposure-selection-caption")).toHaveTextContent("Average of 2 images");
    expect(within(row as HTMLElement).getByRole("button", { name: "Select none" })).toBeTruthy();
    expect(screen.getAllByTestId("exposure-selection-tools")).toHaveLength(1);
  });
});
