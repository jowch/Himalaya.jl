import { describe, it, expect, vi } from "vitest";
import { render, screen, within } from "@testing-library/react";
import { ExposureStrip } from "../../src/ui/ExposureStrip";
import { SampleList } from "../../src/ui/SampleList";
import { beamOffWarning } from "../../src/workspace/model";
import type { CorpusSample } from "../../src/api";
import { makeExposure } from "../helpers/exposure";

describe("beamOffWarning", () => {
  it("warns only for a kept image the server reads as beam-off", () => {
    expect(beamOffWarning(makeExposure({ beam_off: true }))).toBe("Beam likely off");
    expect(beamOffWarning(makeExposure({ beam_off: false }))).toBeUndefined();
    expect(beamOffWarning(makeExposure({ beam_off: true, status: "rejected" }))).toBeUndefined();
  });
});

describe("ExposureStrip warnings", () => {
  it("outlines a warned tile and names the reason under it, leaving other tiles plain", () => {
    const dark = makeExposure({ id: 1, beam_off: true });
    const lit = makeExposure({ id: 2 });
    render(
      <ExposureStrip
        groups={[{ key: 1, title: undefined, exposures: [dark, lit] }]}
        currentId={2}
        onPick={() => {}}
        onToggleDrop={() => {}}
        warningOf={beamOffWarning}
      />,
    );
    const tiles = screen.getAllByRole("listitem");
    expect(tiles[0]!.hasAttribute("data-warned")).toBe(true);
    expect(within(tiles[0]!).getByText("Beam likely off")).toBeTruthy();
    expect(screen.getByRole("button", { name: "Image 1, beam likely off" })).toBeTruthy();
    expect(tiles[1]!.hasAttribute("data-warned")).toBe(false);
    expect(within(tiles[1]!).queryByText("Beam likely off")).toBeNull();
  });
});

describe("SampleList beam-off dot", () => {
  it("marks a row with kept beam-off images and says how many", () => {
    Element.prototype.scrollIntoView = vi.fn();
    const rows = [
      { id: 1, name: "A", beam_off_images: 2 },
      { id: 2, name: "B", beam_off_images: 1 },
      { id: 3, name: "C", beam_off_images: 0 },
    ] as unknown as CorpusSample[];
    render(<SampleList samples={rows} currentId={3} filter="" onFilter={() => {}} onPick={() => {}} />);
    expect(screen.getByRole("img", { name: "2 images likely beam-off" })).toBeTruthy();
    expect(screen.getByRole("img", { name: "1 image likely beam-off" })).toBeTruthy();
    expect(screen.getByRole("button", { name: /^C/ }).querySelector("[data-warned]")).toBeNull();
  });
});
