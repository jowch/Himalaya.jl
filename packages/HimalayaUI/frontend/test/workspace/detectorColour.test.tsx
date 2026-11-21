import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { ExposureStrip } from "../../src/workspace/ui/ExposureStrip";
import { DetectorOverlay } from "../../src/workspace/ui/DetectorOverlay";
import type { Exposure } from "../../src/api";

const exp = (id: number) => ({ id, status: "kept", image_path: `/x/${id}.tif`, filename: `f${id}` } as unknown as Exposure);

describe("workspace detector images", () => {
  it("colour the frame thumbnails with the magma ramp", () => {
    const { container } = render(<ExposureStrip exposures={[exp(1), exp(2)]} currentId={1} onPick={() => {}} onToggleDrop={() => {}} />);
    const imgs = [...container.querySelectorAll("img")];
    expect(imgs.length).toBeGreaterThan(0);
    for (const img of imgs) expect(img.style.filter).toBe("url(#detector-lut-magma)");
  });

  it("colour the enlarged detector view with the same ramp", () => {
    render(<DetectorOverlay exposures={[exp(1)]} currentId={1} onStep={() => {}} onToggleDrop={undefined} onClose={() => {}} />);
    const imgs = [...document.querySelectorAll<HTMLImageElement>("img[data-exposure-id]")];
    expect(imgs.length).toBeGreaterThan(0);
    for (const img of imgs) expect(img.style.filter).toBe("url(#detector-lut-magma)");
  });
});
