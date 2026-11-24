import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { ExposureStrip } from "../../src/ui/ExposureStrip";
import { DetectorOverlay } from "../../src/ui/DetectorOverlay";
import { DetectorLutFilters } from "../../src/detector/DetectorLutFilters";
import { buildPrintDetectorLut } from "../../src/detector/detectorLut";
import type { Exposure } from "../../src/api";

const exp = (id: number) => ({ id, status: "kept", image_path: `/x/${id}.tif`, image_version: "v1", filename: `f${id}` } as unknown as Exposure);

function referencedFilter(img: HTMLImageElement): Element | null {
  const id = /^url\(["']?#([^"')]+)["']?\)$/.exec(img.style.filter)?.[1];
  return id ? document.getElementById(id) : null;
}

function expectMagma(imgs: HTMLImageElement[]): void {
  expect(imgs.length).toBeGreaterThan(0);
  const magma = buildPrintDetectorLut("magma");
  const top = (magma[255 * 3]! / 255).toFixed(4);
  for (const img of imgs) {
    const filter = referencedFilter(img);
    expect(filter, `filter for ${img.style.filter} is mounted`).not.toBeNull();
    expect(filter!.querySelector("feFuncR")!.getAttribute("tableValues")!.split(" ").at(-1)).toBe(top);
  }
}

describe("workspace detector images", () => {
  it("colour the frame thumbnails with a mounted magma filter", () => {
    const { container } = render(
      <>
        <DetectorLutFilters />
        <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1), exp(2)] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}} />
      </>,
    );
    expectMagma([...container.querySelectorAll("img")]);
  });

  it("colour the enlarged detector view with the same mounted filter", () => {
    render(
      <>
        <DetectorLutFilters />
        <DetectorOverlay exposures={[exp(1)]} currentId={1} onStep={() => {}} onToggleDrop={undefined} onClose={() => {}} />
      </>,
    );
    expectMagma([...document.querySelectorAll<HTMLImageElement>("img[data-exposure-id]")]);
  });
});
