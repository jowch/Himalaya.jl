import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { ExposureStrip } from "../../src/ui/ExposureStrip";
import { DetectorOverlay } from "../../src/ui/DetectorOverlay";
import { detectorGeometry } from "../../src/lib/detectorRings";
import type { Exposure, Setup } from "../../src/api";

const exp = (id: number) => ({ id, status: "kept", image_path: `/x/${id}.tif`, image_version: "v1", filename: `f${id}` } as unknown as Exposure);
const geometry = detectorGeometry({
  beam_center_x: 421.409, beam_center_y: 836.946, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095, q_units: "A-1",
} as Setup);

let fetched: string[] = [];
beforeEach(() => {
  fetched = [];
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    fetched.push(url);
    return { ok: true, headers: new Headers({ "X-Image-Width": "981", "X-Image-Height": "1043" }), body: null };
  }));
});

function wrap(ui: ReactNode) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

const rings = (root: ParentNode) => [...root.querySelectorAll<SVGSVGElement>('[data-testid="detector-rings"]')];

describe("predicted rings on the detector images", () => {
  it("draws the hovered candidate's rings on every image card, cropped the way the card crops the image", async () => {
    const { container } = wrap(
      <ExposureStrip
        groups={[{ key: 1, title: undefined, exposures: [exp(1), exp(2)] }]}
        currentId={1} onPick={() => {}} onToggleDrop={() => {}}
        rings={{ color: "red", qs: [0.05, 0.1, 0.2] }} geometryOf={() => geometry}
      />,
    );
    await waitFor(() => expect(rings(container)).toHaveLength(2));
    for (const svg of rings(container)) {
      expect(svg.getAttribute("viewBox")).toBe("0 0 981 1043");
      expect(svg.getAttribute("preserveAspectRatio")).toBe("xMidYMid slice");
      const circles = [...svg.querySelectorAll("circle")];
      expect(circles).toHaveLength(3);
      expect(Number(circles[0]!.getAttribute("cx"))).toBeCloseTo(421.909, 3);
      expect(Number(circles[0]!.getAttribute("cy"))).toBeCloseTo(205.554, 3);
      expect(circles[1]!.getAttribute("stroke")).toBe("red");
    }
    expect(fetched.every((u) => u.includes("thumb=1"))).toBe(true);
  });

  it("draws them on the enlarged image, letterboxed the way the image is", async () => {
    wrap(
      <DetectorOverlay exposures={[exp(1)]} currentId={1} onStep={() => {}} onToggleDrop={undefined} onClose={() => {}}
        rings={{ color: "red", qs: [0.1] }} geometryOf={() => geometry} />,
    );
    await waitFor(() => expect(rings(document)).toHaveLength(1));
    expect(rings(document)[0]!.getAttribute("preserveAspectRatio")).toBe("xMidYMid meet");
  });

  it("draws nothing without a hovered candidate, or for a setup without its geometry", async () => {
    const { container, rerender } = wrap(
      <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1)] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}}
        rings={null} geometryOf={() => geometry} />,
    );
    expect(rings(container)).toHaveLength(0);
    rerender(
      <QueryClientProvider client={new QueryClient()}>
        <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1)] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}}
          rings={{ color: "red", qs: [0.1] }} geometryOf={() => null} />
      </QueryClientProvider>,
    );
    await new Promise((r) => setTimeout(r, 20));
    expect(rings(container)).toHaveLength(0);
    expect(fetched).toEqual([]);
  });
});
