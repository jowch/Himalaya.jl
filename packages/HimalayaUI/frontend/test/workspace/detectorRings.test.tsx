import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { ExposureStrip } from "../../src/ui/ExposureStrip";
import { DetectorOverlay } from "../../src/ui/DetectorOverlay";
import type { DetectorRingsView } from "../../src/ui/DetectorRings";
import { useDetectorSizes } from "../../src/queries";
import type { Exposure } from "../../src/api";

const exp = (id: number) => ({ id, status: "kept", image_path: `/x/${id}.tif`, image_version: "v1", filename: `f${id}` } as unknown as Exposure);
const view: DetectorRingsView = { size: { width: 981, height: 1043 }, cx: 421.9, cy: 205.6, rings: [115, 230, 461].map((r) => ({ r, color: "red" })) };
const rings = (root: ParentNode) => [...root.querySelectorAll<SVGSVGElement>('[data-testid="detector-rings"]')];

describe("predicted rings on the detector images", () => {
  it("draws the rings on every image card, cropped the way the card crops the image", () => {
    const { container } = render(
      <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1), exp(2)] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}}
        ringsOf={() => view} />,
    );
    expect(rings(container)).toHaveLength(2);
    for (const svg of rings(container)) {
      expect(svg.getAttribute("viewBox")).toBe("0 0 981 1043");
      expect(svg.getAttribute("preserveAspectRatio")).toBe("xMidYMid slice");
      const circles = [...svg.querySelectorAll("circle")];
      expect(circles.map((c) => c.getAttribute("r"))).toEqual(["115", "230", "461"]);
      expect(circles[0]!.getAttribute("cx")).toBe("421.9");
      expect(circles[0]!.getAttribute("cy")).toBe("205.6");
      expect(circles[0]!.getAttribute("stroke")).toBe("red");
    }
  });

  it("draws them on the enlarged image, letterboxed the way the image is", () => {
    render(<DetectorOverlay exposures={[exp(1)]} currentId={1} onStep={() => {}} onToggleDrop={undefined} onClose={() => {}} ringsOf={() => view} />);
    expect(rings(document)).toHaveLength(1);
    expect(rings(document)[0]!.getAttribute("preserveAspectRatio")).toBe("xMidYMid meet");
  });

  it("draws nothing on an image the page gives no rings, or on one without a TIFF", () => {
    const noTiff = { ...exp(3), image_path: null } as Exposure;
    const { container } = render(
      <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1), noTiff] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}}
        ringsOf={(e) => (e.id === 1 ? null : view)} />,
    );
    expect(rings(container)).toHaveLength(0);
  });
});

describe("useDetectorSizes", () => {
  let fetched: string[] = [];
  beforeEach(() => {
    fetched = [];
    vi.stubGlobal("fetch", vi.fn(async (url: string) => {
      fetched.push(url);
      const h = url.includes("/2/") ? {} : { "X-Image-Width": "981", "X-Image-Height": "1043" };
      return { ok: true, headers: new Headers(h), body: null };
    }));
  });

  it("reads each image's raw detector size from its response headers, once per image", async () => {
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const wrapper = ({ children }: { children: ReactNode }) => <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
    const urls = ["/api/exposures/1/image?thumb=1", "/api/exposures/2/image?thumb=1"];
    const { result, rerender } = renderHook(() => useDetectorSizes(urls), { wrapper });
    await waitFor(() => expect(result.current.get(urls[0]!)).toEqual({ width: 981, height: 1043 }));
    expect(result.current.has(urls[1]!)).toBe(false);
    rerender();
    expect(fetched).toEqual(urls);
  });
});
