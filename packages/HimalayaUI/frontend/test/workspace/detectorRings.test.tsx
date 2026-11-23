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
const view: DetectorRingsView = { size: { width: 981, height: 1043 }, cx: 421.9, cy: 205.6, side: "right", rings: [115, 230, 461].map((r) => ({ r, color: "red" })) };
const rings = (root: ParentNode) => [...root.querySelectorAll<SVGSVGElement>('[data-testid="detector-rings"]')];

describe("predicted rings on the detector images", () => {
  it("draws each ring as a half-circle on its side of the beam, so the other half shows the data", () => {
    const halves = (side: "left" | "right" | "top" | "bottom") => {
      const { container, unmount } = render(
        <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1)] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}}
          ringsOf={() => ({ ...view, side, rings: [{ r: 100, color: "red" }] })} />,
      );
      const svg = rings(container)[0]!;
      expect(svg.querySelectorAll("circle")).toHaveLength(0);
      const d = svg.querySelector('[data-role="det-ring"]')!.getAttribute("d");
      unmount();
      return d;
    };
    expect(halves("right")).toBe("M 421.9 105.6 A 100 100 0 0 1 421.9 305.6");
    expect(halves("left")).toBe("M 421.9 105.6 A 100 100 0 0 0 421.9 305.6");
    expect(halves("top")).toBe("M 321.9 205.6 A 100 100 0 0 1 521.9 205.6");
    expect(halves("bottom")).toBe("M 321.9 205.6 A 100 100 0 0 0 521.9 205.6");
  });

  it("draws the rings on every image card, cropped the way the card crops the image", () => {
    const { container } = render(
      <ExposureStrip groups={[{ key: 1, title: undefined, exposures: [exp(1), exp(2)] }]} currentId={1} onPick={() => {}} onToggleDrop={() => {}}
        ringsOf={() => view} />,
    );
    expect(rings(container)).toHaveLength(2);
    for (const svg of rings(container)) {
      expect(svg.getAttribute("viewBox")).toBe("0 0 981 1043");
      expect(svg.getAttribute("preserveAspectRatio")).toBe("xMidYMid slice");
      const arcs = [...svg.querySelectorAll('[data-role="det-ring"]')];
      expect(arcs.map((a) => a.getAttribute("d")!.split(" ")[4])).toEqual(["115", "230", "461"]);
      expect(arcs[0]!.getAttribute("stroke")).toBe("red");
      expect(arcs[0]!.getAttribute("stroke-dasharray")).toBe("4 3");
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

  it("reloads an image the browser cached before the server sent its size", async () => {
    const seen: [string, RequestCache | undefined][] = [];
    vi.stubGlobal("fetch", vi.fn(async (url: string, init?: RequestInit) => {
      seen.push([url, init?.cache]);
      const h = init?.cache === "reload" ? { "X-Image-Width": "981", "X-Image-Height": "1043" } : {};
      return { ok: true, headers: new Headers(h), body: null };
    }));
    const { getDetectorSize } = await import("../../src/api");
    expect(await getDetectorSize("/api/exposures/1/image?thumb=1")).toEqual({ width: 981, height: 1043 });
    expect(seen.map((s) => s[1])).toEqual(["default", "reload"]);
  });

  it("reads each image's raw detector size from its response headers, once per image", async () => {
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const wrapper = ({ children }: { children: ReactNode }) => <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
    const urls = ["/api/exposures/1/image?thumb=1", "/api/exposures/2/image?thumb=1"];
    const { result, rerender } = renderHook(() => useDetectorSizes(urls), { wrapper });
    await waitFor(() => expect(result.current.get(urls[0]!)).toEqual({ width: 981, height: 1043 }));
    expect(result.current.has(urls[1]!)).toBe(false);
    rerender();
    expect(fetched).toEqual([urls[0], urls[1], urls[1]]);
  });
});
