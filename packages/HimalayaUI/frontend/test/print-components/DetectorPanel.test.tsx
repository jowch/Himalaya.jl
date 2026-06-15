import { render, screen, fireEvent, within } from "@testing-library/react";
import { vi, beforeEach } from "vitest";
import { DetectorPanel } from "../../src/print/components/DetectorPanel";

const TINY_PNG =
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwADhQGAWjR9awAAAABJRU5ErkJggg==";

beforeEach(() => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    headers: { get: (k: string) => (k === "X-Image-Width" ? "2048" : k === "X-Image-Height" ? "1024" : null) },
    blob: () => Promise.resolve(new Blob(
      [Uint8Array.from(atob(TINY_PNG), (c) => c.charCodeAt(0))], { type: "image/png" })),
  } as unknown as Response);
  global.createImageBitmap = vi.fn().mockResolvedValue({ width: 1, height: 1, close: vi.fn() } as unknown as ImageBitmap);
  const mockOffscreen = {
    getContext: () => ({ drawImage: vi.fn(), getImageData: () => ({ data: new Uint8ClampedArray(4) }) }),
  };
  // @ts-expect-error JSDOM stub
  global.OffscreenCanvas = vi.fn().mockImplementation(() => mockOffscreen);
});

const RINGS = [
  { q: 0.045, color: "var(--color-pn3m)" },
  { q: 0.064, color: "var(--color-pn3m)" },
];

describe("DetectorPanel", () => {
  it("renders the default label and the detector frame", () => {
    render(<DetectorPanel src={null} />);
    expect(screen.getByTestId("detector-panel")).toBeInTheDocument();
    expect(screen.getByText("Detector image")).toBeInTheDocument();
    expect(screen.getByTestId("detector-frame")).toBeInTheDocument();
  });
  it("renders a custom label", () => {
    render(<DetectorPanel src={null} label="Real source" />);
    expect(screen.getByText("Real source")).toBeInTheDocument();
  });
  it("renders the header tools slot (exposure strip)", () => {
    render(<DetectorPanel src={null} tools={<div data-testid="expo">e</div>} />);
    expect(screen.getByTestId("expo")).toBeInTheDocument();
  });
  it("renders a hint line when provided", () => {
    render(<DetectorPanel src={null} hint="The real source." />);
    expect(screen.getByText("The real source.")).toBeInTheDocument();
  });
  it("forwards a placement-only className", () => {
    render(<DetectorPanel src={null} className="h-full" />);
    expect(screen.getByTestId("detector-panel").className).toContain("h-full");
  });
  it("overlays one phase ring per provided reflection", () => {
    const { container } = render(<DetectorPanel src={null} rings={RINGS} />);
    expect(screen.getByTestId("detector-rings")).toBeInTheDocument();
    expect(container.querySelectorAll('[data-role="det-ring"]')).toHaveLength(2);
  });
  it("renders no ring overlay when no rings are provided", () => {
    render(<DetectorPanel src={null} />);
    expect(screen.queryByTestId("detector-rings")).not.toBeInTheDocument();
  });
  it("lights the ring whose q matches hoveredQ (the triple-link)", () => {
    const { container } = render(
      <DetectorPanel src={null} rings={RINGS} hoveredQ={0.045} />,
    );
    const hot = container.querySelector('[data-role="det-ring"][data-hot="true"]');
    expect(hot).not.toBeNull();
    expect(hot?.getAttribute("data-ring-q")).toBe("0.045");
  });
  it("places rings on a measured off-center beamCenter override", () => {
    const { container } = render(
      <DetectorPanel src={null} rings={RINGS} beamCenter={{ x: 0.43, y: 0.2 }} />,
    );
    const ring = container.querySelector('[data-role="ring-sharp"]');
    expect(ring?.getAttribute("cx")).toBe("0.43");
    expect(ring?.getAttribute("cy")).toBe("0.2");
  });
  it("captions the rings with a phase chip + frame-scope copy when ringPhases is provided", () => {
    render(
      <DetectorPanel src={null} rings={RINGS} ringPhases={["Hexagonal"]} />,
    );
    const caption = screen.getByTestId("detector-ring-caption");
    expect(within(caption).getByTestId("phase-chip")).toHaveTextContent(
      "Hexagonal",
    );
    // Reading order: the label LEADS ("rings: Hexagonal this frame's
    // indexing"), so the a11y string never opens with a bare phase name.
    expect(caption.textContent).toMatch(/^rings:/);
    expect(caption).toHaveTextContent("this frame's indexing");
    // The caption is the WCAG 1.4.1 second channel — it must live in the a11y
    // tree, never hidden.
    expect(caption.closest('[aria-hidden="true"]')).toBeNull();
  });

  it("renders one chip per distinct phase, in the given (rail) order", () => {
    render(
      <DetectorPanel
        src={null}
        rings={RINGS}
        ringPhases={["Pn3m", "Im3m"]}
      />,
    );
    const caption = screen.getByTestId("detector-ring-caption");
    const chips = within(caption).getAllByTestId("phase-chip");
    expect(chips.map((c) => c.textContent)).toEqual(["Pn3m", "Im3m"]);
  });

  it("omits the caption entirely when ringPhases is empty or absent", () => {
    const { rerender } = render(<DetectorPanel src={null} ringPhases={[]} />);
    expect(screen.queryByTestId("detector-ring-caption")).toBeNull();
    rerender(<DetectorPanel src={null} rings={RINGS} />);
    expect(screen.queryByTestId("detector-ring-caption")).toBeNull();
  });

  it("fires onHoverQ with the ring's q when a ring is hovered", () => {
    const onHoverQ = vi.fn();
    const { container } = render(
      <DetectorPanel src={null} rings={RINGS} onHoverQ={onHoverQ} />,
    );
    const hit = container.querySelector('[data-role="ring-hit"]');
    fireEvent.mouseEnter(hit!);
    expect(onHoverQ).toHaveBeenCalledWith(0.045);
  });

  test("passes orient through to DetectorRings", async () => {
    render(<DetectorPanel src="/x.png" rings={[0.1, 0.2]} orient="landscape" />);
    const rings = await screen.findByTestId("detector-rings");
    expect(rings.getAttribute("data-orient")).toBe("landscape");
  });

  test("forwards onRawSize/onOrient to DetectorImage (props accepted, no crash)", () => {
    const onRawSize = vi.fn(); const onOrient = vi.fn();
    expect(() =>
      render(<DetectorPanel src="/x.png" onRawSize={onRawSize} onOrient={onOrient} />),
    ).not.toThrow();
  });
});
