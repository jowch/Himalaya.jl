import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ThumbnailGallery, type GalleryExposure } from "../../src/print/components/ThumbnailGallery";

const THREE_EXPOSURES: GalleryExposure[] = [
  { id: 37, src: null, frameNo: "37" },
  { id: 65, src: null, frameNo: "65" },
  { id: 66, src: null, frameNo: "66" },
];

describe("<ThumbnailGallery> root element", () => {
  it("renders [data-testid='thumbnail-gallery']", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getByTestId("thumbnail-gallery")).toBeInTheDocument();
  });

  it("defaults data-size to 'sm' and data-align to 'start'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    const gallery = screen.getByTestId("thumbnail-gallery");
    expect(gallery).toHaveAttribute("data-size", "sm");
    expect(gallery).toHaveAttribute("data-align", "start");
  });

  it("the overflow-x-auto scroller carries tabindex=-1 (must NOT be an auto Tab-stop)", () => {
    // Chrome auto-adds scrollable containers with overflowing content to the
    // sequential Tab order even without a tabindex; tabindex=-1 opts this
    // scroller out so it never becomes a stray keyboard stop (e.g. a redundant
    // stop on the loupe filmstrip). Children stay independently tabbable.
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("tabindex", "-1");
  });
});

describe("<ThumbnailGallery> exposure rendering", () => {
  it("renders one [data-testid='thumbnail'] per exposure", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getAllByTestId("thumbnail")).toHaveLength(3);
  });

  it("renders zero thumbnails for an empty exposure array", () => {
    render(<ThumbnailGallery exposures={[]} />);
    expect(screen.queryAllByTestId("thumbnail")).toHaveLength(0);
  });
});

describe("<ThumbnailGallery> onSelect", () => {
  it("clicking the second thumb fires onSelect with that exposure's id", () => {
    const onSelect = vi.fn();
    render(
      <ThumbnailGallery
        exposures={THREE_EXPOSURES}
        onSelect={onSelect}
      />
    );
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[1]); // second thumb → id=65
    expect(onSelect).toHaveBeenCalledTimes(1);
    expect(onSelect).toHaveBeenCalledWith(65);
  });

  it("does not throw when clicked without an onSelect handler", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    const thumbs = screen.getAllByTestId("thumbnail");
    expect(() => fireEvent.click(thumbs[0])).not.toThrow();
  });
});

describe("<ThumbnailGallery> selectedId", () => {
  it("marks exactly the matching thumb as selected; the others are not", () => {
    render(
      <ThumbnailGallery
        exposures={THREE_EXPOSURES}
        selectedId={65}
      />
    );
    const thumbs = screen.getAllByTestId("thumbnail");
    // thumb[0] = id 37, thumb[1] = id 65, thumb[2] = id 66
    const state0 = thumbs[0].getAttribute("data-state") ?? "";
    const state1 = thumbs[1].getAttribute("data-state") ?? "";
    const state2 = thumbs[2].getAttribute("data-state") ?? "";

    expect(state0.split(" ")).not.toContain("selected");
    expect(state1.split(" ")).toContain("selected");
    expect(state2.split(" ")).not.toContain("selected");
  });

  it("no thumb is selected when selectedId is undefined", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    const thumbs = screen.getAllByTestId("thumbnail");
    for (const thumb of thumbs) {
      const state = thumb.getAttribute("data-state") ?? "";
      expect(state.split(" ")).not.toContain("selected");
    }
  });
});

describe("<ThumbnailGallery> size", () => {
  it("size='lg' → gallery data-size='lg' AND every child data-size='lg'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} size="lg" />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-size", "lg");
    const thumbs = screen.getAllByTestId("thumbnail");
    for (const thumb of thumbs) {
      expect(thumb).toHaveAttribute("data-size", "lg");
    }
  });

  it("default size → gallery data-size='sm' AND every child data-size='sm'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-size", "sm");
    const thumbs = screen.getAllByTestId("thumbnail");
    for (const thumb of thumbs) {
      expect(thumb).toHaveAttribute("data-size", "sm");
    }
  });
});

describe("<ThumbnailGallery> align", () => {
  it("align='center' → gallery data-align='center'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} align="center" />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-align", "center");
  });

  it("default align='start' → data-align='start'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-align", "start");
  });
});

describe("<ThumbnailGallery> per-exposure flag passthrough", () => {
  it("a representative exposure → its thumb renders [data-role='thumb-rep']", () => {
    const exposures: GalleryExposure[] = [
      { id: 37, src: null, frameNo: "37", representative: true },
      { id: 65, src: null, frameNo: "65" },
    ];
    const { container } = render(<ThumbnailGallery exposures={exposures} />);
    const repMarkers = container.querySelectorAll("[data-role='thumb-rep']");
    expect(repMarkers).toHaveLength(1);
  });

  it("a rejected exposure → its thumb renders [data-testid='reject-overlay']", () => {
    const exposures: GalleryExposure[] = [
      { id: 37, src: null, frameNo: "37" },
      { id: 65, src: null, frameNo: "65", rejected: true },
    ];
    render(<ThumbnailGallery exposures={exposures} />);
    const overlays = screen.getAllByTestId("reject-overlay");
    expect(overlays).toHaveLength(1);
  });

  it("a non-rejected, non-representative exposure has neither marker", () => {
    const { container } = render(
      <ThumbnailGallery exposures={[{ id: 37, src: null }]} />
    );
    expect(container.querySelector("[data-role='thumb-rep']")).not.toBeInTheDocument();
    expect(screen.queryByTestId("reject-overlay")).not.toBeInTheDocument();
  });

  it("a kept exposure → its thumb renders [data-role='thumb-kept']; unscreened thumbs do not", () => {
    const exposures: GalleryExposure[] = [
      { id: 37, src: null, frameNo: "37", kept: true },
      { id: 65, src: null, frameNo: "65" }, // unscreened — no kept flag
    ];
    const { container } = render(<ThumbnailGallery exposures={exposures} />);
    expect(container.querySelectorAll("[data-role='thumb-kept']")).toHaveLength(1);
    const thumbs = screen.getAllByTestId("thumbnail");
    expect((thumbs[0].getAttribute("data-state") ?? "").split(" ")).toContain("kept");
    expect(thumbs[1]).toHaveAttribute("data-state", "normal");
  });
});

describe("<ThumbnailGallery> onActivate (double-click → loupe)", () => {
  it("fires onActivate(id) on double-click of a thumb", () => {
    const onActivate = vi.fn();
    render(
      <ThumbnailGallery
        exposures={[{ id: 7, src: null, frameNo: 1 }]}
        onActivate={onActivate}
      />,
    );
    fireEvent.doubleClick(screen.getAllByTestId("thumbnail")[0]);
    expect(onActivate).toHaveBeenCalledWith(7);
  });

  it("does not throw when double-clicked without an onActivate handler", () => {
    render(<ThumbnailGallery exposures={[{ id: 7, src: null, frameNo: 1 }]} />);
    expect(() =>
      fireEvent.doubleClick(screen.getAllByTestId("thumbnail")[0]),
    ).not.toThrow();
  });
});

describe("<ThumbnailGallery> selectedIds (multi-select cull model)", () => {
  it("marks every thumb whose id is in selectedIds as selected", () => {
    render(
      <ThumbnailGallery exposures={THREE_EXPOSURES} selectedIds={new Set([37, 66])} />
    );
    const thumbs = screen.getAllByTestId("thumbnail");
    // thumb[0]=37 (in set), thumb[1]=65 (out), thumb[2]=66 (in set)
    expect((thumbs[0].getAttribute("data-state") ?? "").split(" ")).toContain("selected");
    expect((thumbs[1].getAttribute("data-state") ?? "").split(" ")).not.toContain("selected");
    expect((thumbs[2].getAttribute("data-state") ?? "").split(" ")).toContain("selected");
  });

  it("ORs with selectedId — a thumb selected by either source is marked", () => {
    render(
      <ThumbnailGallery
        exposures={THREE_EXPOSURES}
        selectedId={65}
        selectedIds={new Set([37])}
      />
    );
    const thumbs = screen.getAllByTestId("thumbnail");
    expect((thumbs[0].getAttribute("data-state") ?? "").split(" ")).toContain("selected"); // via set
    expect((thumbs[1].getAttribute("data-state") ?? "").split(" ")).toContain("selected"); // via selectedId
    expect((thumbs[2].getAttribute("data-state") ?? "").split(" ")).not.toContain("selected");
  });
});
