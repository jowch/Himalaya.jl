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

  it("defaults data-variant to 'sheet'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-variant", "sheet");
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

describe("<ThumbnailGallery> variant", () => {
  it("variant='loupe' → gallery data-variant='loupe' AND children data-size='loupe'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} variant="loupe" />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-variant", "loupe");
    const thumbs = screen.getAllByTestId("thumbnail");
    for (const thumb of thumbs) {
      expect(thumb).toHaveAttribute("data-size", "loupe");
    }
  });

  it("default variant → gallery data-variant='sheet' AND children data-size='sheet'", () => {
    render(<ThumbnailGallery exposures={THREE_EXPOSURES} />);
    expect(screen.getByTestId("thumbnail-gallery")).toHaveAttribute("data-variant", "sheet");
    const thumbs = screen.getAllByTestId("thumbnail");
    for (const thumb of thumbs) {
      expect(thumb).toHaveAttribute("data-size", "sheet");
    }
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
});
