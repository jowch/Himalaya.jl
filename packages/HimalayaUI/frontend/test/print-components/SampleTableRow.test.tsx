import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";
import type { GalleryExposure } from "../../src/print/components/ThumbnailGallery";
import type { Tag } from "../../src/print/ui/tag";

const EXPOSURES: GalleryExposure[] = [
  { id: 37, src: null, frameNo: "37", representative: true },
  { id: 65, src: null, frameNo: "65" },
  { id: 66, src: null, frameNo: "66", rejected: true },
];

const TAGS: Tag[] = [
  { key: "LL37" },
  { key: "temperature", value: "37C" },
];

function renderRow(overrides: Partial<Parameters<typeof SampleTableRow>[0]> = {}) {
  return render(
    <SampleTableRow
      name="Sample A"
      sampleId="s-001"
      exposures={EXPOSURES}
      kept={2}
      total={3}
      tags={TAGS}
      {...overrides}
    />
  );
}

describe("<SampleTableRow> structure", () => {
  it("(a) renders the root data-testid=sample-table-row", () => {
    renderRow();
    expect(screen.getByTestId("sample-table-row")).toBeInTheDocument();
  });

  it("(b) contains the SpecCell name and sampleId text", () => {
    renderRow();
    expect(screen.getByText("Sample A")).toBeInTheDocument();
    expect(screen.getByText("s-001")).toBeInTheDocument();
  });

  it("(c) contains a thumbnail-gallery with one thumbnail per exposure", () => {
    renderRow();
    const gallery = screen.getByTestId("thumbnail-gallery");
    expect(gallery).toBeInTheDocument();
    const thumbs = screen.getAllByTestId("thumbnail");
    expect(thumbs).toHaveLength(EXPOSURES.length);
  });

  it("renders exactly 5 grid-child cell wrappers", () => {
    renderRow();
    const row = screen.getByTestId("sample-table-row");
    const grid = row.firstElementChild!;
    expect(grid.children).toHaveLength(5);
  });

  it("(d) contains the KeptCell counts (kept and total)", () => {
    const { container } = renderRow();
    expect(container.querySelector("[data-role='kept-count']")).toHaveTextContent("2");
    expect(container.querySelector("[data-role='kept-total']")).toHaveTextContent("/ 3");
  });

  it("(h) renders the tags via TagList", () => {
    renderRow();
    expect(screen.getByText("LL37")).toBeInTheDocument();
    expect(screen.getByText("temperature")).toBeInTheDocument();
  });
});

describe("<SampleTableRow> status column", () => {
  it("(e) phase set → PhaseChip present and status-unset absent", () => {
    const { container } = renderRow({ phase: "Pn3m" });
    expect(screen.getByTestId("phase-chip")).toBeInTheDocument();
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
    expect(container.querySelector("[data-role='status-unset']")).not.toBeInTheDocument();
  });

  it("(e) phase null → status-unset present", () => {
    const { container } = renderRow({ phase: null });
    expect(container.querySelector("[data-role='status-unset']")).toBeInTheDocument();
  });
});

describe("<SampleTableRow> screened reflection", () => {
  it("(f) data-screened is 'true' when screened", () => {
    renderRow({ screened: true });
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("data-screened", "true");
  });

  it("(f) data-screened is 'false' when not screened", () => {
    renderRow({ screened: false });
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("data-screened", "false");
  });

  it("(f) data-screened is 'false' when screened is omitted", () => {
    renderRow();
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("data-screened", "false");
  });
});

describe("<SampleTableRow> exposure selection", () => {
  it("(g) selecting a thumb bubbles onSelectExposure with that exposure id", () => {
    const onSelectExposure = vi.fn();
    renderRow({ onSelectExposure });
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[1]);
    expect(onSelectExposure).toHaveBeenCalledTimes(1);
    expect(onSelectExposure).toHaveBeenCalledWith(65);
  });

  it("selectedExposureIds highlights EVERY matching thumb (multi-select cull model)", () => {
    // EXPOSURES order: [37, 65, 66] — select the first and last.
    renderRow({ selectedExposureIds: new Set([37, 66]) });
    const thumbs = screen.getAllByTestId("thumbnail");
    expect((thumbs[0].getAttribute("data-state") ?? "").split(" ")).toContain("selected");
    expect((thumbs[1].getAttribute("data-state") ?? "").split(" ")).not.toContain("selected");
    expect((thumbs[2].getAttribute("data-state") ?? "").split(" ")).toContain("selected");
  });
});

// ---------------------------------------------------------------------------
// Nav seams (Task 3)
// ---------------------------------------------------------------------------

const baseProps = {
  name: "POPC",
  sampleId: "#42",
  exposures: [{ id: 1, src: null, frameNo: 1 }],
  kept: 1,
  total: 1,
  tags: [] as import("../../src/print/ui/tag").Tag[],
};

describe("<SampleTableRow> nav seams", () => {
  it("forwards onOpenLoupe to the SpecCell name button", () => {
    const onOpenLoupe = vi.fn();
    render(<SampleTableRow {...baseProps} onOpenLoupe={onOpenLoupe} />);
    fireEvent.click(screen.getByRole("button", { name: /POPC/ }));
    expect(onOpenLoupe).toHaveBeenCalled();
  });

  it("forwards onActivateExposure to the gallery on thumb double-click", () => {
    const onActivateExposure = vi.fn();
    render(<SampleTableRow {...baseProps} onActivateExposure={onActivateExposure} />);
    fireEvent.doubleClick(screen.getAllByTestId("thumbnail")[0]);
    expect(onActivateExposure).toHaveBeenCalledWith(1);
  });

  it("makes the status cell a Focus door (button) when onOpenFocus is set", () => {
    const onOpenFocus = vi.fn();
    render(<SampleTableRow {...baseProps} phase={null} onOpenFocus={onOpenFocus} />);
    fireEvent.click(screen.getByRole("button", { name: /index/i }));
    expect(onOpenFocus).toHaveBeenCalled();
  });

  it("status cell is NOT a button when onOpenFocus is absent", () => {
    render(<SampleTableRow {...baseProps} phase={null} />);
    expect(screen.queryByRole("button", { name: /index/i })).toBeNull();
  });
});
