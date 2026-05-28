import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import type { CorpusSample, Exposure } from "../src/api";
import { LoupeSidebar } from "../src/components/LoupeSidebar";

function sample(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 7, experiment_id: 3, name: "JC042", display_name: "DOPE 80%",
    notes: null, tags: [], q_units: "A-1", ...over,
  };
}

function exposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 7, filename: "JC042-001.dat", kind: "file",
    selected: false, status: null, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
    ...over,
  };
}

function defaultProps() {
  const e = exposure();
  return {
    exposure: e,
    exposures: [e],
    sample: sample(),
    onDropToggle: vi.fn(),
    onSetRepresentative: vi.fn(),
  };
}

describe("LoupeSidebar — meta-list", () => {
  it("shows filename, kind, frame position and status", () => {
    const exposures = [
      exposure({ id: 1 }), exposure({ id: 2, kind: "averaged", filename: null }),
    ];
    render(
      <LoupeSidebar
        {...defaultProps()}
        exposure={exposures[1]}
        exposures={exposures}
      />,
    );
    expect(screen.getByTestId("loupe-meta-filename")).toHaveTextContent("—");
    expect(screen.getByTestId("loupe-meta-kind")).toHaveTextContent("averaged");
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    expect(screen.getByTestId("loupe-meta-status")).toHaveTextContent("pending");
  });
});

describe("LoupeSidebar — verdict", () => {
  it("offers Drop for a kept exposure and calls onDropToggle", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} />);
    const toggle = screen.getByTestId("loupe-drop-toggle");
    expect(toggle).toHaveTextContent("Drop");
    fireEvent.click(toggle);
    expect(props.onDropToggle).toHaveBeenCalledTimes(1);
  });

  it("offers Restore for a dropped exposure", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} exposure={exposure({ status: "rejected" })} />);
    expect(screen.getByTestId("loupe-drop-toggle")).toHaveTextContent("Restore");
  });

  // T-4: the kept verdict dot is SAGE (the success status token), not the
  // terracotta interaction accent. DESIGN.md status block pins success = sage.
  it("paints the kept-dot with the sage success token (not the accent)", () => {
    render(<LoupeSidebar {...defaultProps()} />);
    const dot = screen.getByTestId("loupe-kept-dot");
    expect(dot.className).toContain("bg-success");
    expect(dot.className).not.toContain("bg-accent");
    expect(dot.className).not.toContain("bg-print-accent");
  });

  // T-5 boundary: a dropped exposure's dot uses the terracotta accent.
  it("paints the dropped-dot with the print accent", () => {
    render(
      <LoupeSidebar {...defaultProps()} exposure={exposure({ status: "rejected" })} />,
    );
    const dot = screen.getByTestId("loupe-kept-dot");
    expect(dot.className).toContain("bg-print-accent");
    expect(dot.className).not.toContain("bg-success");
  });

  it("displays an existing rejection reason for a dropped exposure", () => {
    const props = defaultProps();
    render(
      <LoupeSidebar
        {...props}
        exposure={exposure({
          status: "rejected",
          tags: [
            { id: 5, key: "rejection_reason", value: "Beam flare", source: "manual" },
          ],
        })}
      />,
    );
    expect(screen.getByTestId("loupe-verdict")).toHaveTextContent("Beam flare");
  });
});

describe("LoupeSidebar — representative", () => {
  it("offers Set as representative and calls onSetRepresentative", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} />);
    const btn = screen.getByTestId("loupe-set-representative");
    fireEvent.click(btn);
    expect(props.onSetRepresentative).toHaveBeenCalledTimes(1);
  });

  it("marks an already-representative exposure and hides the button", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} exposure={exposure({ selected: true })} />);
    expect(screen.getByTestId("loupe-rep")).toHaveTextContent(/Representative/);
    expect(screen.queryByTestId("loupe-set-representative")).not.toBeInTheDocument();
  });
});

describe("LoupeSidebar — sample tags", () => {
  it("renders tags as read-only chips with no remove control", () => {
    const props = defaultProps();
    render(
      <LoupeSidebar
        {...props}
        sample={sample({
          tags: [{ id: 9, key: "lipid", value: "DOPE", source: "manual" }],
        })}
      />,
    );
    const tags = screen.getByTestId("loupe-tags");
    expect(tags).toHaveTextContent("lipid: DOPE");
    // Read-only: SampleMetadataCard's remove button is aria-labelled
    // "Remove <key> tag" — the loupe must not render it (editing is #159).
    expect(screen.queryByLabelText("Remove lipid tag")).not.toBeInTheDocument();
  });

  it("shows an empty-state hint when the sample has no tags", () => {
    render(<LoupeSidebar {...defaultProps()} />);
    expect(screen.getByTestId("loupe-tags")).toHaveTextContent("No tags yet");
  });
});
