import { beforeEach, describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { ReactNode } from "react";
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
    signalLevel: 0,
    onDropToggle: vi.fn(),
    onSetRepresentative: vi.fn(),
  };
}

// LoupeSidebar now wires the corpus tag mutations (#207 / #159), so it needs
// a QueryClientProvider in scope. Mock fetch so the mutators don't try to hit
// a real backend during the unrelated suites; tag-specific tests below stub it.
function wrap(node: ReactNode) {
  const client = makeClient();
  return (
    <QueryClientProvider client={client}>{node}</QueryClientProvider>
  );
}

beforeEach(() => {
  vi.restoreAllMocks();
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }),
  );
});

describe("LoupeSidebar — meta-list", () => {
  // R2-M13: the redundant Status row drops here; Kind stays.
  it("shows filename, kind and frame position", () => {
    const exposures = [
      exposure({ id: 1 }), exposure({ id: 2, kind: "averaged", filename: null }),
    ];
    render(wrap(
      <LoupeSidebar
        {...defaultProps()}
        exposure={exposures[1]}
        exposures={exposures}
      />,
    ));
    expect(screen.getByTestId("loupe-meta-filename")).toHaveTextContent("—");
    expect(screen.getByTestId("loupe-meta-kind")).toHaveTextContent("averaged");
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    // The Status row is gone — the verdict card carries the same fact.
    expect(screen.queryByTestId("loupe-meta-status")).not.toBeInTheDocument();
  });
});

describe("LoupeSidebar — signal meter (M-8)", () => {
  it("renders a 5-bar signal meter in the meta list", () => {
    render(wrap(<LoupeSidebar {...defaultProps()} signalLevel={3} />));
    const meter = screen.getByTestId("loupe-meta-signal");
    expect(meter).toBeInTheDocument();
    expect(meter.querySelectorAll("[data-bar]")).toHaveLength(5);
  });

  it("fills exactly signalLevel bars", () => {
    render(wrap(<LoupeSidebar {...defaultProps()} signalLevel={3} />));
    const meter = screen.getByTestId("loupe-meta-signal");
    expect(meter.querySelectorAll('[data-bar="on"]')).toHaveLength(3);
    expect(meter.querySelectorAll('[data-bar="off"]')).toHaveLength(2);
  });

  it("clamps signalLevel into the 0..5 range", () => {
    const { rerender } = render(
      wrap(<LoupeSidebar {...defaultProps()} signalLevel={9} />),
    );
    expect(
      screen.getByTestId("loupe-meta-signal").querySelectorAll('[data-bar="on"]'),
    ).toHaveLength(5);
    rerender(wrap(<LoupeSidebar {...defaultProps()} signalLevel={-2} />));
    expect(
      screen.getByTestId("loupe-meta-signal").querySelectorAll('[data-bar="on"]'),
    ).toHaveLength(0);
  });
});

describe("LoupeSidebar — verdict", () => {
  it("offers Drop for a kept exposure and calls onDropToggle", () => {
    const props = defaultProps();
    render(wrap(<LoupeSidebar {...props} />));
    const toggle = screen.getByTestId("loupe-drop-toggle");
    expect(toggle).toHaveTextContent("Drop");
    fireEvent.click(toggle);
    expect(props.onDropToggle).toHaveBeenCalledTimes(1);
  });

  it("offers Restore for a dropped exposure", () => {
    const props = defaultProps();
    render(wrap(
      <LoupeSidebar {...props} exposure={exposure({ status: "rejected" })} />,
    ));
    expect(screen.getByTestId("loupe-drop-toggle")).toHaveTextContent("Restore");
  });

  // T-4: the kept verdict dot is SAGE (the success status token), not the
  // terracotta interaction accent. DESIGN.md status block pins success = sage.
  it("paints the kept-dot with the sage success token (not the accent)", () => {
    render(wrap(<LoupeSidebar {...defaultProps()} />));
    const dot = screen.getByTestId("loupe-kept-dot");
    expect(dot.className).toContain("bg-success");
    expect(dot.className).not.toContain("bg-accent");
    expect(dot.className).not.toContain("bg-print-accent");
  });

  // T-5 boundary: a dropped exposure's dot uses the terracotta accent.
  it("paints the dropped-dot with the print accent", () => {
    render(wrap(
      <LoupeSidebar {...defaultProps()} exposure={exposure({ status: "rejected" })} />,
    ));
    const dot = screen.getByTestId("loupe-kept-dot");
    expect(dot.className).toContain("bg-print-accent");
    expect(dot.className).not.toContain("bg-success");
  });

  it("displays an existing rejection reason for a dropped exposure", () => {
    const props = defaultProps();
    render(wrap(
      <LoupeSidebar
        {...props}
        exposure={exposure({
          status: "rejected",
          tags: [
            { id: 5, key: "rejection_reason", value: "Beam flare", source: "manual" },
          ],
        })}
      />,
    ));
    expect(screen.getByTestId("loupe-verdict")).toHaveTextContent("Beam flare");
  });
});

describe("LoupeSidebar — representative", () => {
  it("offers Set as representative and calls onSetRepresentative", () => {
    const props = defaultProps();
    render(wrap(<LoupeSidebar {...props} />));
    const btn = screen.getByTestId("loupe-set-representative");
    fireEvent.click(btn);
    expect(props.onSetRepresentative).toHaveBeenCalledTimes(1);
  });

  it("marks an already-representative exposure and hides the button", () => {
    const props = defaultProps();
    render(wrap(
      <LoupeSidebar {...props} exposure={exposure({ selected: true })} />,
    ));
    expect(screen.getByTestId("loupe-rep")).toHaveTextContent(/Representative/);
    expect(screen.queryByTestId("loupe-set-representative")).not.toBeInTheDocument();
  });

  // R2-M12: when active, the rep-box swaps the neutral border for a
  // terracotta-tinged border, mirroring the mockup's `.rep-box.is-rep`.
  it("marks the rep-box with data-is-rep when representative", () => {
    const props = defaultProps();
    render(wrap(
      <LoupeSidebar {...props} exposure={exposure({ selected: true })} />,
    ));
    expect(screen.getByTestId("loupe-rep")).toHaveAttribute("data-is-rep", "true");
  });

  it("omits the data-is-rep attribute when not representative", () => {
    render(wrap(<LoupeSidebar {...defaultProps()} />));
    expect(screen.getByTestId("loupe-rep")).not.toHaveAttribute("data-is-rep");
  });
});

describe("LoupeSidebar — sample tags", () => {
  // L-10: tags render as free-form single tokens (bare value), not
  // `key: value`. The loupe now hosts read+write tag editing (#159 / #207).
  it("renders tags as bare value tokens", () => {
    const props = defaultProps();
    render(wrap(
      <LoupeSidebar
        {...props}
        sample={sample({
          tags: [{ id: 9, key: "lipid", value: "DOPE", source: "manual" }],
        })}
      />,
    ));
    const tags = screen.getByTestId("loupe-tags");
    expect(tags).toHaveTextContent("DOPE");
    expect(tags).not.toHaveTextContent("lipid: DOPE");
  });

  it("renders a remove-button alongside each tag chip (#159)", () => {
    const props = defaultProps();
    render(wrap(
      <LoupeSidebar
        {...props}
        sample={sample({
          tags: [{ id: 9, key: "lipid", value: "DOPE", source: "manual" }],
        })}
      />,
    ));
    // The corpus tag editor mirrors SampleMetadataCard's aria-label convention:
    // "Remove <key> tag" — fall back to the value when the key is empty.
    expect(screen.getByLabelText("Remove lipid tag")).toBeInTheDocument();
  });

  it("shows a `+ tag` invite when the sample has no tags", () => {
    render(wrap(<LoupeSidebar {...defaultProps()} />));
    const tags = screen.getByTestId("loupe-tags");
    expect(tags).toHaveTextContent("+ tag");
  });

  it("opens an inline form when the `+ tag` invite is clicked", () => {
    render(wrap(<LoupeSidebar {...defaultProps()} />));
    fireEvent.click(screen.getByTestId("loupe-tag-add"));
    expect(screen.getByTestId("loupe-tag-form")).toBeInTheDocument();
    expect(screen.getByPlaceholderText("value")).toBeInTheDocument();
  });
});
