import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SegmentedControl } from "../../src/components/ui/SegmentedControl";

type Mode = "log" | "linear";
const TWO = [
  { value: "log" as Mode, label: "log q", testId: "scale-log" },
  { value: "linear" as Mode, label: "linear q", testId: "scale-linear" },
];

type GMode = "bySample" | "byPhase" | "distinct";
const THREE = [
  { value: "bySample" as GMode, label: "By sample" },
  { value: "byPhase" as GMode, label: "By phase" },
  { value: "distinct" as GMode, label: "Distinct" },
];

describe("SegmentedControl — semantics (group)", () => {
  it("renders one button per option with the option label as its name", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    expect(screen.getByRole("button", { name: "log q" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "linear q" })).toBeInTheDocument();
  });

  it("container carries the required aria-label and role=group (default)", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    const root = screen.getByTestId("scale-toggle");
    expect(root).toHaveAttribute("role", "group");
    expect(root).toHaveAttribute("aria-label", "q-axis scale");
  });

  it("group role drives aria-pressed; active=true, others=false", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("button", { name: "log q" })).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByRole("button", { name: "linear q" })).toHaveAttribute("aria-pressed", "false");
  });

  it("clicking an unselected segment fires onChange once with its value", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={onChange}
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: "linear q" }));
    expect(onChange).toHaveBeenCalledTimes(1);
    expect(onChange).toHaveBeenCalledWith("linear");
  });

  it("clicking the already-active segment re-fires onChange with its value", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={onChange}
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: "log q" }));
    expect(onChange).toHaveBeenCalledWith("log");
  });
});

describe("SegmentedControl — data contract (E2E selectors)", () => {
  it("each segment exposes data-value and data-active reflecting value", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="linear"
        onChange={() => {}}
      />,
    );
    const log = screen.getByRole("button", { name: "log q" });
    const lin = screen.getByRole("button", { name: "linear q" });
    expect(log).toHaveAttribute("data-value", "log");
    expect(log).toHaveAttribute("data-active", "false");
    expect(lin).toHaveAttribute("data-value", "linear");
    expect(lin).toHaveAttribute("data-active", "true");
  });

  it("container reflects active value via data-state", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="linear"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    expect(screen.getByTestId("scale-toggle")).toHaveAttribute("data-state", "linear");
  });

  it("applies per-option testId to the segment button", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
      />,
    );
    expect(screen.getByTestId("scale-log")).toBeInTheDocument();
    expect(screen.getByTestId("scale-linear")).toBeInTheDocument();
  });

});

describe("SegmentedControl — radiogroup semantics", () => {
  it("radiogroup role makes segments role=radio with aria-checked", () => {
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("radiogroup", { name: "Trace grouping mode" })).toBeInTheDocument();
    expect(screen.getByRole("radio", { name: "By phase" })).toHaveAttribute("aria-checked", "true");
    expect(screen.getByRole("radio", { name: "By sample" })).toHaveAttribute("aria-checked", "false");
  });

  it("roving tabindex: only the active radio is tabbable", () => {
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("radio", { name: "By phase" })).toHaveAttribute("tabindex", "0");
    expect(screen.getByRole("radio", { name: "By sample" })).toHaveAttribute("tabindex", "-1");
    expect(screen.getByRole("radio", { name: "Distinct" })).toHaveAttribute("tabindex", "-1");
  });

  it("ArrowRight on the active radio selects the next option (onChange)", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={onChange}
      />,
    );
    fireEvent.keyDown(screen.getByRole("radio", { name: "By phase" }), { key: "ArrowRight" });
    expect(onChange).toHaveBeenCalledWith("distinct");
  });

  it("ArrowRight moves DOM focus to the newly-selected radio (focus follows selection)", () => {
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="byPhase"
        onChange={() => {}}
      />,
    );
    fireEvent.keyDown(screen.getByRole("radio", { name: "By phase" }), { key: "ArrowRight" });
    expect(document.activeElement).toBe(screen.getByRole("radio", { name: "Distinct" }));
  });

  it("ArrowLeft wraps from the first option to the last", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<GMode>
        aria-label="Trace grouping mode"
        role="radiogroup"
        variant="plain"
        options={THREE}
        value="bySample"
        onChange={onChange}
      />,
    );
    fireEvent.keyDown(screen.getByRole("radio", { name: "By sample" }), { key: "ArrowLeft" });
    expect(onChange).toHaveBeenCalledWith("distinct");
  });

  it("arrow keys do NOT move selection for role=group", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={onChange}
      />,
    );
    fireEvent.keyDown(screen.getByRole("button", { name: "log q" }), { key: "ArrowRight" });
    expect(onChange).not.toHaveBeenCalled();
  });
});

describe("SegmentedControl — disabled segment + variants/size as data attrs", () => {
  it("a disabled segment is not clickable and is marked disabled", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={[TWO[0], { ...TWO[1], disabled: true, title: "Open a sample" }]}
        value="log"
        onChange={onChange}
      />,
    );
    const lin = screen.getByRole("button", { name: "linear q" });
    expect(lin).toBeDisabled();
    expect(lin).toHaveAttribute("title", "Open a sample");
    fireEvent.click(lin);
    expect(onChange).not.toHaveBeenCalled();
  });

  it("reflects variant + size as data-variant / data-size on the container", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        variant="plain"
        size="md"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    const root = screen.getByTestId("scale-toggle");
    expect(root).toHaveAttribute("data-variant", "plain");
    expect(root).toHaveAttribute("data-size", "md");
  });

  it("defaults to data-variant=bordered data-size=sm", () => {
    render(
      <SegmentedControl<Mode>
        aria-label="q-axis scale"
        options={TWO}
        value="log"
        onChange={() => {}}
        testId="scale-toggle"
      />,
    );
    const root = screen.getByTestId("scale-toggle");
    expect(root).toHaveAttribute("data-variant", "bordered");
    expect(root).toHaveAttribute("data-size", "sm");
  });
});
