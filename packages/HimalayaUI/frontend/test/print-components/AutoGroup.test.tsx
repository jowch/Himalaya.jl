import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { AutoGroup } from "../../src/print/components/AutoGroup";

describe("<AutoGroup> rendering", () => {
  it("renders the root callout box", () => {
    render(<AutoGroup>grouped six samples</AutoGroup>);
    expect(screen.getByTestId("auto-group")).toBeInTheDocument();
  });

  it("renders the star icon", () => {
    const { container } = render(<AutoGroup>body</AutoGroup>);
    expect(container.querySelector("[data-role='autogroup-star']")).toBeInTheDocument();
  });

  it("renders the body text", () => {
    render(<AutoGroup>read six samples as one series</AutoGroup>);
    expect(screen.getByText("read six samples as one series")).toBeInTheDocument();
  });

  it("renders body with embedded emphasis", () => {
    render(
      <AutoGroup>
        You selected <strong className="text-ink font-semibold">6 samples</strong>.
      </AutoGroup>,
    );
    expect(screen.getByText("6 samples")).toBeInTheDocument();
  });
});

describe("<AutoGroup> variant", () => {
  it("defaults to the summary variant", () => {
    render(<AutoGroup>body</AutoGroup>);
    expect(screen.getByTestId("auto-group")).toHaveAttribute("data-variant", "summary");
  });

  it("reflects the compose variant", () => {
    render(<AutoGroup variant="compose">body</AutoGroup>);
    expect(screen.getByTestId("auto-group")).toHaveAttribute("data-variant", "compose");
  });

  it("reflects the summary variant explicitly", () => {
    render(<AutoGroup variant="summary">body</AutoGroup>);
    expect(screen.getByTestId("auto-group")).toHaveAttribute("data-variant", "summary");
  });
});

describe("<AutoGroup> title", () => {
  it("does not render a title by default", () => {
    render(<AutoGroup>body copy here</AutoGroup>);
    expect(screen.queryByText("Auto-grouped")).not.toBeInTheDocument();
  });

  it("renders the title when passed", () => {
    render(
      <AutoGroup variant="compose" title="Auto-grouped">
        body copy here
      </AutoGroup>,
    );
    expect(screen.getByText("Auto-grouped")).toBeInTheDocument();
  });
});

describe("<AutoGroup> actions", () => {
  it("renders no action buttons when none are passed", () => {
    render(<AutoGroup>body</AutoGroup>);
    expect(screen.queryByRole("button")).not.toBeInTheDocument();
  });

  it("renders a button per action with its label", () => {
    render(
      <AutoGroup actions={[{ label: "Confirm series" }, { label: "Adjust", muted: true }]}>
        body
      </AutoGroup>,
    );
    expect(screen.getByRole("button", { name: "Confirm series" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Adjust" })).toBeInTheDocument();
  });

  it("fires onClick when an action is clicked", () => {
    const onConfirm = vi.fn();
    render(
      <AutoGroup actions={[{ label: "Confirm series", onClick: onConfirm }]}>body</AutoGroup>,
    );
    fireEvent.click(screen.getByRole("button", { name: "Confirm series" }));
    expect(onConfirm).toHaveBeenCalledTimes(1);
  });

  it("renders an action with no onClick as DISABLED (controls-don't-lie)", () => {
    render(<AutoGroup actions={[{ label: "Confirm series" }]}>body</AutoGroup>);
    const btn = screen.getByRole("button", { name: "Confirm series" });
    expect(btn).toBeDisabled();
    expect(btn).toHaveAttribute("aria-disabled", "true");
  });

  it("an action WITH onClick stays enabled even alongside an inert sibling", () => {
    const onAdjust = vi.fn();
    render(
      <AutoGroup
        actions={[
          { label: "Confirm series" },
          { label: "Adjust", muted: true, onClick: onAdjust },
        ]}
      >
        body
      </AutoGroup>,
    );
    expect(screen.getByRole("button", { name: "Confirm series" })).toBeDisabled();
    const adjust = screen.getByRole("button", { name: "Adjust" });
    expect(adjust).not.toBeDisabled();
    fireEvent.click(adjust);
    expect(onAdjust).toHaveBeenCalledTimes(1);
  });
});
