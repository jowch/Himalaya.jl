import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
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

  it("carries the summary variant marker", () => {
    render(<AutoGroup>body</AutoGroup>);
    expect(screen.getByTestId("auto-group")).toHaveAttribute("data-variant", "summary");
  });
});
