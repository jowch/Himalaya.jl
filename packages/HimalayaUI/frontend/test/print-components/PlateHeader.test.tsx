import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PlateHeader } from "../../src/print/components/PlateHeader";

describe("<PlateHeader>", () => {
  it("renders the title text", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.getByText("Lipid 1‑2 + LL37 1:0.5")).toBeInTheDocument();
  });

  it("renders the title as an h2 by default", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.getByRole("heading", { level: 2 })).toBeInTheDocument();
  });

  it("renders the title as an h1 when as='h1'", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" as="h1" />);
    expect(screen.getByRole("heading", { level: 1 })).toBeInTheDocument();
  });

  it("renders the title as an h3 when as='h3'", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" as="h3" />);
    expect(screen.getByRole("heading", { level: 3 })).toBeInTheDocument();
  });

  it("renders kicker text when provided", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" kicker="Integration" />);
    expect(screen.getByText("Integration")).toBeInTheDocument();
  });

  it("does not render a kicker when omitted", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.queryByText("Integration")).not.toBeInTheDocument();
  });

  it("renders subtitle text when provided", () => {
    render(
      <PlateHeader
        title="Lipid 1‑2 + LL37 1:0.5"
        subtitle="smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03"
      />
    );
    expect(
      screen.getByText("smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03")
    ).toBeInTheDocument();
  });

  it("does not render a subtitle when omitted", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" />);
    expect(
      screen.queryByText("smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03")
    ).not.toBeInTheDocument();
  });

  it("renders children (tools slot) when provided", () => {
    render(
      <PlateHeader title="Lipid 1‑2 + LL37 1:0.5">
        <button>Export</button>
      </PlateHeader>
    );
    expect(screen.getByRole("button", { name: "Export" })).toBeInTheDocument();
  });

  it("does not render a tools slot region when children are omitted", () => {
    render(<PlateHeader title="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.queryByRole("button", { name: "Export" })).not.toBeInTheDocument();
  });

  it("uses headingText for the heading and keeps an interactive title out of the heading (BU-NOHEAD)", () => {
    render(
      <PlateHeader
        as="h1"
        headingText="LL37 ratio series"
        title={<input aria-label="Series title" defaultValue="LL37 ratio series" />}
      />
    );
    // The named heading comes from headingText, not the interactive control.
    expect(
      screen.getByRole("heading", { level: 1, name: "LL37 ratio series" })
    ).toBeInTheDocument();
    // The control renders, but NOT as a descendant of the heading.
    expect(screen.getByLabelText("Series title").closest("h1")).toBeNull();
  });
});
