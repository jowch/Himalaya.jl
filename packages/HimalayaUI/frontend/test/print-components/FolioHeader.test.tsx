import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { FolioHeader } from "../../src/print/components/FolioHeader";

describe("<FolioHeader>", () => {
  it("renders the root data-testid=folio-header", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={5} countLabel="series in the folio" />);
    expect(screen.getByTestId("folio-header")).toBeInTheDocument();
  });

  it("renders the kicker text", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={5} countLabel="series in the folio" />);
    expect(screen.getByText("Folio")).toBeInTheDocument();
  });

  it("renders the title inside an h1", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={5} countLabel="series in the folio" />);
    const heading = screen.getByRole("heading", { level: 1 });
    expect(heading).toHaveTextContent("Saved series");
  });

  it("renders the count number", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={5} countLabel="series in the folio" />);
    expect(screen.getByText("5")).toBeInTheDocument();
  });

  it("renders the count label", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={5} countLabel="series in the folio" />);
    expect(screen.getByText("series in the folio")).toBeInTheDocument();
  });
});

describe("<FolioHeader> nullable count (FOL-HONEST-DERIVED)", () => {
  it("count={null} renders the em-dash placeholder, not '0'", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={null} countLabel="series in the folio" />);
    const countEl = screen.getByTestId("folio-count");
    expect(countEl).toHaveTextContent("—");
    expect(countEl).not.toHaveTextContent("0");
  });

  it("count={null} exposes the placeholder via the house role=img + label pattern (reliably announced, unlike a bare aria-label on a roleless div)", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={null} countLabel="series in the folio" />);
    const countEl = screen.getByTestId("folio-count");
    expect(countEl).toHaveAttribute("role", "img");
    expect(countEl).toHaveAttribute("aria-label", "series count unavailable");
  });

  it("count={0} still renders '0' (a genuine empty folio stays honest)", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={0} countLabel="series in the folio" />);
    const countEl = screen.getByTestId("folio-count");
    expect(countEl).toHaveTextContent("0");
    expect(countEl).not.toHaveTextContent("—");
  });

  it("count={0} does NOT render the aria-label (0 is a real value)", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={0} countLabel="series in the folio" />);
    expect(screen.getByTestId("folio-count")).not.toHaveAttribute("aria-label");
  });
});

describe("<FolioHeader> subtitle", () => {
  it("renders the subtitle when passed", () => {
    render(
      <FolioHeader
        kicker="Folio"
        title="Saved series"
        subtitle="A wall of saved series across every beamtime."
        count={5}
        countLabel="series in the folio"
      />,
    );
    expect(screen.getByText("A wall of saved series across every beamtime.")).toBeInTheDocument();
  });

  it("does NOT render the subtitle text when omitted", () => {
    render(<FolioHeader kicker="Folio" title="Saved series" count={5} countLabel="series in the folio" />);
    expect(screen.queryByText("A wall of saved series across every beamtime.")).not.toBeInTheDocument();
  });
});
