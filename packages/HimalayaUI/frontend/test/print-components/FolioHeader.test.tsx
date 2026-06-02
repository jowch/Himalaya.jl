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
