import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { TopBar } from "../../src/print/ui/TopBar";

describe("<TopBar>", () => {
  it("renders the wordmark, children, and rightSlot nodes", () => {
    render(
      <TopBar
        wordmark={<span data-testid="wordmark">Himalaya</span>}
        rightSlot={<span data-testid="right">Account</span>}
      >
        <span data-testid="children">tabs</span>
      </TopBar>,
    );
    expect(screen.getByTestId("wordmark")).toBeInTheDocument();
    expect(screen.getByTestId("children")).toBeInTheDocument();
    expect(screen.getByTestId("right")).toBeInTheDocument();
  });

  it("has data-testid=topbar", () => {
    render(<TopBar />);
    expect(screen.getByTestId("topbar")).toBeInTheDocument();
  });

  it("renders a caller-supplied data-testid on the outer header", () => {
    render(<TopBar data-testid="corpus-topbar" />);
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
    expect(screen.queryByTestId("topbar")).toBeNull();
  });

  it("exposes the banner landmark role", () => {
    render(<TopBar />);
    expect(screen.getByRole("banner")).toBeInTheDocument();
  });
});
