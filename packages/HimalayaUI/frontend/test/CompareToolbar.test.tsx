import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { CompareToolbar } from "../src/components/CompareToolbar";

const wrap = (ui: JSX.Element) => <MemoryRouter>{ui}</MemoryRouter>;

describe("CompareToolbar — Compare UX C-10", () => {
  it("mounts the children controls", () => {
    render(wrap(<CompareToolbar
      groupingControl={<button data-testid="g">G</button>}
      annotationControl={<button data-testid="a">A</button>}
      forksCount={2}
      onCopyLink={() => {}}
      onDelete={() => {}}
      onDiscardChanges={null}
      onFork={() => {}}
      exportControl={<button data-testid="x">X</button>}
      saveControl={<button data-testid="s">S</button>}
    />));
    expect(screen.getByTestId("g")).toBeInTheDocument();
    expect(screen.getByTestId("a")).toBeInTheDocument();
    expect(screen.getByTestId("x")).toBeInTheDocument();
    expect(screen.getByTestId("s")).toBeInTheDocument();
  });

  it("opens the more menu and triggers actions", () => {
    const onCopyLink = vi.fn();
    render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksCount={0} onCopyLink={onCopyLink} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    fireEvent.click(screen.getByText("Copy link"));
    expect(onCopyLink).toHaveBeenCalled();
  });

  it("includes 'Discard changes' only when dirty", () => {
    const { queryByText, rerender } = render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksCount={0} onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    expect(queryByText("Discard changes")).toBeNull();

    // Menu is still open from the click above; rerender preserves component
    // state, so clicking "more" again would toggle it closed. The newly
    // dirty toolbar should surface "Discard changes" in the open menu.
    rerender(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksCount={0} onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={() => {}} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    expect(screen.getByText("Discard changes")).toBeInTheDocument();
  });
});
