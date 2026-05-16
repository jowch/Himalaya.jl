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
      forksList={[{ id: 1, title: "f1", href: "/compare/all/1" }]}
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
      forksList={[]} onCopyLink={onCopyLink} onDelete={() => {}}
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
      forksList={[]} onCopyLink={() => {}} onDelete={() => {}}
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
      forksList={[]} onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={() => {}} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    expect(screen.getByText("Discard changes")).toBeInTheDocument();
  });
});

describe("CompareToolbar — forks list (Compare UX C-17)", () => {
  it("renders forks list inside ⋯ more dropdown when provided", () => {
    render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksList={[
        { id: 9, title: "fork-9", href: "/compare/all/9" },
        { id: 10, title: "fork-10", href: "/compare/all/10" },
      ]}
      onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    expect(screen.getByText("fork-9")).toBeInTheDocument();
    expect(screen.getByText("fork-10")).toBeInTheDocument();
    // Each fork is a link to its review path.
    expect(screen.getByText("fork-9").closest("a"))
      .toHaveAttribute("href", "/compare/all/9");
  });

  it("shows the fork count in the forks section header", () => {
    render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksList={[
        { id: 9, title: "fork-9", href: "/compare/all/9" },
        { id: 10, title: "fork-10", href: "/compare/all/10" },
      ]}
      onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    expect(screen.getByTestId("compare-toolbar-forks").textContent)
      .toMatch(/Forks \(2\)/);
  });

  it("shows an empty state when there are no forks", () => {
    render(wrap(<CompareToolbar
      groupingControl={null} annotationControl={null}
      forksList={[]} onCopyLink={() => {}} onDelete={() => {}}
      onDiscardChanges={null} onFork={() => {}}
      exportControl={null} saveControl={null}
    />));
    fireEvent.click(screen.getByTestId("compare-toolbar-more"));
    const forks = screen.getByTestId("compare-toolbar-forks");
    expect(forks.textContent).toMatch(/Forks \(0\)/);
    expect(forks.textContent?.toLowerCase()).toContain("no forks yet");
  });
});
