import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { CompareTitleStrip } from "../src/components/CompareTitleStrip";

const wrap = (ui: JSX.Element) => <MemoryRouter>{ui}</MemoryRouter>;

describe("CompareTitleStrip — Compare UX C-8", () => {
  it("renders title and meta", () => {
    render(wrap(<CompareTitleStrip
      title="Cubic vs Hex"
      description={null}
      memberCount={4}
      authorUsername="alice"
      isCurrentUserAuthor={false}
      lastEventAt={null}
      forkedFromTitle={null}
      forkedFromHref={null}
      onTitleChange={() => {}}
      onDescChange={() => {}}
      now={Date.parse("2026-05-14T12:00:00Z")}
    />));
    expect(screen.getByText("Cubic vs Hex")).toBeInTheDocument();
    expect(screen.getByText(/by alice/)).toBeInTheDocument();
    expect(screen.getByText(/4 traces/)).toBeInTheDocument();
  });

  it("renders 'by you' when current user is the author", () => {
    render(wrap(<CompareTitleStrip
      title="t" description={null}
      memberCount={1}
      authorUsername="alice"
      isCurrentUserAuthor={true}
      lastEventAt={null} forkedFromTitle={null} forkedFromHref={null}
      onTitleChange={() => {}} onDescChange={() => {}}
      now={Date.now()}
    />));
    expect(screen.getByText(/by you/)).toBeInTheDocument();
  });

  it("renders forked-from link when set", () => {
    render(wrap(<CompareTitleStrip
      title="t" description={null}
      memberCount={1} authorUsername={null} isCurrentUserAuthor={false}
      lastEventAt={null}
      forkedFromTitle="Parent"
      forkedFromHref="/compare/all/9"
      onTitleChange={() => {}} onDescChange={() => {}}
      now={Date.now()}
    />));
    const link = screen.getByText("Parent");
    expect(link.closest("a")).toHaveAttribute("href", "/compare/all/9");
  });

  it("commits a title edit", () => {
    const onTitleChange = vi.fn();
    render(wrap(<CompareTitleStrip
      title="old" description={null}
      memberCount={1} authorUsername={null} isCurrentUserAuthor={false}
      lastEventAt={null} forkedFromTitle={null} forkedFromHref={null}
      onTitleChange={onTitleChange}
      onDescChange={() => {}}
      now={Date.now()}
    />));
    fireEvent.click(screen.getByText("old"));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "new" } });
    fireEvent.keyDown(screen.getByRole("textbox"), { key: "Enter" });
    expect(onTitleChange).toHaveBeenCalledWith("new");
  });
});
