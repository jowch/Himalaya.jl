import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { useState } from "react";
import { RepresentativeBox } from "../../src/print/components/RepresentativeBox";

describe("<RepresentativeBox> not representative", () => {
  it('does NOT show "Representative for indexing" header', () => {
    render(<RepresentativeBox isRepresentative={false} />);
    expect(
      screen.queryByText("Representative for indexing")
    ).not.toBeInTheDocument();
  });

  it("shows the default body text", () => {
    render(<RepresentativeBox isRepresentative={false} />);
    expect(
      screen.getByText(
        "One frame per sample carries forward to the Index stage. Pick the cleanest, strongest frame."
      )
    ).toBeInTheDocument();
  });

  it('button is labeled /set as representative/i', () => {
    render(<RepresentativeBox isRepresentative={false} onSetRepresentative={() => {}} />);
    expect(
      screen.getByRole("button", { name: /set as representative/i })
    ).toBeInTheDocument();
  });

  it("clicking the button calls onSetRepresentative", () => {
    const onSetRepresentative = vi.fn();
    render(
      <RepresentativeBox
        isRepresentative={false}
        onSetRepresentative={onSetRepresentative}
      />
    );
    fireEvent.click(screen.getByRole("button", { name: /set as representative/i }));
    expect(onSetRepresentative).toHaveBeenCalledTimes(1);
  });

  it("Card does NOT have data-selected when not representative", () => {
    const { container } = render(<RepresentativeBox isRepresentative={false} />);
    expect(container.querySelector("[data-selected='true']")).not.toBeInTheDocument();
  });
});

describe("<RepresentativeBox> is representative", () => {
  it('shows "Representative for indexing" header', () => {
    render(<RepresentativeBox isRepresentative={true} />);
    expect(
      screen.getByText("Representative for indexing")
    ).toBeInTheDocument();
  });

  it("Card has data-selected=true when isRepresentative=true", () => {
    const { container } = render(<RepresentativeBox isRepresentative={true} />);
    expect(container.querySelector("[data-selected='true']")).toBeInTheDocument();
  });

  it("shows default body text even when representative", () => {
    render(<RepresentativeBox isRepresentative={true} />);
    expect(
      screen.getByText(
        "One frame per sample carries forward to the Index stage. Pick the cleanest, strongest frame."
      )
    ).toBeInTheDocument();
  });
});

describe("<RepresentativeBox> children override", () => {
  it("replaces the default body text when children provided", () => {
    render(
      <RepresentativeBox isRepresentative={false}>
        Custom body text here.
      </RepresentativeBox>
    );
    expect(screen.getByText("Custom body text here.")).toBeInTheDocument();
    expect(
      screen.queryByText(
        "One frame per sample carries forward to the Index stage. Pick the cleanest, strongest frame."
      )
    ).not.toBeInTheDocument();
  });
});

describe("<RepresentativeBox> no onSetRepresentative", () => {
  it("does not throw when clicked without onSetRepresentative", () => {
    render(<RepresentativeBox isRepresentative={false} />);
    expect(() =>
      fireEvent.click(screen.getByRole("button", { name: /set as representative/i }))
    ).not.toThrow();
  });
});

describe("<RepresentativeBox> dropped-representative warning (LO-REPDROP)", () => {
  const WARNING_COPY =
    "Warning. The representative frame is dropped. It still carries to the Index stage. Restore it or set another frame as representative.";

  it("renders the warning when representativeDropped, leading with the severity word", () => {
    render(
      <RepresentativeBox isRepresentative={false} representativeDropped />
    );
    const warning = screen.getByTestId("rep-dropped-warning");
    expect(warning).toHaveTextContent(WARNING_COPY);
    expect(warning.textContent?.startsWith("Warning.")).toBe(true);
  });

  it("renders the warning even on the representative frame itself", () => {
    render(<RepresentativeBox isRepresentative={true} representativeDropped />);
    expect(screen.getByTestId("rep-dropped-warning")).toBeInTheDocument();
  });

  it("no warning when the prop is absent", () => {
    render(<RepresentativeBox isRepresentative={false} />);
    expect(screen.queryByTestId("rep-dropped-warning")).not.toBeInTheDocument();
  });

  it("no warning when representativeDropped is false", () => {
    render(
      <RepresentativeBox isRepresentative={false} representativeDropped={false} />
    );
    expect(screen.queryByTestId("rep-dropped-warning")).not.toBeInTheDocument();
  });
});

describe("<RepresentativeBox> focus retention on set (button unmounts on success)", () => {
  // Mirrors production: selectExposureMutator's optimistic onMutate flips
  // `selected` instantly, so the clicked button unmounts in the same
  // interaction (controls-don't-lie) and would otherwise dump focus on body.
  function Harness(): JSX.Element {
    const [isRep, setIsRep] = useState(false);
    return (
      <RepresentativeBox
        isRepresentative={isRep}
        onSetRepresentative={() => setIsRep(true)}
      />
    );
  }

  it("moves focus to the box when the activated button unmounts", () => {
    render(<Harness />);
    const btn = screen.getByRole("button", { name: /set as representative/i });
    btn.focus();
    fireEvent.click(btn);
    // The button is gone…
    expect(
      screen.queryByRole("button", { name: /set as representative/i })
    ).not.toBeInTheDocument();
    // …and focus did NOT fall back to document.body.
    expect(screen.getByTestId("representative-box")).toHaveFocus();
  });

  it("a flip arriving WITHOUT interaction inside the box (e.g. SSE from another user) does not steal focus", () => {
    const { rerender } = render(
      <RepresentativeBox isRepresentative={false} onSetRepresentative={() => {}} />
    );
    rerender(
      <RepresentativeBox isRepresentative={true} onSetRepresentative={() => {}} />
    );
    expect(screen.getByTestId("representative-box")).not.toHaveFocus();
    expect(document.body).toHaveFocus();
  });
});

describe("<RepresentativeBox> controls-don't-lie (LO-REPLIES)", () => {
  it("OMITS the set-as-representative button when this frame IS the representative", () => {
    render(
      <RepresentativeBox isRepresentative={true} onSetRepresentative={() => {}} />
    );
    expect(
      screen.queryByRole("button", { name: /set as representative/i })
    ).not.toBeInTheDocument();
  });

  it("keeps the button when this frame is NOT the representative", () => {
    render(
      <RepresentativeBox isRepresentative={false} onSetRepresentative={() => {}} />
    );
    expect(
      screen.getByRole("button", { name: /set as representative/i })
    ).toBeInTheDocument();
  });
});
