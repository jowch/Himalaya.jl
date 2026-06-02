import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
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
        "One exposure per sample carries forward to the Index stage. Pick the cleanest, strongest frame."
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
        "One exposure per sample carries forward to the Index stage. Pick the cleanest, strongest frame."
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
        "One exposure per sample carries forward to the Index stage. Pick the cleanest, strongest frame."
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
