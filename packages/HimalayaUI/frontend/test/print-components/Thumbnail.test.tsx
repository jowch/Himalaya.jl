import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Thumbnail } from "../../src/print/components/Thumbnail";

describe("<Thumbnail> root element", () => {
  it("renders [data-testid='thumbnail'] as a <button>", () => {
    render(<Thumbnail src={null} />);
    const btn = screen.getByTestId("thumbnail");
    expect(btn.tagName).toBe("BUTTON");
  });

  it("defaults data-size to 'sm'", () => {
    render(<Thumbnail src={null} />);
    expect(screen.getByTestId("thumbnail")).toHaveAttribute("data-size", "sm");
  });

  it("opts into the inset focus ring so the gallery's overflow can't clip it (UI-RINGCLIP)", () => {
    render(<Thumbnail src={null} />);
    expect(screen.getByTestId("thumbnail")).toHaveAttribute(
      "data-focus-ring",
      "inset",
    );
  });

  it("sets data-size='lg' when size='lg'", () => {
    render(<Thumbnail src={null} size="lg" />);
    expect(screen.getByTestId("thumbnail")).toHaveAttribute("data-size", "lg");
  });

  it("Shift+Enter activates the thumbnail (opens the loupe at this frame) — SA-THUMBKEY", () => {
    const onDoubleClick = vi.fn();
    render(<Thumbnail src={null} onClick={() => {}} onDoubleClick={onDoubleClick} />);
    fireEvent.keyDown(screen.getByTestId("thumbnail"), { key: "Enter", shiftKey: true });
    expect(onDoubleClick).toHaveBeenCalledOnce();
  });

  it("plain Enter does NOT activate the loupe (it stays the native select-toggle)", () => {
    const onDoubleClick = vi.fn();
    render(<Thumbnail src={null} onClick={() => {}} onDoubleClick={onDoubleClick} />);
    fireEvent.keyDown(screen.getByTestId("thumbnail"), { key: "Enter" });
    expect(onDoubleClick).not.toHaveBeenCalled();
  });
});

describe("<Thumbnail> frameNo", () => {
  it("renders [data-role='thumb-fno'] with the frame number text when frameNo is provided", () => {
    render(<Thumbnail src={null} frameNo="65" />);
    const fno = screen.getByText("65");
    expect(fno).toHaveAttribute("data-role", "thumb-fno");
  });

  it("renders frameNo when provided as a number", () => {
    render(<Thumbnail src={null} frameNo={65} />);
    const fno = screen.getByText("65");
    expect(fno).toHaveAttribute("data-role", "thumb-fno");
  });

  it("does NOT render [data-role='thumb-fno'] when frameNo is omitted", () => {
    const { container } = render(<Thumbnail src={null} />);
    expect(container.querySelector("[data-role='thumb-fno']")).not.toBeInTheDocument();
  });

  it("renders [data-role='thumb-fno'] with text '0' when frameNo={0} (valid first frame)", () => {
    render(<Thumbnail src={null} frameNo={0} />);
    const fno = screen.getByText("0");
    expect(fno).toHaveAttribute("data-role", "thumb-fno");
  });

  it("pairs the frame number with a legibility scrim over busy frames (SA-RESCORE3 F10)", () => {
    const { container } = render(<Thumbnail src={null} frameNo={65} />);
    expect(container.querySelector("[data-role='thumb-fno-scrim']")).toBeInTheDocument();
  });

  it("renders no scrim when there is no frame number", () => {
    const { container } = render(<Thumbnail src={null} />);
    expect(container.querySelector("[data-role='thumb-fno-scrim']")).not.toBeInTheDocument();
  });
});

describe("<Thumbnail> marker meaning (SA-RESCORE3 F8)", () => {
  it("the representative marker carries a hover title naming its meaning (not aria-label only)", () => {
    const { container } = render(<Thumbnail src={null} representative />);
    const rep = container.querySelector("[data-role='thumb-rep']");
    expect(rep).toHaveAttribute("title", "Representative exposure");
  });

});

describe("<Thumbnail> rejected dimming", () => {
  it("marks the image wrapper [data-dimmed='true'] when rejected=true", () => {
    const { container } = render(<Thumbnail src={null} rejected />);
    expect(container.querySelector("[data-dimmed='true']")).toBeInTheDocument();
  });

  it("marks the image wrapper [data-dimmed='false'] when not rejected", () => {
    const { container } = render(<Thumbnail src={null} />);
    expect(container.querySelector("[data-dimmed='false']")).toBeInTheDocument();
  });
});

describe("<Thumbnail> representative marker", () => {
  it("renders [data-role='thumb-rep'] when representative=true", () => {
    const { container } = render(<Thumbnail src={null} representative />);
    expect(container.querySelector("[data-role='thumb-rep']")).toBeInTheDocument();
  });

  it("does NOT render [data-role='thumb-rep'] by default", () => {
    const { container } = render(<Thumbnail src={null} />);
    expect(container.querySelector("[data-role='thumb-rep']")).not.toBeInTheDocument();
  });

  it("does NOT render [data-role='thumb-rep'] when representative=false", () => {
    const { container } = render(<Thumbnail src={null} representative={false} />);
    expect(container.querySelector("[data-role='thumb-rep']")).not.toBeInTheDocument();
  });
});

describe("<Thumbnail> kept frame (the default, not-dropped)", () => {
  it("a kept frame stays data-state 'normal' with no status marker", () => {
    const { container } = render(<Thumbnail src={null} frameNo={65} />);
    expect(screen.getByTestId("thumbnail")).toHaveAttribute("data-state", "normal");
    expect(container.querySelector("[data-role='thumb-kept']")).not.toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Frame 65" })).toBeInTheDocument();
  });
});

describe("<Thumbnail> reject overlay", () => {
  it("renders [data-testid='reject-overlay'] when rejected=true", () => {
    render(<Thumbnail src={null} rejected />);
    expect(screen.getByTestId("reject-overlay")).toBeInTheDocument();
  });

  it("does NOT render [data-testid='reject-overlay'] by default", () => {
    render(<Thumbnail src={null} />);
    expect(screen.queryByTestId("reject-overlay")).not.toBeInTheDocument();
  });

  it("does NOT render [data-testid='reject-overlay'] when rejected=false", () => {
    render(<Thumbnail src={null} rejected={false} />);
    expect(screen.queryByTestId("reject-overlay")).not.toBeInTheDocument();
  });
});

describe("<Thumbnail> selected state", () => {
  it("includes 'selected' in data-state when selected=true", () => {
    render(<Thumbnail src={null} selected />);
    const state = screen.getByTestId("thumbnail").getAttribute("data-state") ?? "";
    expect(state.split(" ")).toContain("selected");
  });

  it("does NOT include 'selected' in data-state by default", () => {
    render(<Thumbnail src={null} />);
    const state = screen.getByTestId("thumbnail").getAttribute("data-state") ?? "";
    expect(state.split(" ")).not.toContain("selected");
  });

  it("includes 'normal' in data-state when no state flags are set", () => {
    render(<Thumbnail src={null} />);
    const state = screen.getByTestId("thumbnail").getAttribute("data-state") ?? "";
    expect(state.split(" ")).toContain("normal");
  });
});

describe("<Thumbnail> onClick", () => {
  it("fires onClick when the button is clicked", () => {
    const onClick = vi.fn();
    render(<Thumbnail src={null} onClick={onClick} />);
    fireEvent.click(screen.getByTestId("thumbnail"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("does not throw when clicked without an onClick handler", () => {
    render(<Thumbnail src={null} />);
    expect(() => fireEvent.click(screen.getByTestId("thumbnail"))).not.toThrow();
  });
});

describe("<Thumbnail> combined states", () => {
  it("representative + selected + rejected renders all three markers simultaneously", () => {
    const { container } = render(
      <Thumbnail src={null} representative selected rejected frameNo="65" />
    );
    // representative marker
    expect(container.querySelector("[data-role='thumb-rep']")).toBeInTheDocument();
    // reject overlay
    expect(screen.getByTestId("reject-overlay")).toBeInTheDocument();
    // frame number
    expect(container.querySelector("[data-role='thumb-fno']")).toBeInTheDocument();
    // data-state contains all three tokens
    const state = screen.getByTestId("thumbnail").getAttribute("data-state") ?? "";
    const tokens = state.split(" ");
    expect(tokens).toContain("selected");
    expect(tokens).toContain("rejected");
    expect(tokens).toContain("representative");
  });
});

describe("<Thumbnail> accessible name (WCAG 4.1.2)", () => {
  it("names the button 'Frame N' when frameNo is set", () => {
    render(<Thumbnail src={null} frameNo={65} />);
    expect(screen.getByRole("button", { name: "Frame 65" })).toBeInTheDocument();
  });

  it("names the button 'Exposure' when frameNo is omitted", () => {
    render(<Thumbnail src={null} />);
    expect(screen.getByRole("button", { name: "Exposure" })).toBeInTheDocument();
  });

  it("appends ', representative' to the name", () => {
    render(<Thumbnail src={null} frameNo={65} representative />);
    expect(
      screen.getByRole("button", { name: "Frame 65, representative" }),
    ).toBeInTheDocument();
  });

  it("appends ', dropped' to the name when rejected", () => {
    render(<Thumbnail src={null} frameNo={65} rejected />);
    expect(
      screen.getByRole("button", { name: "Frame 65, dropped" }),
    ).toBeInTheDocument();
  });

  it("appends both suffixes in order (representative then dropped)", () => {
    render(<Thumbnail src={null} frameNo={65} representative rejected />);
    expect(
      screen.getByRole("button", { name: "Frame 65, representative, dropped" }),
    ).toBeInTheDocument();
  });
});

describe("<Thumbnail> aria-pressed (selection toggle)", () => {
  it("reflects selected=true via aria-pressed when onClick is provided", () => {
    render(<Thumbnail src={null} frameNo={65} selected onClick={() => {}} />);
    expect(screen.getByRole("button", { name: /Frame 65/ })).toHaveAttribute(
      "aria-pressed",
      "true",
    );
  });

  it("reflects selected=false via aria-pressed when onClick is provided", () => {
    render(<Thumbnail src={null} frameNo={65} onClick={() => {}} />);
    expect(screen.getByRole("button", { name: /Frame 65/ })).toHaveAttribute(
      "aria-pressed",
      "false",
    );
  });

  it("omits aria-pressed entirely when no onClick (not a toggle)", () => {
    render(<Thumbnail src={null} frameNo={65} selected />);
    expect(screen.getByRole("button", { name: /Frame 65/ })).not.toHaveAttribute(
      "aria-pressed",
    );
  });
});

describe("<Thumbnail> DetectorImage integration", () => {
  it("renders the image placeholder when src is null", () => {
    render(<Thumbnail src={null} />);
    expect(screen.getByTestId("detector-image-placeholder")).toBeInTheDocument();
  });
});

describe("<Thumbnail> xs size (dense exposure strip)", () => {
  it("renders the xs (dense exposure-strip) size", () => {
    render(<Thumbnail src={null} size="xs" />);
    expect(screen.getByTestId("thumbnail").dataset.size).toBe("xs");
  });

  it("sets data-size='xs' when size='xs'", () => {
    render(<Thumbnail src={null} size="xs" />);
    expect(screen.getByTestId("thumbnail")).toHaveAttribute("data-size", "xs");
  });
});
