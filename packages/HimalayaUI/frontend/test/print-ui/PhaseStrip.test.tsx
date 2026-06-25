import { render, screen, within } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PhaseStrip } from "../../src/print/ui/PhaseStrip";
import type { PhaseSegment } from "../../src/print/ui/PhaseStrip";

function segs(...segments: PhaseSegment[]): PhaseSegment[] {
  return segments;
}

describe("<PhaseStrip> single phase", () => {
  it("renders a solid background (no gradient) and no data-coexist-count", () => {
    render(<PhaseStrip segments={segs({ phase: "Pn3m" })} />);
    const cell = screen.getByTestId("ps-seg");
    const bg = cell.style.background;
    expect(bg).not.toContain("linear-gradient");
    expect(bg.length).toBeGreaterThan(0);
    expect(cell.getAttribute("data-coexist-count")).toBeNull();
    expect(cell.getAttribute("aria-label")).toBe("Pn3m");
    expect(cell.getAttribute("title")).toBe("Pn3m");
  });
});

describe("<PhaseStrip> 2-phase coexistence", () => {
  it("labels both phases, sets data-coexist-count=2, and paints an N-band gradient", () => {
    render(
      <PhaseStrip
        segments={segs({ phase: "Pn3m", coexistWith: ["Lamellar"] })}
      />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.getAttribute("aria-label")).toBe("Pn3m + Lamellar (coexistence)");
    expect(cell.getAttribute("title")).toBe("Pn3m + Lamellar (coexistence)");
    expect(cell.getAttribute("data-coexist-count")).toBe("2");

    const bg = cell.style.background;
    expect(bg.startsWith("linear-gradient(100deg")).toBe(true);
    // Equal hard bands: the 1/2 boundary lands at 50.00%.
    expect(bg).toContain("50.00%");
    // Two distinct color stops (the dominant + the coexisting phase).
    const colorStops = bg.match(/oklch|rgb|hsl|#|var\(/g) ?? [];
    expect(colorStops.length).toBeGreaterThanOrEqual(2);
  });
});

describe("<PhaseStrip> 3-phase coexistence", () => {
  it("labels all three phases, sets data-coexist-count=3, and paints a gradient", () => {
    render(
      <PhaseStrip
        segments={segs({
          phase: "Pn3m",
          coexistWith: ["Lamellar", "Hexagonal"],
        })}
      />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.getAttribute("aria-label")).toBe(
      "Pn3m + Lamellar + Hexagonal (coexistence)",
    );
    expect(cell.getAttribute("title")).toBe(
      "Pn3m + Lamellar + Hexagonal (coexistence)",
    );
    expect(cell.getAttribute("data-coexist-count")).toBe("3");
    expect(cell.style.background.startsWith("linear-gradient(100deg")).toBe(
      true,
    );
  });
});

describe("<PhaseStrip> assistive-tech exposure", () => {
  it("exposes every segment as role=img with the right accessible names", () => {
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Pn3m" },
          { phase: null },
          { phase: null, state: "form_factor" },
          { phase: null, state: "null" },
          { phase: "Pn3m", coexistWith: ["Lamellar"] },
        )}
      />,
    );
    const cells = screen.getAllByRole("img");
    expect(cells).toHaveLength(5);
    expect(screen.getByRole("img", { name: "Pn3m" })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "Unindexed" })).toBeInTheDocument();
    expect(
      screen.getByRole("img", { name: "Form factor (no Bragg peaks)" }),
    ).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "No phase" })).toBeInTheDocument();
    expect(
      screen.getByRole("img", { name: "Pn3m + Lamellar (coexistence)" }),
    ).toBeInTheDocument();
  });

  it("speaks the transition relation: arrow stays decorative, 'to' is accessible", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m" }, { phase: "Lamellar" })} />,
    );
    const cap = screen.getByTestId("ps-cap");
    // The visual arrow is hidden from AT...
    const arrow = cap.querySelector('[aria-hidden="true"]');
    expect(arrow).not.toBeNull();
    expect(arrow).toHaveTextContent("→");
    // ...and the relation word "to" is exposed between the phase names.
    const to = within(cap).getByText("to");
    expect(to).not.toHaveAttribute("aria-hidden");
  });
});

describe("<PhaseStrip> caption truthfulness (SC-PREVCLAIM)", () => {
  it("full coverage + one phase + no coexistence → '<phase> throughout' (unchanged)", () => {
    render(
      <PhaseStrip
        segments={segs({ phase: "Pn3m" }, { phase: "Pn3m" }, { phase: "Pn3m" })}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent(/Pn3m\s*throughout/);
    expect(cap).not.toHaveTextContent(/\bof\b/);
  });

  it("partial coverage + one phase → carries the fraction, never 'throughout'", () => {
    render(
      <PhaseStrip
        segments={segs({ phase: "Pn3m" }, { phase: null }, { phase: null })}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent(/Pn3m\s*in 1 of 3/);
    expect(cap).not.toHaveTextContent(/throughout/i);
  });

  it("form_factor and null-state cells count as uncovered in the fraction", () => {
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Im3m" },
          { phase: null, state: "form_factor" },
          { phase: null, state: "null" },
        )}
      />,
    );
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/in 1 of 3/);
  });

  it("SC-PREVCLAIM repro: one coexistence cell of three never reads 'throughout'", () => {
    // The lying caption was "Im3m throughout" for [Im3m+Lamellar, null, null].
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Im3m", coexistWith: ["Lamellar"] },
          { phase: null },
          { phase: null },
        )}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Im3m");
    expect(cap).toHaveTextContent("Lamellar");
    expect(cap).toHaveTextContent(/in 1 of 3/);
    expect(cap).not.toHaveTextContent(/throughout/i);
  });

  it("uniform coexistence in every cell → 'A + B throughout' (the pair IS throughout)", () => {
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Im3m", coexistWith: ["Lamellar"] },
          { phase: "Im3m", coexistWith: ["Lamellar"] },
        )}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Im3m");
    expect(cap).toHaveTextContent("Lamellar");
    expect(cap).toHaveTextContent(/throughout/);
  });

  it("non-uniform coexistence (some cells pure) → fraction, not 'throughout'", () => {
    // Lamellar is NOT throughout here, so the caption may not claim it is.
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Im3m", coexistWith: ["Lamellar"] },
          { phase: "Im3m" },
          { phase: "Im3m" },
        )}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Lamellar");
    expect(cap).toHaveTextContent(/in 3 of 3/);
    expect(cap).not.toHaveTextContent(/throughout/i);
  });

  it("full-coverage transition keeps its shape, with no fraction noise", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m" }, { phase: "Lamellar" })} />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("→");
    expect(cap).not.toHaveTextContent(/\bof\b/);
  });

  it("transition with unindexed gaps appends the fraction softly", () => {
    render(
      <PhaseStrip
        segments={segs({ phase: "Pn3m" }, { phase: null }, { phase: "Lamellar" })}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("→");
    expect(cap).toHaveTextContent(/in 2 of 3/);
  });

  it("the coexistence '+' glyph is decorative with sr 'with' (parity with '→'/'to')", () => {
    render(
      <PhaseStrip
        segments={segs(
          { phase: "Im3m", coexistWith: ["Lamellar"] },
          { phase: null },
        )}
      />,
    );
    const cap = screen.getByTestId("ps-cap");
    const plus = cap.querySelector('[aria-hidden="true"]');
    expect(plus).not.toBeNull();
    expect(plus).toHaveTextContent("+");
    const withWord = within(cap).getByText("with");
    expect(withWord).not.toHaveAttribute("aria-hidden");
  });
});

describe("<PhaseStrip> multi-segment count and layout", () => {
  it("renders one ps-seg per segment", () => {
    render(<PhaseStrip segments={segs({ phase: "Pn3m" }, { phase: "Lamellar" }, { phase: "Im3m" })} />);
    expect(screen.getAllByTestId("ps-seg")).toHaveLength(3);
  });

  it("treats a non-monotone strip as a transition (distinct-count, not first===last)", () => {
    // [Pn3m, Lamellar, Pn3m] has two distinct phases → transition, NOT "throughout".
    render(<PhaseStrip segments={segs({ phase: "Pn3m" }, { phase: "Lamellar" }, { phase: "Pn3m" })} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent(/→/);
    expect(cap).not.toHaveTextContent(/throughout/i);
  });

  it("labels each segment with its phase via aria-label and title (mixed indexed/null)", () => {
    render(<PhaseStrip segments={segs({ phase: "Pn3m" }, { phase: null })} />);
    const segsEls = screen.getAllByTestId("ps-seg");
    expect(segsEls[0]).toHaveAttribute("aria-label", "Pn3m");
    expect(segsEls[0]).toHaveAttribute("title", "Pn3m");
    expect(segsEls[1]).toHaveAttribute("aria-label", "Unindexed");
    expect(segsEls[1]).toHaveAttribute("title", "Unindexed");
  });
});

describe("<PhaseStrip> empty label", () => {
  it("renders the default empty label when no segment is indexed", () => {
    render(<PhaseStrip segments={segs({ phase: null }, { phase: null })} />);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/no clear phase/i);
  });

  it("honors an emptyLabel override", () => {
    render(
      <PhaseStrip
        segments={segs({ phase: null })}
        emptyLabel="Members not yet indexed; phase preview unavailable."
      />,
    );
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/not yet indexed/i);
  });
});

describe("<PhaseStrip> size and placement props", () => {
  it("exposes a data-size attribute reflecting the size prop (default md)", () => {
    const { rerender, container } = render(<PhaseStrip segments={segs({ phase: "Pn3m" })} />);
    expect(container.firstChild).toHaveAttribute("data-size", "md");
    rerender(<PhaseStrip segments={segs({ phase: "Pn3m" })} size="sm" />);
    expect(container.firstChild).toHaveAttribute("data-size", "sm");
  });

  it("applies the placement className to the root", () => {
    const { container } = render(<PhaseStrip segments={segs({ phase: "Pn3m" })} className="mt-5" />);
    expect(container.firstChild).toHaveClass("mt-5");
  });
});

describe("<PhaseStrip> data-state attribute", () => {
  it("renders a hollow dashed cell with data-state=form_factor for a form_factor segment", () => {
    render(<PhaseStrip segments={segs({ phase: null, state: "form_factor" })} />);
    const segment = screen.getAllByTestId("ps-seg")[0]!;
    expect(segment).toHaveAttribute("data-state", "form_factor");
    expect(segment).toHaveAttribute("aria-label", expect.stringMatching(/form factor/i));
  });

  it("renders a null-state cell with data-state=null distinct from form_factor", () => {
    render(<PhaseStrip segments={segs({ phase: null, state: "null" })} />);
    const segment = screen.getAllByTestId("ps-seg")[0]!;
    expect(segment).toHaveAttribute("data-state", "null");
    expect(segment.getAttribute("data-state")).not.toBe("form_factor");
  });
});

describe("<PhaseStrip> coexistWith empty/null", () => {
  it("treats an empty coexistWith array as single-phase", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m", coexistWith: [] })} />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.style.background).not.toContain("linear-gradient");
    expect(cell.getAttribute("data-coexist-count")).toBeNull();
    expect(cell.getAttribute("aria-label")).toBe("Pn3m");
  });

  it("treats a null coexistWith as single-phase", () => {
    render(
      <PhaseStrip segments={segs({ phase: "Pn3m", coexistWith: null })} />,
    );
    const cell = screen.getByTestId("ps-seg");
    expect(cell.style.background).not.toContain("linear-gradient");
    expect(cell.getAttribute("data-coexist-count")).toBeNull();
  });
});
