import { act, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { CompareBar, type CompareChip } from "../../src/ui/CompareBar";
import { resetTooltipWarmth, TOOLTIP_OPEN_DELAY_MS } from "../../src/ui/Tooltip";

const chip = (over: Partial<CompareChip> = {}): CompareChip => ({ id: 2, name: "MO_PEG5", color: "#888", scale: undefined, ...over });

function bar(chips: CompareChip[], subtractId?: number) {
  const noop = vi.fn();
  render(<CompareBar chips={chips} subtractId={subtractId} hidden={false} modKey="⌘" onToggleShown={noop} onSubtract={noop} onRemove={noop} onAdd={noop} />);
}

function hover(name: string): HTMLElement {
  fireEvent.mouseEnter(screen.getByRole("button", { name }));
  act(() => {
    vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
  });
  return screen.getByRole("tooltip");
}

describe("CompareBar chip tooltip", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("names the chip and its setup and says how to hide and subtract", () => {
    bar([chip({ setup: "setup2" })]);
    expect(hover("MO_PEG5, setup2")).toHaveTextContent("MO_PEG5 · setup2Click to hide. ⌘-click to subtract it from the current trace.");
  });

  it("says only how to hide when the open sample cannot be subtracted from", () => {
    const noop = vi.fn();
    render(<CompareBar chips={[chip()]} subtractId={undefined} hidden={false} modKey="⌘" canSubtract={false} onToggleShown={noop} onSubtract={noop} onRemove={noop} onAdd={noop} />);
    expect(hover("MO_PEG5")).toHaveTextContent("MO_PEG5Click to hide.");
    expect(screen.getByRole("tooltip")).not.toHaveTextContent("subtract");
  });

  it("says how to stop subtracting on the subtracted chip, with no setup when it is the current one", () => {
    bar([chip()], 2);
    expect(hover("MO_PEG5")).toHaveTextContent("MO_PEG5Click to hide. ⌘-click to stop subtracting.");
  });

  it("repeats the scale's note on the chip, so keyboard focus reaches it", () => {
    bar([chip({ scale: { scale: 1, source: "monitor" } as CompareChip["scale"] })]);
    expect(hover("MO_PEG5")).toHaveTextContent("Scaled by the transmitted beam");
  });

  it("carries no native title on the chip button", () => {
    bar([chip({ setup: "setup2" })]);
    expect(screen.getByRole("button", { name: "MO_PEG5, setup2" })).not.toHaveAttribute("title");
  });
});

describe("CompareBar scale tooltip", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("explains a high-q match on the scale, with no native title", () => {
    bar([chip({ scale: { scale: 1, source: "matched" } as CompareChip["scale"] })]);
    const scale = document.querySelector("[data-scale-source]") as HTMLElement;
    expect(scale).not.toHaveAttribute("title");
    fireEvent.mouseEnter(scale);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent("No transmitted-beam reading, so matched at high q");
  });

  it("tells a capillary setup apart from one whose type is not set yet", () => {
    for (const [source, note] of [["capillary", "Capillary setup, so matched at high q"], ["unknown_mode", "Setup type not set, so matched at high q"]] as const) {
      resetTooltipWarmth();
      const { unmount } = render(<CompareBar chips={[chip({ scale: { scale: 1, source } })]} subtractId={undefined} hidden={false} modKey="⌘" onToggleShown={() => {}} onSubtract={() => {}} onRemove={() => {}} onAdd={() => {}} />);
      fireEvent.mouseEnter(document.querySelector("[data-scale-source]") as HTMLElement);
      act(() => {
        vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
      });
      expect(screen.getByRole("tooltip")).toHaveTextContent(note);
      unmount();
    }
  });

  it("tells a screen reader why a chip was matched at high q, and says nothing extra for a monitor-scaled one", () => {
    bar([
      chip({ scale: { scale: 1, source: "other_setup" } as CompareChip["scale"] }),
      chip({ id: 3, name: "MO_PEG6", scale: { scale: 1, source: "monitor" } as CompareChip["scale"] }),
    ]);
    expect(screen.getByRole("button", { name: "MO_PEG5" })).toHaveAccessibleDescription("From another setup, so matched at high q");
    expect(screen.getByRole("button", { name: "MO_PEG6" })).toHaveAccessibleDescription("");
  });

  it("shows only the factor on a chip from another setup, and says in its tooltip that it was matched at high q", () => {
    bar([chip({ scale: { scale: 1, source: "other_setup" } as CompareChip["scale"] })]);
    const scale = document.querySelector("[data-scale-source]") as HTMLElement;
    expect(scale.textContent).toBe("×1.00");
    fireEvent.mouseEnter(scale);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent("From another setup, so matched at high q");
  });
});

describe("CompareBar in Minus view", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("names what is subtracted from a drawn chip, and keeps the pair's scale", () => {
    bar([chip({ minus: { kind: "drawn", background: "lipidA 01A (B0100)" }, scale: { scale: 1, source: "normalised" } })]);
    expect(hover("MO_PEG5")).toHaveTextContent("Minus lipidA 01A (B0100)");
    expect(document.querySelector("[data-scale-source='normalised']")).not.toBeNull();
    expect(document.querySelector("[data-minus-note]")).toBeNull();
  });

  it("says in the tooltip that a pair was matched at high q, and shows only its factor", () => {
    bar([chip({ minus: { kind: "drawn", background: "buf" }, scale: { scale: 0.9, source: "beam_off" } })]);
    expect(hover("MO_PEG5")).toHaveTextContent("Taken with the beam off, so matched at high q");
    expect(document.querySelector("[data-scale-source='beam_off']")!.textContent).toBe("×0.900");
  });

  it.each([
    ["buffer", "Not subtracted: this is a buffer"],
    ["none", "Not subtracted: its background is set to none"],
    ["hidden", "Not subtracted: the setup of its background is hidden"],
  ] as const)("says on the chip and in its tooltip why a %s chip is not drawn", (reason, note) => {
    bar([chip({ minus: { kind: "skipped", reason } })]);
    expect(document.querySelector(`[data-minus-note='${reason}']`)).toHaveTextContent("not subtracted");
    expect(hover("MO_PEG5")).toHaveTextContent(note);
    expect(document.querySelector("[data-scale-source]")).toBeNull();
  });

  it("reads a background still loading as loading, not as an error", () => {
    bar([chip({ minus: { kind: "skipped", reason: "loading" } })]);
    expect(document.querySelector("[data-minus-note='loading']")).toHaveTextContent("loading");
    expect(document.querySelector("[data-minus-note='loading']")).not.toHaveTextContent("not subtracted");
    expect(hover("MO_PEG5")).toHaveTextContent("Loading the background trace");
  });

  it("gives a chip that is not drawn its full reason as an accessible description, without hover", () => {
    bar([chip({ minus: { kind: "skipped", reason: "buffer" } })]);
    expect(screen.getByRole("button", { name: "MO_PEG5" })).toHaveAccessibleDescription("Not subtracted: this is a buffer");
  });

  it("gives a loading chip the loading note, and a drawn chip what is subtracted", () => {
    bar([chip({ minus: { kind: "skipped", reason: "loading" } }), chip({ id: 3, name: "MO_PEG6", minus: { kind: "drawn", background: "buf" } })]);
    expect(screen.getByRole("button", { name: "MO_PEG5" })).toHaveAccessibleDescription("Loading the background trace");
    expect(screen.getByRole("button", { name: "MO_PEG6" })).toHaveAccessibleDescription("Minus buf");
  });

  it("keeps the short note whole and apart from the name, which is the only part that truncates", () => {
    bar([chip({ name: "a very long sample name that cannot fit in the chip", minus: { kind: "skipped", reason: "buffer" } })]);
    const note = document.querySelector("[data-minus-note='buffer']") as HTMLElement;
    const name = document.querySelector("[data-chip-name]") as HTMLElement;
    expect(note.textContent).toBe("not subtracted");
    expect(name.contains(note)).toBe(false);
    expect(note.hasAttribute("data-chip-name")).toBe(false);
    expect(document.querySelectorAll("[data-chip-truncates]")).toHaveLength(1);
    expect(name).toHaveAttribute("data-chip-truncates");
  });

  it("still says how to hide and subtract on a chip that is not drawn", () => {
    bar([chip({ minus: { kind: "skipped", reason: "buffer" } })]);
    expect(hover("MO_PEG5")).toHaveTextContent("Click to hide. ⌘-click to subtract it from the current trace.");
  });
});
