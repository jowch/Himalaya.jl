import { describe, it, expect, vi } from "vitest";
import { act, render, screen, fireEvent } from "@testing-library/react";
import { SampleList } from "../../src/ui/SampleList";
import type { CorpusSample } from "../../src/api";

const sample = (id: number, name: string) => ({ id, name } as unknown as CorpusSample);

describe("SampleList filter keys", () => {
  it("leaves Enter to an input method that is still composing", () => {
    const onPick = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001")]} currentId={undefined} filter="JC" onFilter={() => {}} onPick={onPick} />,
    );
    const box = screen.getByLabelText("Filter samples");
    fireEvent.keyDown(box, { key: "Enter", isComposing: true });
    expect(onPick).not.toHaveBeenCalled();
    fireEvent.keyDown(box, { key: "Enter" });
    expect(onPick).toHaveBeenCalledWith(1);
  });

  it("Tab and Shift+Tab move a highlight through the matches with wrap, focus stays in the field, and Enter opens the highlighted one", () => {
    Element.prototype.scrollIntoView = vi.fn();
    const onPick = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001"), sample(2, "JC002"), sample(3, "JC003")]} currentId={undefined} filter="JC" onFilter={() => {}} onPick={onPick} />,
    );
    const box = screen.getByLabelText("Filter samples");
    const highlighted = (): string | null | undefined =>
      document.querySelector("[data-testid=workspace-samples] [data-highlighted]")?.getAttribute("data-sample-id");
    expect(highlighted()).toBeUndefined();
    act(() => box.focus());
    expect(highlighted()).toBe("1");
    expect(box).toHaveAttribute("aria-activedescendant", "workspace-sample-1");
    fireEvent.keyDown(box, { key: "Tab" });
    expect(highlighted()).toBe("2");
    expect(document.activeElement).toBe(box);
    fireEvent.keyDown(box, { key: "Tab", shiftKey: true });
    fireEvent.keyDown(box, { key: "Tab", shiftKey: true });
    expect(highlighted()).toBe("3");
    fireEvent.keyDown(box, { key: "Tab" });
    expect(highlighted()).toBe("1");
    fireEvent.keyDown(box, { key: "ArrowDown" });
    expect(highlighted()).toBe("2");
    fireEvent.keyDown(box, { key: "Enter" });
    expect(onPick).toHaveBeenCalledWith(2);
  });

  it("starts the highlight over on the first match when the filter changes, and Esc leaves the field without picking", () => {
    Element.prototype.scrollIntoView = vi.fn();
    const onPick = vi.fn();
    const list = (filter: string) => (
      <SampleList samples={filter === "JC" ? [sample(1, "JC001"), sample(2, "JC002")] : [sample(2, "JC002")]} currentId={undefined} filter={filter} onFilter={() => {}} onPick={onPick} />
    );
    const { rerender } = render(list("JC"));
    const box = screen.getByLabelText("Filter samples");
    act(() => box.focus());
    fireEvent.keyDown(box, { key: "Tab" });
    rerender(list("JC002"));
    expect(document.querySelector("[data-highlighted]")?.getAttribute("data-sample-id")).toBe("2");
    fireEvent.keyDown(box, { key: "Escape" });
    expect(document.activeElement).not.toBe(box);
    expect(document.querySelector("[data-highlighted]")).toBeNull();
    expect(onPick).not.toHaveBeenCalled();
  });
});

describe("SampleList weights", () => {
  it("sets only the current sample's name in semibold, so the cursor row stands out", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001"), sample(2, "JC002")]} currentId={2} filter="" onFilter={() => {}} onPick={() => {}} />,
    );
    expect(screen.getByRole("button", { name: /JC002/ }).className).toMatch(/\bfont-semibold\b/);
    expect(screen.getByRole("button", { name: /JC001/ }).className).not.toMatch(/\bfont-semibold\b/);
  });
});

describe("SampleList selection", () => {
  it("tints the current row flat, like an Assignment choice, with no accent bar", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001"), sample(2, "JC002")]} currentId={2} filter="" onFilter={() => {}} onPick={() => {}} />,
    );
    const cur = screen.getByRole("button", { name: /JC002/ }).className;
    const other = screen.getByRole("button", { name: /JC001/ }).className;
    expect(cur).toMatch(/\bbg-accent\/25\b/);
    expect(cur).not.toMatch(/shadow-|border-l|ring-/);
    expect(other).not.toMatch(/\bbg-accent\/25\b/);
    expect(cur.replace(/\b(font-semibold|bg-accent\/25|hover:\S+)\b/g, "").split(/\s+/).sort()).toEqual(
      other.replace(/\b(font-semibold|bg-accent\/25|hover:\S+)\b/g, "").split(/\s+/).sort(),
    );
  });
});

describe("SampleList comparisons", () => {
  it("marks compared rows", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001"), sample(2, "JC002")]} currentId={1} filter="" onFilter={() => {}} onPick={() => {}} compareColors={new Map([[1, "blue"], [2, "red"]])} hiddenCompares={new Set([2])} />,
    );
    expect(screen.getByRole("img", { name: "Compared, hidden" }).getAttribute("class")).toMatch(/opacity-45/);
    expect(screen.getByRole("img", { name: "Compared" }).getAttribute("class")).not.toMatch(/opacity-45/);
    expect(screen.getByRole("button", { name: /JC002/ })).toHaveAttribute("data-compared");
  });

  it("a modifier click toggles the comparison instead of opening the sample", () => {
    const onPick = vi.fn();
    const onToggleCompare = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001")]} currentId={undefined} filter="" onFilter={() => {}} onPick={onPick} onToggleCompare={onToggleCompare} />,
    );
    const row = screen.getByRole("button", { name: /JC001/ });
    fireEvent.click(row, { metaKey: true });
    fireEvent.click(row, { ctrlKey: true });
    expect(onToggleCompare).toHaveBeenCalledTimes(2);
    expect(onPick).not.toHaveBeenCalled();
    fireEvent.click(row);
    expect(onPick).toHaveBeenCalledWith(1);
  });
});

describe("SampleList setup headers", () => {
  const inSetup = (id: number, name: string, setup_id: number) => ({ id, name, setup_id } as unknown as CorpusSample);

  it("heads each run of a setup's samples with its name", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList
        samples={[inSetup(1, "A", 3), inSetup(2, "B", 3), inSetup(3, "C", 4)]}
        setupHeaders={new Map([[3, "Long pipe"], [4, "setup2"]])}
        currentId={1} filter="" onFilter={() => {}} onPick={() => {}}
      />,
    );
    expect(screen.getAllByTestId("setup-header").map((h) => h.textContent)).toEqual(["Long pipe", "setup2"]);
    const items = screen.getByTestId("workspace-samples").querySelectorAll(":scope > li");
    expect([...items].map((li) => li.textContent?.trim().slice(0, 9))).toEqual(["Long pipe", "A○", "B○", "setup2", "C○"]);
  });

  it("shows no headers without them, as in a one-setup workspace", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(<SampleList samples={[inSetup(1, "A", 3)]} currentId={1} filter="" onFilter={() => {}} onPick={() => {}} />);
    expect(screen.queryByTestId("setup-header")).toBeNull();
  });
});

describe("SampleList tooltips", () => {
  it("names the unreviewed mark and the bridge count on hover", async () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001")]} currentId={1} filter="" onFilter={() => {}} onPick={() => {}} bridgeCounts={new Map([[1, 2]])} />,
    );
    fireEvent.mouseEnter(screen.getByRole("img", { name: "Unreviewed" }));
    expect(await screen.findByRole("tooltip")).toHaveTextContent("Unreviewed");
    fireEvent.mouseLeave(screen.getByRole("img", { name: "Unreviewed" }));
    fireEvent.mouseEnter(screen.getByText("⇄2"));
    expect(await screen.findByRole("tooltip")).toHaveTextContent("2 bridges");
    expect(screen.getByText("⇄2").parentElement).not.toHaveAttribute("aria-describedby");
  });
});

describe("SampleList analysis gaps", () => {
  it("says on the row why nothing is analyzed", () => {
    Element.prototype.scrollIntoView = vi.fn();
    const rows = [
      { ...sample(1, "JC001"), analysis_gap: "all_dropped" },
      { ...sample(2, "JC002"), analysis_gap: "none_selected" },
      sample(3, "JC003"),
    ] as CorpusSample[];
    render(<SampleList samples={rows} currentId={3} filter="" onFilter={() => {}} onPick={() => {}} />);
    expect(screen.getByRole("img", { name: "Unreviewed. Every image is dropped, so nothing is analyzed." })).toHaveTextContent("All dropped");
    expect(screen.getByRole("img", { name: "Unreviewed. No image is selected, so nothing is analyzed." })).toHaveTextContent("None selected");
    expect(screen.getByRole("img", { name: "Unreviewed" })).toHaveTextContent("○");
  });
});
