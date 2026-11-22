import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
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
  });
});
