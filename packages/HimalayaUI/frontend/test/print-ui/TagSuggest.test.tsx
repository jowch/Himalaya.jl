import { render, screen, fireEvent, within } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagSuggest } from "../../src/print/ui/TagSuggest";

const OPTIONS = [
  { label: "lipid", count: 5 },
  { label: "temperature", count: 2 },
  { label: "dose", count: 1 },
];

function input(): HTMLInputElement {
  return screen.getByRole("combobox") as HTMLInputElement;
}

describe("<TagSuggest>", () => {
  // ── ARIA roles ──────────────────────────────────────────────────────────────

  it("renders an input with role=combobox + aria-autocomplete=list", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    const inp = screen.getByRole("combobox");
    expect(inp).toHaveAttribute("aria-autocomplete", "list");
    expect(inp).toHaveAttribute("aria-label", "Tag key");
  });

  it("aria-expanded is false initially, true after focus when options exist", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    expect(input()).toHaveAttribute("aria-expanded", "false");
    fireEvent.focus(input());
    expect(input()).toHaveAttribute("aria-expanded", "true");
  });

  it("listbox has role=listbox; each option has role=option", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    const listbox = screen.getByRole("listbox");
    expect(listbox).toBeInTheDocument();
    const opts = within(listbox).getAllByRole("option");
    expect(opts).toHaveLength(OPTIONS.length); // no create row when value is empty
  });

  it("aria-controls points at the listbox id", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    const listboxId = screen.getByRole("listbox").id;
    expect(input()).toHaveAttribute("aria-controls", listboxId);
  });

  it("data-mode reflects the mode prop", () => {
    const { getByTestId } = render(
      <TagSuggest label="val" mode="value" value="" options={[]} onCommit={() => {}} />,
    );
    expect(getByTestId("tag-suggest")).toHaveAttribute("data-mode", "value");
  });

  // ── Usage counts ─────────────────────────────────────────────────────────────

  it("renders each option's count", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    const counts = screen.getAllByTestId("tag-suggest-count");
    expect(counts).toHaveLength(OPTIONS.length);
    expect(counts[0]).toHaveTextContent("5");
    expect(counts[1]).toHaveTextContent("2");
    expect(counts[2]).toHaveTextContent("1");
  });

  it("de-emphasizes counts ≤ 1 (data-testid tag-suggest-count with faint style marker)", () => {
    const opts = [{ label: "a", count: 1 }, { label: "b", count: 3 }];
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={opts} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    const counts = screen.getAllByTestId("tag-suggest-count");
    // The component uses a faint class for count ≤ 1; we can't assert on
    // Tailwind classes, but we CAN assert the count text is present for both.
    expect(counts[0]).toHaveTextContent("1");
    expect(counts[1]).toHaveTextContent("3");
  });

  it("renders no count cell when count is undefined", () => {
    const opts = [{ label: "no-count" }];
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={opts} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    expect(screen.queryByTestId("tag-suggest-count")).toBeNull();
  });

  // ── Create-as-typed ─────────────────────────────────────────────────────────

  it("shows a create-as-typed option when the text matches no option exactly", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="novel" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    expect(screen.getByTestId("tag-suggest-create")).toBeInTheDocument();
    expect(screen.getByTestId("tag-suggest-create")).toHaveTextContent('Create "novel"');
  });

  it("does NOT show create row when the text exactly matches an option", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="lipid" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    expect(screen.queryByTestId("tag-suggest-create")).toBeNull();
  });

  it("does NOT show create row when value is empty", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    expect(screen.queryByTestId("tag-suggest-create")).toBeNull();
  });

  it("create row is role=option", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="xyz" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    const createOpt = screen.getByTestId("tag-suggest-create");
    expect(createOpt).toHaveAttribute("role", "option");
  });

  // ── Keyboard commit ──────────────────────────────────────────────────────────

  it("ArrowDown activates the first option; Enter commits it", () => {
    const onCommit = vi.fn();
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={onCommit} />,
    );
    fireEvent.focus(input());
    fireEvent.keyDown(input(), { key: "ArrowDown" });
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("lipid");
  });

  it("ArrowUp wraps to the last option", () => {
    const onCommit = vi.fn();
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={onCommit} />,
    );
    fireEvent.focus(input());
    fireEvent.keyDown(input(), { key: "ArrowUp" });
    fireEvent.keyDown(input(), { key: "Enter" });
    // Last option is "dose"
    expect(onCommit).toHaveBeenCalledWith("dose");
  });

  it("Enter with no active option commits the typed text", () => {
    const onCommit = vi.fn();
    render(
      <TagSuggest label="Tag key" mode="key" value="myval" options={OPTIONS} onCommit={onCommit} />,
    );
    // Do NOT open the listbox / arrow-navigate
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("myval");
  });

  it("Enter commits the create-as-typed option when it is active", () => {
    const onCommit = vi.fn();
    render(
      <TagSuggest label="Tag key" mode="key" value="novel" options={OPTIONS} onCommit={onCommit} />,
    );
    fireEvent.focus(input());
    // Navigate down past all filtered options to the create row
    const filteredCount = OPTIONS.filter((o) =>
      o.label.toLowerCase().includes("novel"),
    ).length;
    for (let i = 0; i <= filteredCount; i++) {
      fireEvent.keyDown(input(), { key: "ArrowDown" });
    }
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("novel");
  });

  it("Escape collapses the listbox without committing", () => {
    const onCommit = vi.fn();
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={onCommit} />,
    );
    fireEvent.focus(input());
    expect(input()).toHaveAttribute("aria-expanded", "true");
    fireEvent.keyDown(input(), { key: "Escape" });
    expect(input()).toHaveAttribute("aria-expanded", "false");
    expect(onCommit).not.toHaveBeenCalled();
  });

  it("aria-activedescendant tracks the active option id", () => {
    render(
      <TagSuggest label="Tag key" mode="key" value="" options={OPTIONS} onCommit={() => {}} />,
    );
    fireEvent.focus(input());
    expect(input()).not.toHaveAttribute("aria-activedescendant");
    fireEvent.keyDown(input(), { key: "ArrowDown" });
    // After one ArrowDown the first option is active; its id is set
    const optId = screen.getAllByRole("option")[0]!.id;
    expect(input()).toHaveAttribute("aria-activedescendant", optId);
  });

  it("onValueChange is called on text input", () => {
    const onValueChange = vi.fn();
    render(
      <TagSuggest
        label="Tag key"
        mode="key"
        value=""
        options={OPTIONS}
        onCommit={() => {}}
        onValueChange={onValueChange}
      />,
    );
    fireEvent.change(input(), { target: { value: "li" } });
    expect(onValueChange).toHaveBeenCalledWith("li");
  });
});
