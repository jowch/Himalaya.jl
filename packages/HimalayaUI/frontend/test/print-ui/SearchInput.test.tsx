import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SearchInput } from "../../src/print/ui/SearchInput";

describe("<SearchInput>", () => {
  it("renders an input with the given placeholder", () => {
    render(<SearchInput value="" onChange={() => {}} placeholder="Search series" />);
    expect(screen.getByPlaceholderText("Search series")).toBeInTheDocument();
  });

  it("reflects value", () => {
    render(<SearchInput value="lipid" onChange={() => {}} placeholder="Search series" />);
    const input = screen.getByPlaceholderText("Search series") as HTMLInputElement;
    expect(input.value).toBe("lipid");
  });

  it("fires onChange with the new value on typing", () => {
    const onChange = vi.fn();
    render(<SearchInput value="" onChange={onChange} placeholder="Search series" />);
    const input = screen.getByPlaceholderText("Search series");
    fireEvent.change(input, { target: { value: "abc" } });
    expect(onChange).toHaveBeenCalledWith("abc");
  });

  it("marks the magnifier icon aria-hidden", () => {
    const { container } = render(
      <SearchInput value="" onChange={() => {}} placeholder="Search series" />,
    );
    const svg = container.querySelector("svg");
    expect(svg).not.toBeNull();
    expect(svg?.getAttribute("aria-hidden")).toBe("true");
  });
});

describe("<SearchInput> inline clear (FOL P2-4)", () => {
  it("offers no clear affordance while empty (a control that cannot act is removed)", () => {
    render(<SearchInput value="" onChange={() => {}} placeholder="Search series" />);
    expect(screen.queryByRole("button", { name: "Clear search" })).toBeNull();
  });

  it("shows a one-click clear once non-empty", () => {
    render(<SearchInput value="lipid" onChange={() => {}} placeholder="Search series" />);
    expect(screen.getByRole("button", { name: "Clear search" })).toBeInTheDocument();
  });

  it("clicking the clear fires onChange('')", () => {
    const onChange = vi.fn();
    render(<SearchInput value="lipid" onChange={onChange} placeholder="Search series" />);
    fireEvent.click(screen.getByRole("button", { name: "Clear search" }));
    expect(onChange).toHaveBeenCalledWith("");
  });

  it("focus stays in the input after clearing (the user is mid-search)", () => {
    const onChange = vi.fn();
    render(<SearchInput value="lipid" onChange={onChange} placeholder="Search series" />);
    fireEvent.click(screen.getByRole("button", { name: "Clear search" }));
    expect(document.activeElement).toBe(screen.getByPlaceholderText("Search series"));
  });

  // The >=24x24 CSS-px hit box (WCAG 2.5.8 / the LO-TAGTARGET precedent) is a
  // geometry fact jsdom cannot measure; the affordance is pinned semantically
  // here and the box uses the same in-flow h-6/w-6 idiom the tag-remove x
  // already pins in a real browser (e2e/corpus-culling.spec.ts).
});
