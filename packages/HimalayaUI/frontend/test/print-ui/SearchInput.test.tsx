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
