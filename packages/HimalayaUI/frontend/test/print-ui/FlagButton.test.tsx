import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { FlagButton } from "../../src/print/ui/FlagButton";

describe("<FlagButton>", () => {
  it("renders the value and is a button", () => {
    render(<FlagButton value="1 : 0.25" />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.tagName).toBe("BUTTON");
    expect(btn).toHaveTextContent("1 : 0.25");
  });
  it("defaults to the ok state with a re-open title", () => {
    render(<FlagButton value="1 : 1" />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.getAttribute("data-state")).toBe("ok");
    expect(btn.getAttribute("title")).toBe("Click to re-open this value");
    expect(screen.queryByText(/check the read/i)).toBeNull();
  });
  it("renders the flagged state with the check-the-read caption", () => {
    render(<FlagButton value="1 : 0" flagged />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.getAttribute("data-state")).toBe("flagged");
    expect(screen.getByText(/check the read/i)).toBeInTheDocument();
  });
  it("fires onClick", () => {
    const onClick = vi.fn();
    render(<FlagButton value="1 : 0" flagged onClick={onClick} />);
    fireEvent.click(screen.getByTestId("flag-button"));
    expect(onClick).toHaveBeenCalledOnce();
  });
});
