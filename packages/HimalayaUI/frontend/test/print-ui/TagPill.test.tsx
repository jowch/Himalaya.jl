import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagPill } from "../../src/print/ui/TagPill";

describe("TagPill", () => {
  it('has data-testid="tag-pill"', () => {
    render(<TagPill tag={{ key: "LL37" }} />);
    expect(screen.getByTestId("tag-pill")).toBeTruthy();
  });

  it('defaults to dense size: data-size="sm"', () => {
    render(<TagPill tag={{ key: "LL37" }} />);
    expect(screen.getByTestId("tag-pill").getAttribute("data-size")).toBe("sm");
  });

  it('size="md" overrides to data-size="md"', () => {
    render(<TagPill tag={{ key: "LL37" }} size="md" />);
    expect(screen.getByTestId("tag-pill").getAttribute("data-size")).toBe("md");
  });

  it("renders a key-only tag as the bare key, with no value span", () => {
    render(<TagPill tag={{ key: "LL37" }} />);
    const keyEl = screen.getByTestId("tag-pill").querySelector(
      '[data-role="tag-key"]',
    );
    expect(keyEl?.textContent).toBe("LL37");
    expect(
      screen.getByTestId("tag-pill").querySelector('[data-role="tag-value"]'),
    ).toBeNull();
  });

  it("renders a key+value tag as key + value", () => {
    render(<TagPill tag={{ key: "temperature", value: "37C" }} />);
    const pill = screen.getByTestId("tag-pill");
    expect(pill.querySelector('[data-role="tag-key"]')?.textContent).toBe(
      "temperature",
    );
    expect(pill.querySelector('[data-role="tag-value"]')?.textContent).toBe(
      "37C",
    );
  });

  it("renders NO remove button when onRemove is omitted", () => {
    render(<TagPill tag={{ key: "LL37" }} />);
    expect(screen.queryByRole("button", { name: "Remove" })).toBeNull();
  });

  it("renders a remove button when onRemove is provided and fires it on click", () => {
    const onRemove = vi.fn();
    render(<TagPill tag={{ key: "LL37" }} onRemove={onRemove} />);
    const button = screen.getByRole("button", { name: "Remove" });
    expect(button).toBeTruthy();
    fireEvent.click(button);
    expect(onRemove).toHaveBeenCalledTimes(1);
  });
});
