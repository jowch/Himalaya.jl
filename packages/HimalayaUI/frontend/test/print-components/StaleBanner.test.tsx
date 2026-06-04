import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { StaleBanner } from "../../src/print/components/StaleBanner";

describe("StaleBanner", () => {
  it("renders nothing when staleCount is 0", () => {
    const { container } = render(
      <StaleBanner staleCount={0} onReanalyze={vi.fn()} />,
    );
    expect(container.firstChild).toBeNull();
    expect(screen.queryByRole("alert")).toBeNull();
  });

  it("renders nothing when staleCount is negative", () => {
    const { container } = render(
      <StaleBanner staleCount={-1} onReanalyze={vi.fn()} />,
    );
    expect(container.firstChild).toBeNull();
  });

  it("renders alert with singular copy for staleCount=1", () => {
    render(<StaleBanner staleCount={1} onReanalyze={vi.fn()} />);
    const alert = screen.getByRole("alert");
    expect(alert).toBeTruthy();
    expect(alert.textContent).toMatch(/1 index is stale/);
  });

  it("renders alert with plural copy for staleCount=3", () => {
    render(<StaleBanner staleCount={3} onReanalyze={vi.fn()} />);
    const alert = screen.getByRole("alert");
    expect(alert.textContent).toMatch(/3 indices are stale/);
  });

  it("fires onReanalyze when button is clicked", async () => {
    const onReanalyze = vi.fn();
    render(<StaleBanner staleCount={2} onReanalyze={onReanalyze} />);
    await userEvent.click(screen.getByRole("button", { name: /re-analyze/i }));
    expect(onReanalyze).toHaveBeenCalledOnce();
  });

  it("disables the button and shows 'Re-analyzing…' when pending=true", () => {
    render(
      <StaleBanner staleCount={2} pending={true} onReanalyze={vi.fn()} />,
    );
    const btn = screen.getByRole("button");
    expect(btn).toBeDisabled();
    expect(btn.textContent).toMatch(/Re-analyzing/);
  });
});
