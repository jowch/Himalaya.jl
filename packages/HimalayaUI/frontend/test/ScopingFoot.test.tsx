import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingFoot } from "../src/components/ScopingFoot";

describe("ScopingFoot", () => {
  it("warns and disables the build when values are flagged", () => {
    render(<ScopingFoot flagCount={1} memberCount={6} keyLabel="ratio" canBuild={false} onBuild={() => {}} />);
    expect(screen.getByTestId("scoping-foot-state")).toHaveTextContent(/to check before you can build/i);
    expect(screen.getByTestId("scoping-open-confirm")).toBeDisabled();
  });
  it("reads ready and enables the build when clear", () => {
    render(<ScopingFoot flagCount={0} memberCount={6} keyLabel="ratio" canBuild onBuild={() => {}} />);
    expect(screen.getByTestId("scoping-foot-state")).toHaveTextContent(/ready to build/i);
    expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled();
  });
  it("fires onBuild when the enabled button is clicked", () => {
    const onBuild = vi.fn();
    render(<ScopingFoot flagCount={0} memberCount={1} keyLabel="ratio" canBuild onBuild={onBuild} />);
    fireEvent.click(screen.getByTestId("scoping-open-confirm"));
    expect(onBuild).toHaveBeenCalledTimes(1);
  });
  it("the ready dot reads through the success (sage) token, not an inline oklch literal", () => {
    render(<ScopingFoot flagCount={0} memberCount={1} keyLabel="ratio" canBuild onBuild={() => {}} />);
    const dot = screen.getByTestId("scoping-foot-state").querySelector("span");
    expect(dot).not.toBeNull();
    expect(dot!.getAttribute("style")).toContain("var(--color-success)");
  });
  it("the build button carries Print ink tokens, not ice-blue accent (S-B/S-C)", () => {
    render(<ScopingFoot flagCount={0} memberCount={1} keyLabel="ratio" canBuild onBuild={() => {}} />);
    const btn = screen.getByTestId("scoping-open-confirm");
    expect(btn.className).toContain("bg-ink");
    expect(btn.className).toContain("text-paper");
    expect(btn.className).not.toContain("bg-accent");
  });
});
