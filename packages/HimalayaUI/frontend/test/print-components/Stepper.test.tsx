import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Stepper } from "../../src/print/components/Stepper";

describe("<Stepper>", () => {
  it("renders the label text", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.getByText("Lipid 1‑2 + LL37 1:0.5")).toBeInTheDocument();
  });

  it("renders the position subtitle when provided", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" position="sample 4 of 9" />);
    expect(screen.getByText("sample 4 of 9")).toBeInTheDocument();
  });

  it("does not render a subtitle when position is omitted", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.queryByText("sample 4 of 9")).not.toBeInTheDocument();
  });

  it("renders two buttons with accessible names matching /previous/i and /next/i", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.getByRole("button", { name: /previous/i })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /next/i })).toBeInTheDocument();
  });

  it("calls onPrev when prev button is clicked", () => {
    const onPrev = vi.fn();
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" onPrev={onPrev} />);
    fireEvent.click(screen.getByRole("button", { name: /previous/i }));
    expect(onPrev).toHaveBeenCalledTimes(1);
  });

  it("calls onNext when next button is clicked", () => {
    const onNext = vi.fn();
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" onNext={onNext} />);
    fireEvent.click(screen.getByRole("button", { name: /next/i }));
    expect(onNext).toHaveBeenCalledTimes(1);
  });

  it("disables the prev button when prevDisabled is true", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" prevDisabled />);
    expect(screen.getByRole("button", { name: /previous/i })).toBeDisabled();
  });

  it("disables the next button when nextDisabled is true", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" nextDisabled />);
    expect(screen.getByRole("button", { name: /next/i })).toBeDisabled();
  });

  it("does not disable either button by default", () => {
    render(<Stepper label="Lipid 1‑2 + LL37 1:0.5" />);
    expect(screen.getByRole("button", { name: /previous/i })).not.toBeDisabled();
    expect(screen.getByRole("button", { name: /next/i })).not.toBeDisabled();
  });
});
