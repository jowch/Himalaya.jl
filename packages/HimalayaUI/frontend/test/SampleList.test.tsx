import { render, screen, fireEvent } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, it, expect, vi } from "vitest";
import { SampleList } from "../src/components/SampleList";
import type { Sample } from "../src/api";

const SAMPLES: Sample[] = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
  { id: 3, experiment_id: 1, label: "D3", name: "UL1", notes: null,
    tags: [{ id: 2, key: "peptide", value: "melittin", source: "manual" }] },
];

describe("<SampleList>", () => {
  it("renders one row per sample", () => {
    render(<SampleList samples={SAMPLES} activeId={undefined} onSelect={() => {}} />);
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(3);
  });

  it("marks the active sample", () => {
    render(<SampleList samples={SAMPLES} activeId={2} onSelect={() => {}} />);
    const active = document.querySelector('[data-active="true"]');
    expect(active?.getAttribute("data-sample-id")).toBe("2");
  });

  it("calls onSelect with sample id on click", () => {
    const onSelect = vi.fn();
    render(<SampleList samples={SAMPLES} activeId={undefined} onSelect={onSelect} />);
    const row = document.querySelector<HTMLElement>('[data-sample-id="3"]')!;
    fireEvent.click(row);
    expect(onSelect).toHaveBeenCalledWith(3);
  });

  it("filter narrows by label, name, and tag", async () => {
    const user = userEvent.setup();
    render(<SampleList samples={SAMPLES} activeId={undefined} onSelect={() => {}} />);
    const filter = screen.getByPlaceholderText(/filter samples/i);

    await user.type(filter, "D1");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(1);

    await user.clear(filter);
    await user.type(filter, "UX");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(2);

    await user.clear(filter);
    await user.type(filter, "lipid");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(1);

    await user.clear(filter);
    await user.type(filter, "melittin");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(1);

    await user.clear(filter);
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(3);
  });
});
