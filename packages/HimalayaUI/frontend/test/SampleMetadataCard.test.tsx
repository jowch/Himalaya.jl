import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { vi } from "vitest";
import { SampleMetadataCard } from "../src/components/SampleMetadataCard";
import type { Sample } from "../src/api";

const makeSample = (overrides: Partial<Sample> = {}): Sample => ({
  id: 1,
  experiment_id: 1,
  display_name: "DOPC 37C",
  name: "D1",
  notes: "run 2",
  tags: [
    { id: 10, key: "concentration", value: "10mM", source: "manifest" },
    { id: 11, key: "temp", value: "37C", source: "manual" },
  ],
  ...overrides,
});

test("renders display_name in input, stable name in breadcrumb, and tags", () => {
  render(
    <SampleMetadataCard
      sample={makeSample()}
      exposureSummary={{ total: 5, accepted: 3, rejected: 2 }}
      onUpdateSample={vi.fn()}
      onAddTag={vi.fn()}
      onRemoveTag={vi.fn()}
    />,
  );
  expect(screen.getByDisplayValue("DOPC 37C")).toBeInTheDocument();
  expect(screen.getByText("D1")).toBeInTheDocument(); // stable name in breadcrumb
  expect(screen.getByText("concentration: 10mM")).toBeInTheDocument();
  expect(
    screen.getByText("5 exposures · 3 accepted · 2 rejected"),
  ).toBeInTheDocument();
});

test("manifest tags have no delete button", () => {
  render(
    <SampleMetadataCard
      sample={makeSample()}
      exposureSummary={{ total: 5, accepted: 3, rejected: 2 }}
      onUpdateSample={vi.fn()}
      onAddTag={vi.fn()}
      onRemoveTag={vi.fn()}
    />,
  );
  const manifestTag = screen.getByText("concentration: 10mM").closest("span")!;
  expect(manifestTag.querySelector("button")).toBeNull();
});

test("calls onUpdateSample on name blur", async () => {
  const onUpdate = vi.fn();
  render(
    <SampleMetadataCard
      sample={makeSample()}
      exposureSummary={{ total: 1, accepted: 1, rejected: 0 }}
      onUpdateSample={onUpdate}
      onAddTag={vi.fn()}
      onRemoveTag={vi.fn()}
    />,
  );
  const nameInput = screen.getByDisplayValue("DOPC 37C");  // display_name value
  await userEvent.clear(nameInput);
  await userEvent.type(nameInput, "New Name");
  await userEvent.tab();
  expect(onUpdate).toHaveBeenCalledWith({ display_name: "New Name" });
});
