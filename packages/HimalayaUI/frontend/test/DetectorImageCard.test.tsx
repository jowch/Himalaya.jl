import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { vi } from "vitest";
import { DetectorImageCard } from "../src/components/DetectorImageCard";
import type { Exposure } from "../src/api";

vi.mock("../src/components/DetectorImage", () => ({
  DetectorImage: () => <canvas data-testid="mock-detector-image" />,
}));

const makeExposure = (overrides: Partial<Exposure>): Exposure => ({
  id: 1,
  sample_id: 1,
  filename: "pos1.dat",
  kind: "file",
  selected: false,
  status: null,
  image_path: "/tmp/foo.tiff",
  tags: [],
  sources: [],
  ...overrides,
});

test("shows Accept and Reject buttons", () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({})}
      onSetStatus={vi.fn()}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />,
  );
  expect(screen.getByRole("button", { name: /accept/i })).toBeInTheDocument();
  expect(screen.getByRole("button", { name: /reject/i })).toBeInTheDocument();
});

test("Reject shows rejection note input after click", async () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({})}
      onSetStatus={vi.fn()}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />,
  );
  await userEvent.click(screen.getByRole("button", { name: /reject/i }));
  expect(screen.getByPlaceholderText(/reason/i)).toBeInTheDocument();
});

test("Use for indexing is disabled when rejected", () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({ status: "rejected" })}
      onSetStatus={vi.fn()}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />,
  );
  expect(screen.getByRole("button", { name: /indexing/i })).toBeDisabled();
});

test("Accept calls onSetStatus with 'accepted'", async () => {
  const onSetStatus = vi.fn();
  render(
    <DetectorImageCard
      exposure={makeExposure({})}
      onSetStatus={onSetStatus}
      onSetIndexing={vi.fn()}
      onAddTag={vi.fn()}
    />,
  );
  await userEvent.click(screen.getByRole("button", { name: /accept/i }));
  expect(onSetStatus).toHaveBeenCalledWith("accepted");
});
