import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { vi } from "vitest";
import { ThumbnailGallery } from "../src/components/ThumbnailGallery";
import type { Exposure } from "../src/api";

vi.mock("../src/components/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

const makeExposure = (overrides: Partial<Exposure>): Exposure => ({
  id: 1,
  sample_id: 1,
  filename: "pos1.dat",
  kind: "file",
  selected: false,
  status: null,
  image_path: null,
  tags: [],
  sources: [],
  ...overrides,
});

test("dims rejected exposures", () => {
  const exposures = [
    makeExposure({ id: 1, filename: "good.dat" }),
    makeExposure({ id: 2, filename: "bad.dat", status: "rejected" }),
  ];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={1}
      onSelect={vi.fn()}
    />,
  );
  const rejected = screen.getByTestId("thumb-cell-2");
  expect(rejected).toHaveClass("opacity-40");
  const good = screen.getByTestId("thumb-cell-1");
  expect(good).not.toHaveClass("opacity-40");
});

test("shows indexing chip on selected=true exposure", () => {
  const exposures = [makeExposure({ id: 1, selected: true })];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={1}
      onSelect={vi.fn()}
    />,
  );
  expect(screen.getByText("⊙ Indexing")).toBeInTheDocument();
});

test("no indexing chip when selected=false", () => {
  const exposures = [makeExposure({ id: 1, selected: false })];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={undefined}
      onSelect={vi.fn()}
    />,
  );
  expect(screen.queryByText("⊙ Indexing")).toBeNull();
});

test("calls onSelect when thumbnail clicked", async () => {
  const onSelect = vi.fn();
  const exposures = [makeExposure({ id: 5, filename: "pos5.dat" })];
  render(
    <ThumbnailGallery
      exposures={exposures}
      selectedId={undefined}
      onSelect={onSelect}
    />,
  );
  await userEvent.click(screen.getByTestId("thumb-cell-5"));
  expect(onSelect).toHaveBeenCalledWith(5);
});
