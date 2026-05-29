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
  image_version: "",
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
  // Assert on the data-rejected attribute the component intentionally
  // exposes, not on the `opacity-40` Tailwind class — see CLAUDE.md and
  // the matching E2E fix in commit 3eddbd9.
  const rejected = screen.getByTestId("thumb-cell-2");
  expect(rejected).toHaveAttribute("data-rejected", "true");
  const good = screen.getByTestId("thumb-cell-1");
  expect(good).not.toHaveAttribute("data-rejected");
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

test("the thumb cell is a keyboard-operable button (role, tabIndex, aria-pressed)", () => {
  const exposures = [makeExposure({ id: 5 })];
  render(
    <ThumbnailGallery exposures={exposures} selectedId={5} onSelect={vi.fn()} />,
  );
  const cell = screen.getByTestId("thumb-cell-5");
  expect(cell).toHaveAttribute("role", "button");
  expect(cell).toHaveAttribute("tabindex", "0");
  // The viewed cell announces its selected state to assistive tech.
  expect(cell).toHaveAttribute("aria-pressed", "true");
});

test("Enter and Space select the focused thumbnail", async () => {
  const onSelect = vi.fn();
  const exposures = [makeExposure({ id: 5 })];
  render(
    <ThumbnailGallery exposures={exposures} selectedId={undefined} onSelect={onSelect} />,
  );
  const cell = screen.getByTestId("thumb-cell-5");
  cell.focus();
  await userEvent.keyboard("{Enter}");
  await userEvent.keyboard(" ");
  expect(onSelect).toHaveBeenCalledTimes(2);
  expect(onSelect).toHaveBeenCalledWith(5);
});
