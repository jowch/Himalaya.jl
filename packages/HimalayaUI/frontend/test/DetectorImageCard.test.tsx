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
  image_version: "v1-1700000000",
  tags: [],
  sources: [],
  ...overrides,
});

const noopHandlers = () => ({
  onSetStatus: vi.fn(),
  onSetIndexing: vi.fn(),
  onAddTag: vi.fn(),
  onRemoveTag: vi.fn(),
});

test("shows Accept and Reject buttons", () => {
  render(<DetectorImageCard exposure={makeExposure({})} {...noopHandlers()} />);
  expect(screen.getByRole("button", { name: /accept/i })).toBeInTheDocument();
  expect(screen.getByRole("button", { name: /reject/i })).toBeInTheDocument();
});

test("Reject → Other reveals custom reason input", async () => {
  render(<DetectorImageCard exposure={makeExposure({})} {...noopHandlers()} />);
  await userEvent.click(screen.getByRole("button", { name: /reject/i }));
  // After clicking Reject, quick-pick chips appear (Flare, Other) — no input yet.
  expect(screen.queryByPlaceholderText(/reason/i)).not.toBeInTheDocument();
  await userEvent.click(screen.getByRole("button", { name: /other/i }));
  expect(screen.getByPlaceholderText(/reason/i)).toBeInTheDocument();
});

test("Reject → Flare immediately rejects with reason tag", async () => {
  const handlers = noopHandlers();
  render(<DetectorImageCard exposure={makeExposure({})} {...handlers} />);
  await userEvent.click(screen.getByRole("button", { name: /reject/i }));
  await userEvent.click(screen.getByRole("button", { name: /flare/i }));
  expect(handlers.onSetStatus).toHaveBeenCalledWith("rejected");
  expect(handlers.onAddTag).toHaveBeenCalledWith("rejection_reason", "Flare");
});

test("Re-rejecting deletes prior rejection_reason tags before adding the new one", async () => {
  // Regression: PR #9 review issue 1 — rejection_reason tags accumulated
  // across reject → un-reject → re-reject because submitReject only added.
  // Now it removes any pre-existing rejection_reason tags first.
  //
  // Realistic flow: the exposure was rejected once with reason "Flare".
  // Backend round-trip leaves the tag in place. The user un-rejects
  // (status: null, tag remains because we don't auto-clean the tag), then
  // re-rejects with a new reason. Without dedup, the backend ends up with
  // two `rejection_reason` tags and the UI silently surfaces only the first.
  const handlers = noopHandlers();
  const exposure = makeExposure({
    // status: null after the un-reject step; the prior tag is still present
    status: null,
    tags: [
      { id: 11, exposure_id: 1, key: "rejection_reason", value: "Flare", source: "manual" },
      { id: 12, exposure_id: 1, key: "rejection_reason", value: "stale", source: "manual" },
      { id: 13, exposure_id: 1, key: "other_label",      value: "keep me", source: "manual" },
    ],
  });
  render(<DetectorImageCard exposure={exposure} {...handlers} />);
  // Reject again → picking → Flare
  await userEvent.click(screen.getByRole("button", { name: /reject/i }));
  await userEvent.click(screen.getByRole("button", { name: /flare/i }));
  // Both prior rejection_reason tags should be removed; the unrelated tag stays.
  expect(handlers.onRemoveTag).toHaveBeenCalledWith(11);
  expect(handlers.onRemoveTag).toHaveBeenCalledWith(12);
  expect(handlers.onRemoveTag).not.toHaveBeenCalledWith(13);
  // And the new tag is added after the dedup.
  expect(handlers.onAddTag).toHaveBeenCalledWith("rejection_reason", "Flare");
});

test("Use for indexing is disabled when rejected", () => {
  render(
    <DetectorImageCard
      exposure={makeExposure({ status: "rejected" })}
      {...noopHandlers()}
    />,
  );
  expect(screen.getByRole("button", { name: /indexing/i })).toBeDisabled();
});

test("Accept calls onSetStatus with 'accepted'", async () => {
  const handlers = noopHandlers();
  render(<DetectorImageCard exposure={makeExposure({})} {...handlers} />);
  await userEvent.click(screen.getByRole("button", { name: /accept/i }));
  expect(handlers.onSetStatus).toHaveBeenCalledWith("accepted");
});
