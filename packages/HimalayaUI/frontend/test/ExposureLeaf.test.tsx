import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ExposureLeaf } from "../src/print/components/ExposureLeaf";
import type { LoadExposure } from "../src/api";

// §8.8: LoadExposure = { id, filename, horizontal_position, timestamp } — NO frame_no/status.
const exp: LoadExposure = {
  id: 100, filename: "HA_0241_001.tif", horizontal_position: 8.4, timestamp: "10:02:00",
};

describe("ExposureLeaf", () => {
  it("shows filename, H-position, and time", () => {
    render(<ExposureLeaf exposure={exp} thumbSrc={null} onMove={() => {}} />);
    expect(screen.getByText("HA_0241_001.tif")).toBeInTheDocument();
    expect(screen.getByText(/8\.4/)).toBeInTheDocument();
    expect(screen.getByText(/10:02/)).toBeInTheDocument();
  });
  it("Move button calls onMove with the exposure id and the button element as anchor", () => {
    const onMove = vi.fn();
    render(<ExposureLeaf exposure={exp} thumbSrc={null} onMove={onMove} />);
    const btn = screen.getByRole("button", { name: /move/i });
    fireEvent.click(btn);
    expect(onMove).toHaveBeenCalledTimes(1);
    const [calledId, calledEl] = onMove.mock.calls[0]!;
    expect(calledId).toBe(100);
    // anchorEl must be the button element itself
    expect(calledEl).toBe(btn);
  });
  it("renders an em-dash when horizontal_position is null", () => {
    render(<ExposureLeaf exposure={{ ...exp, horizontal_position: null }} thumbSrc={null} onMove={() => {}} />);
    expect(screen.getByTestId("exposure-leaf").textContent).toContain("—");
  });
});
