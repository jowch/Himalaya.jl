import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { CompareStatusSurface } from "../src/components/CompareStatusSurface";

describe("CompareStatusSurface — Compare UX C-9", () => {
  it("renders nothing when there are no banners", () => {
    const { container } = render(
      <CompareStatusSurface
        needsReview={null} serverUpdate={null} savedAt={null}
      />,
    );
    expect(container.querySelector("[data-testid='compare-status-surface']")).toBeNull();
  });

  it("renders a needs-review banner with action", () => {
    const onReanalyze = vi.fn();
    render(<CompareStatusSurface
      needsReview={{ onReanalyze }}
      serverUpdate={null}
      savedAt={null}
    />);
    expect(screen.getByText(/Needs review/)).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("compare-status-resnapshot"));
    expect(onReanalyze).toHaveBeenCalled();
  });

  it("renders a server-update banner with acknowledge action", () => {
    const onAcknowledge = vi.fn();
    render(<CompareStatusSurface
      needsReview={null}
      serverUpdate={{ previousHash: "h_prev", onAcknowledge }}
      savedAt={null}
    />);
    expect(screen.getByText(/updated since you last viewed/i)).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("compare-status-acknowledge"));
    expect(onAcknowledge).toHaveBeenCalled();
  });

  it("renders 'Saved' after a recent save", () => {
    render(<CompareStatusSurface needsReview={null} serverUpdate={null} savedAt={Date.now()}/>);
    expect(screen.getByText("Saved")).toBeInTheDocument();
  });

  it("hides the Saved pill after 2s", async () => {
    vi.useFakeTimers();
    const savedAt = Date.now();
    render(<CompareStatusSurface needsReview={null} serverUpdate={null} savedAt={savedAt} />);
    expect(screen.getByText("Saved")).toBeInTheDocument();
    // Flush the retire-after-2s setTimeout's state update inside act().
    act(() => { vi.advanceTimersByTime(2100); });
    expect(screen.queryByText("Saved")).toBeNull();
    vi.useRealTimers();
  });
});
