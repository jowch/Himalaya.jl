import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { act, render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Mutation } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { InfrastructureBanner } from "../src/shell/InfrastructureBanner";
import { useConnection } from "../src/lib/connection";

function withQC(qc: QueryClient) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

interface AddOpts {
  status?: "pending" | "success" | "error" | "idle";
  submittedAt?: number;
  meta?: Record<string, unknown>;
}

function addMutation(qc: QueryClient, opts: AddOpts = {}): Mutation {
  const m = {
    state: {
      status: opts.status ?? "pending",
      submittedAt: opts.submittedAt ?? Date.now(),
      context: {},
      variables: {},
      data: undefined,
      error: null,
      failureCount: 0,
      failureReason: null,
      isPaused: false,
    },
    options: opts.meta ? { meta: opts.meta } : {},
  } as unknown as Mutation;
  qc.getMutationCache().add(m);
  return m;
}

describe("InfrastructureBanner", () => {
  let qc: QueryClient;

  beforeEach(() => {
    vi.useFakeTimers();
    vi.setSystemTime(new Date("2026-01-01T00:00:00Z"));
    qc = new QueryClient();
  });
  afterEach(() => {
    vi.useRealTimers();
    useConnection.setState({ streamDownSince: null });
  });

  it("renders nothing when no pending mutations", () => {
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
  });

  it("does not show banner when pending mutation is <500ms old", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 100 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
  });

  it("shows banner once pending mutation has been in-flight >500ms", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 100 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
    act(() => {
      vi.advanceTimersByTime(1000); // tick past 500ms threshold
    });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "showing");
    expect(banner).toHaveTextContent(/Saving your changes/);
  });

  it("ignores a pending request that saves nothing, such as the check for new files when a workspace opens", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 2000, meta: { savesWork: false } });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    act(() => {
      vi.advanceTimersByTime(1000);
    });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
  });

  it("after 30s says how many changes are not saved yet and that it keeps trying, without asking for a refresh", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 31000 });
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 2000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "stuck");
    expect(banner).toHaveTextContent("2 changes not saved yet. Still trying.");
    expect(screen.queryByRole("button", { name: /refresh/i })).not.toBeInTheDocument();
  });

  it("says the connection is lost once the live stream has been down for 5s, with no edits pending", () => {
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    act(() => { useConnection.getState().setStreamDown(); });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
    act(() => { vi.advanceTimersByTime(5500); });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "offline");
    expect(banner).toHaveTextContent("Connection lost. Retrying…");
    act(() => { useConnection.getState().setStreamOpen(); });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
  });

  it("while the connection is lost, counts the edits that will save when it is back", () => {
    act(() => { useConnection.setState({ streamDownSince: Date.now() - 10_000 }); });
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 40_000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "offline");
    expect(banner).toHaveTextContent("Connection lost. Retrying… 1 change will save when it’s back.");
  });

  it("uses a full hairline border, not a left-edge severity stripe (showing)", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 1000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner.className).toContain("border-hair");
    expect(banner.className).not.toContain("border-l-4");
    expect(banner.className).not.toContain("border-warning");
  });

  it("labels the showing state with a Saving severity word + icon", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 1000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "showing");
    expect(banner).toHaveTextContent("Saving");
    expect(banner.querySelector('[aria-label="Saving"]')).not.toBeNull();
  });

  it("labels the stuck state with a Not saved icon and no em-dash", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 31000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "stuck");
    expect(banner.querySelector('[aria-label="Not saved"]')).not.toBeNull();
    expect(banner).toHaveTextContent("1 change not saved yet. Still trying.");
    expect(banner.textContent ?? "").not.toContain("—");
  });

  it("disappears when the mutation settles", () => {
    const m = addMutation(qc, { status: "pending", submittedAt: Date.now() - 1000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    expect(screen.getByTestId("infrastructure-banner")).toBeInTheDocument();
    act(() => {
      // Mark mutation as settled and remove it from the cache.
      (m.state as { status: string }).status = "success";
      qc.getMutationCache().remove(m);
      vi.advanceTimersByTime(1000);
    });
    expect(screen.queryByTestId("infrastructure-banner")).not.toBeInTheDocument();
  });
});
