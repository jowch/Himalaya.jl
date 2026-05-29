import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { act, render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Mutation } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { InfrastructureBanner } from "../src/components/InfrastructureBanner";

function withQC(qc: QueryClient) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

interface AddOpts {
  status?: "pending" | "success" | "error" | "idle";
  submittedAt?: number;
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
    options: {},
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

  it("upgrades to stuck state after 30s", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 31000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "stuck");
    expect(banner).toHaveTextContent(/Couldn.t save/);
    expect(screen.getByRole("button", { name: /refresh/i })).toBeInTheDocument();
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

  it("labels the stuck state with an Error severity word + icon and corrected prose", () => {
    addMutation(qc, { status: "pending", submittedAt: Date.now() - 31000 });
    render(<InfrastructureBanner />, { wrapper: withQC(qc) });
    const banner = screen.getByTestId("infrastructure-banner");
    expect(banner).toHaveAttribute("data-state", "stuck");
    expect(banner.querySelector('[aria-label="Error"]')).not.toBeNull();
    // Em-dash retired: two sentences, no " — ".
    expect(banner).toHaveTextContent("Couldn’t save. Try refreshing.");
    expect(banner.textContent ?? "").not.toContain("—");
    expect(screen.getByRole("button", { name: /refresh/i })).toBeInTheDocument();
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
