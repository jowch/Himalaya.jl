import { describe, it, expect, vi, afterEach, beforeEach } from "vitest";
import { act, render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider, useMutation } from "@tanstack/react-query";
import { useEffect } from "react";
import { InfrastructureBanner } from "../src/print/shell/InfrastructureBanner";
import { useCheckForNewFilesOnOpen, useTriggerScan } from "../src/queries";
import { triggerScan } from "../src/api";

vi.mock("../src/api", async (orig) => ({ ...(await orig<typeof import("../src/api")>()), triggerScan: vi.fn(() => new Promise(() => {})) }));

function OpenWorkspace(): null {
  useCheckForNewFilesOnOpen(7);
  return null;
}

function Rescan(): null {
  const { mutate } = useTriggerScan(7);
  useEffect(() => mutate(), [mutate]);
  return null;
}

function SaveSomething(): null {
  const { mutate } = useMutation({ mutationFn: () => triggerScan(7) });
  useEffect(() => mutate(), [mutate]);
  return null;
}

async function bannerAfterTwoSeconds(Trigger: () => null): Promise<HTMLElement | null> {
  render(
    <QueryClientProvider client={new QueryClient()}>
      <Trigger />
      <InfrastructureBanner />
    </QueryClientProvider>,
  );
  await act(async () => { await vi.advanceTimersByTimeAsync(1000); });
  act(() => { vi.advanceTimersByTime(1000); });
  expect(triggerScan).toHaveBeenCalled();
  return screen.queryByTestId("infrastructure-banner");
}

describe("a scan in flight", () => {
  beforeEach(() => { vi.useFakeTimers(); vi.mocked(triggerScan).mockClear(); });
  afterEach(() => vi.useRealTimers());

  it("shows the saving banner for an ordinary pending request", async () => {
    expect(await bannerAfterTwoSeconds(SaveSomething)).toBeInTheDocument();
  });

  it.each([
    ["the check for new files when a workspace opens", OpenWorkspace],
    ["a rescan", Rescan],
  ])("never shows the saving banner during %s", async (_name, Trigger) => {
    expect(await bannerAfterTwoSeconds(Trigger)).not.toBeInTheDocument();
  });
});
