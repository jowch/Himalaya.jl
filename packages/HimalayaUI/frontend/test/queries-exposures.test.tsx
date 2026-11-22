import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useSetExposureStatus,
  queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import type { Exposure } from "../src/api";
import { makeExposure } from "./helpers/exposure";

const SAMPLE_ID = 10;
const EXPOSURE_ID = 100;

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(status === 204 ? null : JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

function frame(over: Partial<Exposure>): Exposure {
  return makeExposure({ id: EXPOSURE_ID, sample_id: SAMPLE_ID, filename: "frame_001.tif", ...over });
}

describe("queries — exposure mutations (queue-driven)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  // -------------------------------------------------------------------------
  // useSetExposureStatus
  // -------------------------------------------------------------------------

  describe("useSetExposureStatus", () => {
    it("flips status optimistically in the exposures list", () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [frame({})]);
      vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
      const { result } = renderHook(
        () => useSetExposureStatus(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ exposureId: EXPOSURE_ID, status: "rejected" }); });
      const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
      expect(list[0].status).toBe("rejected");
    });

    it("retains the optimistic status after success (server response is a no-op)", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [frame({})]);
      mockOnce(200, { id: EXPOSURE_ID, status: "rejected" });
      const { result } = renderHook(
        () => useSetExposureStatus(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ exposureId: EXPOSURE_ID, status: "rejected" }); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list[0].status).toBe("rejected");
      });
    });

    it("rolls back status on 4xx", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [frame({ status: null })]);
      mockOnce(422, { error: "invalid" });
      const { result } = renderHook(
        () => useSetExposureStatus(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ exposureId: EXPOSURE_ID, status: "rejected" }); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list[0].status).toBeNull();
      });
    });
  });
});
