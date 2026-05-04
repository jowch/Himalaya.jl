import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useAddExposureTag,
  useRemoveExposureTag,
  useSetExposureStatus,
  useSelectExposure,
  queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import type { Exposure } from "../src/api";

const SAMPLE_ID = 10;
const EXPOSURE_ID = 100;
const OTHER_EXPOSURE_ID = 101;

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

function makeExposure(over: Partial<Exposure>): Exposure {
  return {
    id: EXPOSURE_ID,
    sample_id: SAMPLE_ID,
    filename: "frame_001.tif",
    kind: "file",
    selected: false,
    status: null,
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
    ...over,
  };
}

describe("queries — exposure mutations (queue-driven)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  // -------------------------------------------------------------------------
  // useAddExposureTag
  // -------------------------------------------------------------------------

  describe("useAddExposureTag", () => {
    it("appends an optimistic placeholder synchronously", () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({})]);
      // Stub fetch with a never-resolving promise so optimistic state is observed.
      vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
      const { result } = renderHook(
        () => useAddExposureTag(SAMPLE_ID, EXPOSURE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ key: "buffer", value: "PBS" }); });
      const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
      expect(list[0].tags).toHaveLength(1);
      expect(list[0].tags[0].id).toBeLessThan(0);
      expect(list[0].tags[0].key).toBe("buffer");
      expect(list[0].tags[0].value).toBe("PBS");
    });

    it("replaces the optimistic placeholder with the server tag on success", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({})]);
      mockOnce(201, { id: 7, key: "buffer", value: "PBS", source: "manual" });
      const { result } = renderHook(
        () => useAddExposureTag(SAMPLE_ID, EXPOSURE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ key: "buffer", value: "PBS" }); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list[0].tags).toHaveLength(1);
        expect(list[0].tags[0].id).toBe(7);
        expect(list[0].tags[0].key).toBe("buffer");
      });
    });

    it("rolls back to the seeded state on 4xx", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({})]);
      mockOnce(400, { error: "bad" });
      const { result } = renderHook(
        () => useAddExposureTag(SAMPLE_ID, EXPOSURE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ key: "buffer", value: "PBS" }); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list[0].tags).toHaveLength(0);
      });
    });
  });

  // -------------------------------------------------------------------------
  // useRemoveExposureTag
  // -------------------------------------------------------------------------

  describe("useRemoveExposureTag", () => {
    it("removes the tag optimistically", () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [
        makeExposure({ tags: [{ id: 7, key: "buffer", value: "PBS", source: "manual" }] }),
      ]);
      vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
      const { result } = renderHook(
        () => useRemoveExposureTag(SAMPLE_ID, EXPOSURE_ID), { wrapper },
      );
      act(() => { result.current.mutate(7); });
      const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
      expect(list[0].tags).toHaveLength(0);
    });

    it("keeps the tag removed after a 204 success", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [
        makeExposure({ tags: [{ id: 7, key: "buffer", value: "PBS", source: "manual" }] }),
      ]);
      mockOnce(204, null);
      const { result } = renderHook(
        () => useRemoveExposureTag(SAMPLE_ID, EXPOSURE_ID), { wrapper },
      );
      act(() => { result.current.mutate(7); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list[0].tags).toHaveLength(0);
      });
    });

    it("restores the tag on 4xx", async () => {
      const { client, wrapper } = withClient();
      const seeded = [
        makeExposure({ tags: [{ id: 7, key: "buffer", value: "PBS", source: "manual" }] }),
      ];
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), seeded);
      mockOnce(404, { error: "missing" });
      const { result } = renderHook(
        () => useRemoveExposureTag(SAMPLE_ID, EXPOSURE_ID), { wrapper },
      );
      act(() => { result.current.mutate(7); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list[0].tags).toHaveLength(1);
        expect(list[0].tags[0].id).toBe(7);
      });
    });
  });

  // -------------------------------------------------------------------------
  // useSetExposureStatus
  // -------------------------------------------------------------------------

  describe("useSetExposureStatus", () => {
    it("flips status optimistically in the exposures list", () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({})]);
      vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
      const { result } = renderHook(
        () => useSetExposureStatus(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ exposureId: EXPOSURE_ID, status: "rejected" }); });
      const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
      expect(list[0].status).toBe("rejected");
    });

    it("also updates the single-exposure cache when present", () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({})]);
      client.setQueryData(queryKeys.exposure(EXPOSURE_ID), makeExposure({}));
      vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
      const { result } = renderHook(
        () => useSetExposureStatus(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate({ exposureId: EXPOSURE_ID, status: "accepted" }); });
      const single = client.getQueryData(queryKeys.exposure(EXPOSURE_ID)) as Exposure;
      expect(single.status).toBe("accepted");
    });

    it("retains the optimistic status after success (server response is a no-op)", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({})]);
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
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [makeExposure({ status: null })]);
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

  // -------------------------------------------------------------------------
  // useSelectExposure
  // -------------------------------------------------------------------------

  describe("useSelectExposure", () => {
    it("flips selected to the targeted exposure optimistically", () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [
        makeExposure({ id: EXPOSURE_ID, selected: true }),
        makeExposure({ id: OTHER_EXPOSURE_ID, selected: false }),
      ]);
      vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
      const { result } = renderHook(
        () => useSelectExposure(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate(OTHER_EXPOSURE_ID); });
      const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
      expect(list.find((e) => e.id === EXPOSURE_ID)!.selected).toBe(false);
      expect(list.find((e) => e.id === OTHER_EXPOSURE_ID)!.selected).toBe(true);
    });

    it("retains the selection after success", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [
        makeExposure({ id: EXPOSURE_ID, selected: true }),
        makeExposure({ id: OTHER_EXPOSURE_ID, selected: false }),
      ]);
      mockOnce(200, { id: OTHER_EXPOSURE_ID, selected: true });
      const { result } = renderHook(
        () => useSelectExposure(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate(OTHER_EXPOSURE_ID); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list.find((e) => e.id === OTHER_EXPOSURE_ID)!.selected).toBe(true);
      });
    });

    it("rolls back selection on 4xx", async () => {
      const { client, wrapper } = withClient();
      client.setQueryData(queryKeys.exposures(SAMPLE_ID), [
        makeExposure({ id: EXPOSURE_ID, selected: true }),
        makeExposure({ id: OTHER_EXPOSURE_ID, selected: false }),
      ]);
      mockOnce(409, { error: "conflict" });
      const { result } = renderHook(
        () => useSelectExposure(SAMPLE_ID), { wrapper },
      );
      act(() => { result.current.mutate(OTHER_EXPOSURE_ID); });
      await waitFor(() => {
        const list = client.getQueryData(queryKeys.exposures(SAMPLE_ID)) as Exposure[];
        expect(list.find((e) => e.id === EXPOSURE_ID)!.selected).toBe(true);
        expect(list.find((e) => e.id === OTHER_EXPOSURE_ID)!.selected).toBe(false);
      });
    });
  });
});
