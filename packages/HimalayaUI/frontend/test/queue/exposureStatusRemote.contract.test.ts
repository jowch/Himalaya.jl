import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import type { Exposure } from "../../src/api";
import { queryKeys } from "../../src/queries";

function exposure(id: number, status: Exposure["status"] = null): Exposure {
  return { id, status } as Exposure;
}

describe("a foreign set_exposure_status reaches the sample's exposure list", () => {
  it("patches the list holding the exposure and leaves other samples alone", () => {
    const qc = new QueryClient();
    const mine = [exposure(5), exposure(6)];
    const other = [exposure(7)];
    qc.setQueryData(queryKeys.exposures(1), mine);
    qc.setQueryData(queryKeys.exposures(2), other);

    applyRemoteToCache(
      { id: 9, kind: "set_exposure_status", entity_type: "exposure", entity_id: 5, payload: { status: "rejected" } },
      qc,
    );

    expect(qc.getQueryData<Exposure[]>(queryKeys.exposures(1))?.map((e) => e.status)).toEqual(["rejected", null]);
    expect(qc.getQueryData(queryKeys.exposures(2))).toBe(other);
  });

  it("restores a dropped exposure", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.exposures(1), [exposure(5, "rejected")]);
    applyRemoteToCache(
      { id: 10, kind: "set_exposure_status", entity_type: "exposure", entity_id: 5, payload: { status: null } },
      qc,
    );
    expect(qc.getQueryData<Exposure[]>(queryKeys.exposures(1))?.[0]?.status).toBeNull();
  });
});
