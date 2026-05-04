/**
 * Rollback symmetry test.
 *
 * Invariant: for every mutator, calling `restore()` from the RollbackContext
 * returned by `onMutate` MUST return the cache to bit-exact prior state.
 *
 * Why: a partial-rollback bug leaves the cache in a fictional state — neither
 * the user's optimistic intent nor the server-confirmed state — which surfaces
 * later as flicker, ghost rows, or "the row is back even after I undid it."
 * Three rollback bugs of this shape have been fixed in this PR (forgotten
 * variant query key, missing optional snapshot, leaked partial state).
 *
 * Method: snapshot every cache key the mutator might touch (deep-cloned),
 * run onMutate, run restore, deep-equal-assert the cache against the snapshot.
 *
 * Discipline: when adding a new mutator, add a row to MUTATORS below. The
 * structure forces you to enumerate the keys the mutator touches; if your
 * onMutate touches a key not in `cacheKeysTouched`, that is itself a bug
 * (the rollback can't restore what it didn't snapshot).
 */
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator, removeIndexFromGroupMutator, deleteIndexMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import { reanalyzeExposureMutator } from "../../src/lib/queue/mutators/reanalyzeExposure";
import {
  updateSampleMutator,
  addSampleTagMutator, removeSampleTagMutator,
  addExposureTagMutator, removeExposureTagMutator,
  postSampleMessageMutator,
  setExposureStatusMutator, selectExposureMutator,
} from "../../src/lib/queue/mutators/trivial";
import { queryKeys } from "../../src/queries";
import type { QueryKey } from "@tanstack/react-query";

// Fixtures used to populate cache before onMutate runs.
const PEAK = {
  id: 7, exposure_id: 5, q: 0.5, intensity: 1.0, prominence: 0.8,
  sharpness: 30, source: "auto" as const, excluded: false,
};
const EXPOSURE = {
  id: 5, sample_id: 1, filename: null, kind: "file" as const, selected: true,
  status: null, image_path: null, image_version: "",
  trace_hash: null, analysis_inputs_hash: "h0", tags: [], sources: [],
};
const SAMPLE = {
  id: 10, experiment_id: 1, label: "D1", name: "n", notes: null,
  tags: [{ id: 1, key: "k", value: "v", source: "manual" }],
};
const GROUP = {
  id: 1, exposure_id: 5, kind: "custom" as const, active: true, members: [10] as number[],
};
const INDEX = {
  id: 10, exposure_id: 5, phase: "Pn3m", basis: 0.1, score: 0.9,
  r_squared: 0.99, lattice_d: 50, ngc: 0.5, status: "candidate" as const,
  kind: "auto" as const, inputs_hash: "h1", peaks: [], predicted_q: [0.1],
};
const EXP_LIST_KEY = ["sample", 1, "exposures", { excludeRejected: false }] as const;

interface Spec {
  name: string;
  // Keys to snapshot (deep-clone) before onMutate.
  keys: QueryKey[];
  // Seed the cache. Caller writes whatever is necessary so prev !== undefined.
  seed: (qc: QueryClient) => void;
  // Mutator + flat payload merged exactly as useQueueMutation would build it.
  run: (qc: QueryClient) => void;
}

const SPECS: Spec[] = [
  {
    name: "peakAdd",
    keys: [queryKeys.peaks(5)],
    seed: (qc) => qc.setQueryData(queryKeys.peaks(5), [PEAK]),
    run: (qc) => {
      const ctx = peakAddMutator.onMutate({
        kind: "peak_added", clientOpId: "op", payload: { q: 0.42 },
        exposureId: 5, username: "alice", clientId: "tab",
        q: 0.42,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "peakRemove",
    keys: [queryKeys.peaks(5)],
    seed: (qc) => qc.setQueryData(queryKeys.peaks(5), [PEAK]),
    run: (qc) => {
      const ctx = peakRemoveMutator.onMutate({
        kind: "peak_removed", clientOpId: "op", payload: { peakId: 7 },
        exposureId: 5, username: "alice", clientId: "tab",
        peakId: 7,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "peakExclude",
    keys: [queryKeys.peaks(5)],
    seed: (qc) => qc.setQueryData(queryKeys.peaks(5), [PEAK]),
    run: (qc) => {
      const ctx = peakExcludeMutator.onMutate({
        kind: "peak_excluded", clientOpId: "op",
        payload: { peakId: 7, q: 0.5 },
        exposureId: 5, username: "alice", clientId: "tab",
        peakId: 7, q: 0.5,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "peakUnexclude",
    keys: [queryKeys.peaks(5)],
    seed: (qc) => qc.setQueryData(queryKeys.peaks(5), [{ ...PEAK, excluded: true }]),
    run: (qc) => {
      const ctx = peakUnexcludeMutator.onMutate({
        kind: "peak_unexcluded", clientOpId: "op",
        payload: { peakId: 7, q: 0.5 },
        exposureId: 5, username: "alice", clientId: "tab",
        peakId: 7, q: 0.5,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "addIndexToGroup",
    keys: [queryKeys.groups(5)],
    seed: (qc) => qc.setQueryData(queryKeys.groups(5), [GROUP]),
    run: (qc) => {
      const ctx = addIndexToGroupMutator.onMutate({
        kind: "index_confirmed", clientOpId: "op",
        payload: { groupId: 1, indexId: 99 },
        exposureId: 5, groupId: 1, username: "alice", clientId: "tab",
        indexId: 99,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "removeIndexFromGroup",
    keys: [queryKeys.groups(5)],
    seed: (qc) => qc.setQueryData(queryKeys.groups(5), [GROUP]),
    run: (qc) => {
      const ctx = removeIndexFromGroupMutator.onMutate({
        kind: "index_unconfirmed", clientOpId: "op",
        payload: { groupId: 1, indexId: 10 },
        exposureId: 5, groupId: 1, username: "alice", clientId: "tab",
        indexId: 10,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "deleteIndex",
    keys: [queryKeys.indices(5), queryKeys.groups(5)],
    seed: (qc) => {
      qc.setQueryData(queryKeys.indices(5), [INDEX]);
      qc.setQueryData(queryKeys.groups(5), [GROUP]);
    },
    run: (qc) => {
      const ctx = deleteIndexMutator.onMutate({
        kind: "delete_index", clientOpId: "op",
        payload: { indexId: 10 },
        exposureId: 5, username: "alice", clientId: "tab",
        indexId: 10,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "createSpeculative",
    keys: [queryKeys.indices(5), queryKeys.groups(5)],
    seed: (qc) => {
      qc.setQueryData(queryKeys.indices(5), [INDEX]);
      qc.setQueryData(queryKeys.groups(5), [GROUP]);
    },
    run: (qc) => {
      const ctx = createSpeculativeMutator.onMutate({
        kind: "speculative_created", clientOpId: "op",
        payload: { phase: "Pn3m", anchor_peak_id: 7, anchor_ratio: 1, additional: [] },
        exposureId: 5, username: "alice", clientId: "tab",
        phase: "Pn3m", anchor_peak_id: 7, anchor_ratio: 1, additional: [],
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "reanalyzeExposure",
    keys: [queryKeys.exposure(5), queryKeys.peaks(5), queryKeys.indices(5)],
    seed: (qc) => {
      qc.setQueryData(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData(queryKeys.peaks(5), [PEAK]);
      qc.setQueryData(queryKeys.indices(5), [INDEX]);
    },
    run: (qc) => {
      const ctx = reanalyzeExposureMutator.onMutate({
        kind: "reanalyze_exposure", clientOpId: "op", payload: {},
        exposureId: 5, username: "alice", clientId: "tab",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "updateSample",
    keys: [queryKeys.samples(1), queryKeys.sample(10)],
    seed: (qc) => {
      qc.setQueryData(queryKeys.samples(1), [SAMPLE]);
      qc.setQueryData(queryKeys.sample(10), SAMPLE);
    },
    run: (qc) => {
      const ctx = updateSampleMutator.onMutate({
        kind: "update_sample", clientOpId: "op",
        payload: { sampleId: 10 },
        sampleId: 10, experimentId: 1, username: "alice", clientId: "tab",
        name: "renamed",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "addSampleTag",
    keys: [queryKeys.samples(1)],
    seed: (qc) => qc.setQueryData(queryKeys.samples(1), [SAMPLE]),
    run: (qc) => {
      const ctx = addSampleTagMutator.onMutate({
        kind: "add_tag", clientOpId: "op",
        payload: { key: "x", value: "y" },
        sampleId: 10, experimentId: 1, username: "alice", clientId: "tab",
        key: "x", value: "y",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "removeSampleTag",
    keys: [queryKeys.samples(1)],
    seed: (qc) => qc.setQueryData(queryKeys.samples(1), [SAMPLE]),
    run: (qc) => {
      const ctx = removeSampleTagMutator.onMutate({
        kind: "remove_tag", clientOpId: "op",
        payload: { tagId: 1 },
        sampleId: 10, experimentId: 1, username: "alice", clientId: "tab",
        tagId: 1,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "addExposureTag",
    keys: [EXP_LIST_KEY as unknown as QueryKey],
    seed: (qc) => qc.setQueryData(EXP_LIST_KEY as unknown as QueryKey, [EXPOSURE]),
    run: (qc) => {
      const ctx = addExposureTagMutator.onMutate({
        kind: "add_tag", clientOpId: "op",
        payload: { key: "x", value: "y" },
        sampleId: 1, exposureId: 5, username: "alice", clientId: "tab",
        key: "x", value: "y",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "removeExposureTag",
    keys: [EXP_LIST_KEY as unknown as QueryKey],
    seed: (qc) => qc.setQueryData(
      EXP_LIST_KEY as unknown as QueryKey,
      [{ ...EXPOSURE, tags: [{ id: 99, key: "k", value: "v", source: "manual" }] }],
    ),
    run: (qc) => {
      const ctx = removeExposureTagMutator.onMutate({
        kind: "remove_tag", clientOpId: "op",
        payload: { tagId: 99 },
        sampleId: 1, exposureId: 5, username: "alice", clientId: "tab",
        tagId: 99,
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "postSampleMessage",
    keys: [queryKeys.messages(10)],
    seed: (qc) => qc.setQueryData(queryKeys.messages(10), [
      { id: 1, sample_id: 10, author_id: 3, author: "alice",
        body: "hello", created_at: "2026-05-03" },
    ]),
    run: (qc) => {
      const ctx = postSampleMessageMutator.onMutate({
        kind: "post_message", clientOpId: "op",
        payload: { body: "world" },
        sampleId: 10, username: "alice", clientId: "tab",
        body: "world",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "setExposureStatus",
    keys: [EXP_LIST_KEY as unknown as QueryKey, queryKeys.exposure(5)],
    seed: (qc) => {
      qc.setQueryData(EXP_LIST_KEY as unknown as QueryKey, [EXPOSURE]);
      qc.setQueryData(queryKeys.exposure(5), EXPOSURE);
    },
    run: (qc) => {
      const ctx = setExposureStatusMutator.onMutate({
        kind: "set_exposure_status", clientOpId: "op",
        payload: { exposureId: 5, status: "rejected" },
        sampleId: 1, username: "alice", clientId: "tab",
        exposureId: 5, status: "rejected",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "selectExposure",
    keys: [EXP_LIST_KEY as unknown as QueryKey],
    seed: (qc) => qc.setQueryData(EXP_LIST_KEY as unknown as QueryKey,
      [EXPOSURE, { ...EXPOSURE, id: 6, selected: false }]),
    run: (qc) => {
      const ctx = selectExposureMutator.onMutate({
        kind: "select_exposure", clientOpId: "op",
        payload: { exposureId: 6 },
        sampleId: 1, username: "alice", clientId: "tab",
        exposureId: 6,
      } as any, qc);
      ctx.restore();
    },
  },
];

describe("Rollback symmetry — restore() returns cache to prior state for every mutator", () => {
  for (const spec of SPECS) {
    it(`${spec.name}: snapshot → onMutate → restore → deep-equal prior`, () => {
      const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
      spec.seed(qc);
      // Deep-clone snapshot — JSON round-trip is fine for our cache shapes
      // (no Dates, no functions, no circular refs).
      const before = spec.keys.map((k) => JSON.parse(JSON.stringify(qc.getQueryData(k))));
      spec.run(qc);
      for (const [i, k] of spec.keys.entries()) {
        const after = qc.getQueryData(k);
        expect({ key: spec.keys[i], state: after }).toEqual(
          { key: spec.keys[i], state: before[i] });
      }
    });
  }
});
