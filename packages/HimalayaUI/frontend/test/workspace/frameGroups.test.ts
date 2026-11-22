import { describe, it, expect } from "vitest";
import type { CorpusSample, Exposure } from "../../src/api";
import { allFrames, frameGroups, frameOf } from "../../src/workspace/frameGroups";

const s = (id: number, name: string, setup: number) => ({ id, name, setup_id: setup }) as CorpusSample;
const f = (id: number) => ({ id }) as Exposure;
const rank = (x: CorpusSample) => x.setup_id ?? 0;
const title = (x: CorpusSample) => `${x.name}/${x.setup_id}`;

describe("frameGroups", () => {
  it("puts the open sample first, then mates by setup order, and names groups only when there are mates", () => {
    const groups = frameGroups(
      { sample: s(1, "A", 2), frames: [f(10)] },
      [{ sample: s(3, "C", 3), frames: [f(30)] }, { sample: s(2, "B", 1), frames: [f(20), f(21)] }],
      rank,
      title,
    );
    expect(groups.map((g) => [g.sampleId, g.title])).toEqual([[1, "A/2"], [2, "B/1"], [3, "C/3"]]);
    expect(allFrames(groups).map((e) => e.id)).toEqual([10, 20, 21, 30]);
    expect(frameOf(groups, 21)).toEqual({ group: groups[1], index: 1 });
    expect(frameOf(groups, 99)).toBeUndefined();
  });

  it("leaves out mates whose frames have not loaded, and a lone sample has no title", () => {
    const groups = frameGroups({ sample: s(1, "A", 1), frames: [f(10)] }, [{ sample: s(2, "B", 1), frames: [] }], rank, title);
    expect(groups).toEqual([{ sampleId: 1, title: undefined, exposures: [f(10)] }]);
  });
});
