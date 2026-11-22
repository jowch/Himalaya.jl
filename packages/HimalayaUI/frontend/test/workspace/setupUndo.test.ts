import { describe, expect, it, vi } from "vitest";
import type { Experiment, Setup } from "../../src/api";
import type { EditOps } from "../../src/queries";
import { setupEntry } from "../../src/workspace/undoEntries";

function opsWith(setup: Partial<Setup>): { ops: EditOps; patchSetup: ReturnType<typeof vi.fn> } {
  const patchSetup = vi.fn();
  const s = { id: 7, folder: "setup1", display_name: null, name: "setup1", included: true, ...setup } as Setup;
  const ops = { patchSetup, read: { experiment: async () => ({ id: 3, setups: [s] }) as Experiment } } as unknown as EditOps;
  return { ops, patchSetup };
}

const where = { label: "x", sampleId: 1, exposureId: undefined, setupId: 7 };

describe("setup undo", () => {
  it("undoes a hide and a rename when the server still shows the edit", async () => {
    const hidden = opsWith({ included: false });
    expect(await setupEntry(hidden.ops, { ...where, prev: { included: true }, next: { included: false } }).undo()).toBe(1);
    expect(hidden.patchSetup).toHaveBeenCalledWith(7, { included: true });

    const renamed = opsWith({ display_name: "Long pipe", name: "Long pipe" });
    expect(await setupEntry(renamed.ops, { ...where, prev: { name: null }, next: { name: "Long pipe" } }).undo()).toBe(1);
    expect(renamed.patchSetup).toHaveBeenCalledWith(7, { name: null });
  });

  it("leaves a setup someone else changed since", async () => {
    const moved = opsWith({ display_name: "Other", name: "Other", included: true });
    const entry = setupEntry(moved.ops, { ...where, prev: { name: null, included: true }, next: { name: "Long pipe", included: false } });
    expect(await entry.undo()).toBe(0);
    expect(moved.patchSetup).not.toHaveBeenCalled();
  });
});
