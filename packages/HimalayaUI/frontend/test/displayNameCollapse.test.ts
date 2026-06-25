import { describe, it, expect } from "vitest";
import { sampleDisplayName } from "../src/lib/sample/displayName";
import { SCHEMA_VERSION } from "../src/lib/queue/persistence";

describe("display_name → name collapse (Phase E1, Task 1b)", () => {
  it("sampleDisplayName reads name only", () => {
    expect(sampleDisplayName({ id: 7, name: "HA85 (S01P15)" })).toBe("HA85 (S01P15)");
    expect(sampleDisplayName({ id: 7, name: "" })).toBe("Sample #7");
  });
  it("queue-op SCHEMA_VERSION is bumped to 5 for the label collapse", () => {
    expect(SCHEMA_VERSION).toBe(5);
  });
});
