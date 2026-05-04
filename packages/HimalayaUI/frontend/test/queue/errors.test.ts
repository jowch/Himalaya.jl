import { describe, it, expect } from "vitest";
import { isValidationError, isInfrastructureError, buildValidationMessage } from "../../src/lib/queue/errors";

describe("error classification", () => {
  it("isValidationError matches a 4xx Response error", () => {
    const err = Object.assign(new Error("400 bad request"), { status: 400 });
    expect(isValidationError(err)).toBe(true);
  });

  it("isValidationError rejects non-error inputs", () => {
    expect(isValidationError(null)).toBe(false);
    expect(isValidationError(undefined)).toBe(false);
    expect(isValidationError("string")).toBe(false);
  });

  it("isValidationError rejects 5xx errors", () => {
    const err = Object.assign(new Error("500"), { status: 500 });
    expect(isValidationError(err)).toBe(false);
  });

  it("isInfrastructureError matches 5xx errors", () => {
    const err = Object.assign(new Error("503"), { status: 503 });
    expect(isInfrastructureError(err)).toBe(true);
  });

  it("isInfrastructureError matches network errors (no status)", () => {
    const err = new TypeError("Failed to fetch");
    expect(isInfrastructureError(err)).toBe(true);
  });

  it("isInfrastructureError rejects 4xx", () => {
    const err = Object.assign(new Error("422"), { status: 422 });
    expect(isInfrastructureError(err)).toBe(false);
  });

  it("buildValidationMessage produces a per-kind copy", () => {
    expect(buildValidationMessage("peak_added", new Error("bad q"))).toMatch(/peak/i);
    expect(buildValidationMessage("add_tag", new Error("missing key"))).toMatch(/tag/i);
    // Unknown kind falls back to a generic message.
    expect(buildValidationMessage("update_sample" as any, new Error("oops"))).toBeTruthy();
  });
});
