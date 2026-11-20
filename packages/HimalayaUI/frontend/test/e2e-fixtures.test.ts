import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";

const E2E = join(__dirname, "..", "e2e");

describe("mocked e2e specs", () => {
  const specs = (readdirSync(E2E, { recursive: true }) as string[])
    .filter((f) => f.endsWith(".spec.ts") && !f.startsWith("live/"));

  it("exist", () => {
    expect(specs.length).toBeGreaterThan(0);
  });

  it.each(specs)("%s takes test from e2e/fixtures.ts, which fails on any /api request the spec does not mock", (file) => {
    const src = readFileSync(join(E2E, file), "utf8");
    expect(src).not.toMatch(/from\s+["']@playwright\/test["']/);
    expect(src).toMatch(/import\s*\{[^}]*\btest\b[^}]*\}\s*from\s*["']\.\/fixtures["']/);
  });
});
