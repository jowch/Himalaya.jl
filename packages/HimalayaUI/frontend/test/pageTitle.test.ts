import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const HTML = readFileSync(resolve(__dirname, "../index.html"), "utf-8");

describe("page title", () => {
  it("the browser tab is titled Himalaya", () => {
    const titles = [...HTML.matchAll(/<title>([^<]*)<\/title>/g)].map((m) => m[1]);
    expect(titles).toEqual(["Himalaya"]);
  });
});
