import { describe, it, expect } from "vitest";
import { Linter } from "eslint";
import tseslint from "typescript-eslint";
import uiCopy from "../scripts/eslint-ui-copy.mjs";

const linter = new Linter({ configType: "flat" });

function messages(code: string): string[] {
  return linter
    .verify(code, [
      {
        files: ["**/*.tsx"],
        languageOptions: { parser: tseslint.parser, parserOptions: { ecmaFeatures: { jsx: true } } },
        plugins: { copy: uiCopy },
        rules: { "copy/ui-copy": "error" },
      },
    ], "Probe.tsx")
    .map((m) => m.message);
}

describe("ui-copy lint rule", () => {
  it("flags a telegraphic click hint in a title", () => {
    expect(messages('const a = <span title="Click hides. ⌘-click subtracts." />;')[0]).toMatch(/Click to hide/);
  });

  it("flags a telegraphic key hint in a palette hint field", () => {
    expect(messages('const a = { hint: "On this sample, Enter removes it" };')[0]).toMatch(/Press Enter to/);
  });

  it("flags telegraphic hints inside template literals and ternaries", () => {
    expect(messages("const a = <span title={`${x}-click subtracts.`} />;")).toHaveLength(1);
    expect(messages('const a = <span title={on ? "Click hides." : "Click to show."} />;')).toHaveLength(1);
  });

  it("accepts imperative hints and label tooltips", () => {
    expect(messages('const a = <span title="Click to hide. ⌘-click to subtract." />;')).toEqual([]);
    expect(messages('const a = <Tooltip content="Keyboard shortcuts" />;')).toEqual([]);
    expect(messages("const a = <span title={`Click to ${off ? \"show\" : \"hide\"}.`} />;")).toEqual([]);
  });

  it("requires a string that holds a sentence to end with a period", () => {
    expect(messages('const a = <button title="Suggested. Click to confirm" />;')[0]).toMatch(/period/);
    expect(messages('const a = <span>Press Enter to add</span>;')[0]).toMatch(/period/);
    expect(messages('const a = <button title="Suggested. Click to confirm." />;')).toEqual([]);
  });

  it("enforces the house vocabulary in user-facing copy only", () => {
    expect(messages("const a = <p>Drop this exposure.</p>;")[0]).toMatch(/image/);
    expect(messages('const a = { what: "Previous or next frame" };')[0]).toMatch(/image/);
    expect(messages('const a = <p>Add a tag.</p>;')[0]).toMatch(/label/);
    expect(messages('const a = <p>Mute this image.</p>;')[0]).toMatch(/Drop/);
    expect(messages('const a = <dt>Exposure time</dt>;')).toEqual([]);
    expect(messages('const a = queryKey("exposure", "frame");')).toEqual([]);
    expect(messages('const a = <div data-testid="frame-strip" className="frame" />;')).toEqual([]);
  });

  it("flags title case", () => {
    expect(messages('const a = <h2>Keyboard Shortcut Help</h2>;')[0]).toMatch(/sentence case/);
    expect(messages('const a = <h2>Keyboard shortcuts</h2>;')).toEqual([]);
  });
});
