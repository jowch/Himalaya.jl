import { describe, expect, it } from "vitest";
import { ESLint } from "eslint";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";
import { Button } from "../../src/ui/Button";
import { Card } from "../../src/ui/Card";
import { Checkbox } from "../../src/ui/Checkbox";
import { Dot } from "../../src/ui/Dot";
import { Input } from "../../src/ui/Input";
import { Kicker } from "../../src/ui/Kicker";
import { PhaseChip } from "../../src/ui/PhaseChip";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "../..");

async function titleErrors(code: string): Promise<number> {
  const eslint = new ESLint({ cwd: ROOT });
  const [result] = await eslint.lintText(code, { filePath: resolve(ROOT, "src/__title_probe.tsx") });
  return result!.messages.filter((m) => m.ruleId === "no-restricted-syntax").length;
}

describe("native title tooltips are banned", () => {
  it("rejects title on DOM elements and on router links, and allows component props named title", async () => {
    expect(await titleErrors(`export const A = () => <span title="x">a</span>;`)).toBe(1);
    expect(await titleErrors(`import { Link } from "react-router";\nexport const A = () => <Link to="/" title="x">a</Link>;`)).toBe(1);
    expect(await titleErrors(`import { NavLink } from "react-router";\nexport const A = () => <NavLink to="/" title="x">a</NavLink>;`)).toBe(1);
    expect(await titleErrors(`import { Section } from "./ui/Section";\nexport const A = () => <Section title="Call">a</Section>;`)).toBe(0);
  });

  it("the primitives that spread props onto the DOM do not accept title", () => {
    const rejected = [
      // @ts-expect-error title would render a native tooltip
      <Button key="b" title="x">a</Button>,
      // @ts-expect-error title would render a native tooltip
      <Card key="c" title="x">a</Card>,
      // @ts-expect-error title would render a native tooltip
      <Checkbox key="k" checked={false} onChange={() => {}} title="x" />,
      // @ts-expect-error title would render a native tooltip
      <Dot key="d" title="x" />,
      // @ts-expect-error title would render a native tooltip
      <Input key="i" title="x" />,
      // @ts-expect-error title would render a native tooltip
      <Kicker key="h" title="x">a</Kicker>,
      // @ts-expect-error title would render a native tooltip
      <PhaseChip key="p" phase="Pn3m" title="x" />,
    ];
    expect(rejected).toHaveLength(7);
  });
});
