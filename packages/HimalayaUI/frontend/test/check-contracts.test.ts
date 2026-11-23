import { describe, it, expect } from "vitest";
import { parseContract, checkContract, DECISIONS } from "../scripts/check-contracts.mjs";

const HEAD = [
  "- Mock-up version: `1-a`",
  "- Mock-up sha256: `" + "a".repeat(64) + "`",
  "",
  "| id | Behaviour | Status | Evidence |",
  "|---|---|---|---|",
];

const files: Record<string, string> = {
  "test/x.test.ts": 'it("does the thing", () => {});',
  "test/skip.test.ts": 'it.skip("is skipped", () => {});\nit.todo("is todo");\nit.only("is only", () => {});\ntest.fixme("is fixme", () => {});',
  "test/long.test.ts": 'it("does the thing and more", () => {});',
  "e2e/x.spec.ts": "test('the user sees it', async () => {});\ntest.describe(\"group\", () => {});",
};
const read = (p: string): string | undefined => files[p];

function problems(rows: string[], head = HEAD): string[] {
  return checkContract(parseContract([...head, ...rows].join("\n")), read);
}

describe("parseContract", () => {
  it("reads the pinned version, hash and every table row", () => {
    const c = parseContract([...HEAD, "| a-b | Does it | built | test/x.test.ts › does the thing |"].join("\n"));
    expect(c.version).toBe("1-a");
    expect(c.sha256).toBe("a".repeat(64));
    expect(c.rows).toEqual([{ id: "a-b", behaviour: "Does it", status: "built", evidence: "test/x.test.ts › does the thing", line: 6 }]);
  });

  it("skips the statuses table, whose first column is not a row id", () => {
    const c = parseContract(["| Status | Meaning | Evidence the lint requires |", "|---|---|---|", "| `built` | x | y |"].join("\n"));
    expect(c.rows).toEqual([]);
  });
});

describe("checkContract", () => {
  it("accepts a built row whose test file holds the named test", () => {
    expect(problems(["| a-b | Does it | built | test/x.test.ts › does the thing |"])).toEqual([]);
  });

  it("rejects a built row whose test is gone or renamed", () => {
    expect(problems(["| a-b | Does it | built | test/x.test.ts › does another thing |"])[0]).toMatch(/a-b.*does not contain the test "does another thing"/);
    expect(problems(["| a-b | Does it | built | test/gone.test.ts › does the thing |"])[0]).toMatch(/a-b.*test\/gone\.test\.ts does not exist/);
    expect(problems(["| a-b | Does it | built | some prose |"])[0]).toMatch(/a-b.*path › test title/);
  });

  it("matches the whole test title, not a longer title that contains it", () => {
    expect(problems(["| a-b | Does it | built | test/long.test.ts › does the thing |"])[0]).toMatch(/a-b.*does not contain the test "does the thing"/);
    expect(problems(["| a-b | Does it | built | e2e/x.spec.ts › the user sees it |"])).toEqual([]);
  });

  it("rejects a test that is skipped, todo, focused or marked fixme", () => {
    for (const [title, mod] of [["is skipped", "skip"], ["is todo", "todo"], ["is only", "only"], ["is fixme", "fixme"]]) {
      expect(problems([`| a-b | Does it | built | test/skip.test.ts › ${title} |`])[0]).toMatch(new RegExp(`a-b.*"${title}" is marked \\.${mod}`));
    }
  });

  it("needs a decision for superseded and an owner for planned", () => {
    expect(problems(["| a-b | Does it | superseded | we changed our minds |"])[0]).toMatch(/a-b.*date, a PR number or a doc path/);
    expect(problems(["| a-b | Does it | superseded | Jonathan 2026-11-21 |"])).toEqual([]);
    expect(problems(["| a-b | Does it | planned | someday |"])[0]).toMatch(/a-b.*PR number, an issue or a thread/);
    expect(problems(["| a-b | Does it | planned | PR 12 |"])).toEqual([]);
    expect(problems(["| a-b | Does it | planned | #34 |"])).toEqual([]);
    expect(problems(['| a-b | Does it | planned | thread "Do it" |'])).toEqual([]);
  });

  it("rejects an unknown status, empty evidence and a repeated id", () => {
    expect(problems(["| a-b | Does it | done | x |"])[0]).toMatch(/a-b.*status "done"/);
    expect(problems(["| a-b | Does it | untested |  |"])[0]).toMatch(/a-b.*no evidence/);
    expect(problems(["| a-b | Does it | untested | x |", "| a-b | Again | untested | y |"])[0]).toMatch(/a-b.*more than once/);
  });

  it("needs the pinned mock-up version and hash", () => {
    expect(problems(["| a-b | Does it | untested | x |"], HEAD.slice(2))).toEqual([
      "docs/mockup-contract.md does not pin the mock-up version",
      "docs/mockup-contract.md does not pin the mock-up sha256",
    ]);
  });
});

describe("checkContract on the decisions file", () => {
  const DHEAD = ["| id | Decision | Decided | Status | Evidence |", "|---|---|---|---|---|"];
  const jl = { "packages/x/test/test_a.jl": '@testset "a decision holds" begin end' };
  const dread = (p: string): string | undefined => ({ ...jl, ...files } as Record<string, string>)[p];
  const dproblems = (rows: string[]): string[] =>
    checkContract(parseContract([...DHEAD, ...rows].join("\n")), dread, DECISIONS);

  it("reads the Status and Evidence columns by name, after the Decided column", () => {
    const c = parseContract([...DHEAD, "| a-b | Do it | 2026-11-23 | planned | PR 9 |"].join("\n"));
    expect(c.rows).toEqual([{ id: "a-b", behaviour: "Do it", decided: "2026-11-23", status: "planned", evidence: "PR 9", line: 3 }]);
  });

  it("accepts a Julia @testset as a built row's test, and needs no mock-up pin", () => {
    expect(dproblems(["| a-b | Do it | 2026-11-23 | built | packages/x/test/test_a.jl › a decision holds |"])).toEqual([]);
    expect(dproblems(["| a-b | Do it | 2026-11-23 | built | packages/x/test/test_a.jl › a decision", ].map((r) => r + " |"))[0]).toMatch(/does not contain the test "a decision"/);
  });

  it("needs the date each decision was made", () => {
    expect(dproblems(["| a-b | Do it | soon | planned | PR 9 |"])[0]).toMatch(/a-b.*Decided column must be a YYYY-MM-DD date/);
  });
});

