import { describe, it, expect } from "vitest";
import { parseContract, checkContract, checkRowCount, decisionRefs, deletedRows, DECISIONS, MOCKUP } from "../scripts/check-contracts.mjs";

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
  "test/commented.test.ts": '// it("is commented out", () => {});\n/*\n  it("is in a block comment", () => {});\n*/\nconst glob = "src/**/*.ts";\n  it("runs after a glob string", () => {});',
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

  it("rejects a test that is commented out, by a line or a block comment", () => {
    expect(problems(["| a-b | Does it | built | test/commented.test.ts › is commented out |"])[0]).toMatch(/does not contain the test "is commented out"/);
    expect(problems(["| a-b | Does it | built | test/commented.test.ts › is in a block comment |"])[0]).toMatch(/does not contain the test "is in a block comment"/);
  });

  it("needs a decision for superseded and an owner for planned", () => {
    expect(problems(["| a-b | Does it | superseded | we changed our minds |"])[0]).toMatch(/a-b.*date, a PR number or a doc path/);
    expect(problems(["| a-b | Does it | superseded | Jonathan 2026-11-21 |"])).toEqual([]);
    expect(problems(["| a-b | Does it | planned | someday |"])[0]).toMatch(/a-b.*a planned row must name its owner: a PR number, an issue or a thread/);
    expect(problems(["| a-b | Does it | planned | PR 12 |"])).toEqual([]);
    expect(problems(["| a-b | Does it | planned | #34 |"])).toEqual([]);
    expect(problems(['| a-b | Does it | planned | thread "Do it" |'])).toEqual([]);
    expect(problems(["| a-b | Does it | contradicted | the list shows it |"])[0]).toMatch(/a-b.*a contradicted row must name its owner/);
    expect(problems(["| a-b | Does it | contradicted | #7: the list shows it |"])).toEqual([]);
  });

  it("rejects an unknown status, empty evidence and a repeated id", () => {
    expect(problems(["| a-b | Does it | done | x |"])[0]).toMatch(/a-b.*status "done"/);
    expect(problems(["| a-b | Does it | untested |  |"])[0]).toMatch(/a-b.*no evidence/);
    expect(problems(["| a-b | Does it | untested | x |", "| a-b | Again | untested | y |"])[0]).toMatch(/a-b.*more than once/);
  });

  it("rejects a superseded row that cites a decision id docs/decisions.md doesn't have", () => {
    const ids = new Set(["new-rule"]);
    const check = (row: string): string[] => checkContract(parseContract([...HEAD, row].join("\n")), read, MOCKUP, ids);
    expect(check("| a-b | Does it | superseded | new-rule in docs/decisions.md, 2026-11-23 |")).toEqual([]);
    expect(check("| a-b | Does it | superseded | gone-rule in docs/decisions.md, 2026-11-23 |")[0]).toMatch(/a-b.*cites gone-rule in docs\/decisions\.md, which has no row with that id/);
    expect(check("| a-b | Does it | superseded | docs/decisions.md new-rule, 2026-11-23 |")).toEqual([]);
    expect(check("| a-b | Does it | superseded | 2026-11-23 docs/decisions.md gone-rule |")[0]).toMatch(/a-b.*cites gone-rule in docs\/decisions\.md, which has no row with that id/);
    expect(check("| a-b | Does it | superseded | see docs/decisions.md, 2026-11-23 |")[0]).toMatch(/a-b.*mentions docs\/decisions\.md without naming a decision/);
  });

  it("reads a decision cited either way round", () => {
    expect(decisionRefs("combine-select in docs/decisions.md, 2026-11-23")).toEqual(["combine-select"]);
    expect(decisionRefs("2026-11-23 docs/decisions.md phase-curvature-legend")).toEqual(["phase-curvature-legend"]);
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
  const jl = {
    "packages/x/test/test_a.jl": '@testset "a decision holds" begin end\n# @testset "was disabled" begin end\n#=\n@testset "sits in a block comment" begin end\n=#',
  };
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

  it("rejects a Julia testset that is commented out, by a line or a block comment", () => {
    expect(dproblems(["| a-b | Do it | 2026-11-23 | built | packages/x/test/test_a.jl › was disabled |"])[0]).toMatch(/does not contain the test "was disabled"/);
    expect(dproblems(["| a-b | Do it | 2026-11-23 | built | packages/x/test/test_a.jl › sits in a block comment |"])[0]).toMatch(/does not contain the test "sits in a block comment"/);
  });

  it("needs the date each decision was made", () => {
    expect(dproblems(["| a-b | Do it | soon | planned | PR 9 |"])[0]).toMatch(/a-b.*Decided column must be a YYYY-MM-DD date/);
  });

  it("needs a superseded decision to start with the id of an existing decision that replaced it", () => {
    const rows = (evidence: string): string[] => [
      `| old-rule | Do it | 2026-11-20 | superseded | ${evidence} |`,
      "| new-rule | Do it better | 2026-11-23 | planned | PR 9 |",
    ];
    expect(dproblems(rows("new-rule, 2026-11-23"))).toEqual([]);
    expect(dproblems(rows("gone-rule, 2026-11-23"))[0]).toMatch(/old-rule.*superseded by gone-rule, which has no row in docs\/decisions\.md/);
    expect(dproblems(rows("old-rule, 2026-11-23"))[0]).toMatch(/old-rule.*can't supersede itself/);
    expect(dproblems(rows("Jonathan changed it on 2026-11-23"))[0]).toMatch(/old-rule.*starts its evidence with the id of the decision that replaced it, optionally after its date/);
    expect(dproblems(rows("2026-11-23: new-rule and another"))).toEqual([]);
    expect(dproblems(rows("2026-11-23: gone-rule"))[0]).toMatch(/old-rule.*superseded by gone-rule/);
  });
});

describe("checkRowCount", () => {
  const doc = { ...DECISIONS, floor: 2 };
  const rows = (n: number) => parseContract(["| id | Decision | Decided | Status | Evidence |", "|---|---|---|---|---|",
    ...Array.from({ length: n }, (_, i) => `| r-${i} | Do it | 2026-11-23 | planned | PR 9 |`)].join("\n"));

  it("fails on a file a merge emptied", () => {
    expect(checkRowCount(parseContract(""), doc)[0]).toMatch(/docs\/decisions\.md has no rows/);
  });

  it("fails when rows drop below the floor, and passes at it", () => {
    expect(checkRowCount(rows(1), doc)[0]).toMatch(/has 1 rows, fewer than its floor of 2/);
    expect(checkRowCount(rows(2), doc)).toEqual([]);
  });

  it("names every row the base commit had that is gone now", () => {
    expect(deletedRows(rows(3), rows(2))).toEqual(["r-2"]);
    expect(deletedRows(rows(2), rows(3))).toEqual([]);
  });

  it("holds both real files to a floor", () => {
    expect(MOCKUP.floor).toBeGreaterThan(0);
    expect(DECISIONS.floor).toBeGreaterThan(0);
  });
});
