# The Print finish, round 3 — Agent-Team orchestration design

**Date:** 2026-05-28
**Milestone:** [#4 — HimalayaUI — The Print finish, round 3](https://github.com/jowch/Himalaya.jl/milestone/4) (8 issues) **+ folded-in** [#252](https://github.com/jowch/Himalaya.jl/issues/252) (figure-export palette invariant; no milestone)
**Status:** wave 1 COMPLETE + merged (2026-05-28, PRs #262/#263/#264); the wave-1 Monitor experiment is **resolved (REFUTED)** and waves 2–4 are confirmed on the orchestrator-nudge model — see [Wave-1 results & confirmed decisions](#wave-1-results--confirmed-decisions-2026-05-28) at the end of this doc.
**Base:** `main` = `origin/main` @ `ef6071f` (incl. #253 round-3 findings doc). Reconciled before wave 1; all implementer worktrees branch off this.

## Purpose

Parallelize the round-3 fidelity remediation across a **persistent Agent Team** — one member per issue (or per file-cluster) — with a paired author/reviewer review loop per PR.

It refines the literal "one member per issue, two independent teams" framing into a shape that respects (a) the real file-level dependency DAG, (b) the human-interactive nature of `brainstorming`, (c) the persistent-teammate wake/Monitor mechanics (unproven — wave 1 is the experiment), and (d) the fact that **visual** acceptance criteria need a rendered surface against real data, not just a green build.

This doc is the orchestrator's + operator's design record. The orchestrator distributes the runbook + per-issue kit to each member **in its spawn prompt** (members branch off `main`, which does not carry this doc).

## Source-of-truth bundle (per issue)

Every implementer's kit (all present on `main` @ `ef6071f`):
- The issue body itself (Context + Scope + Done-when + exact `file:line` — treated as the **approved spec**).
- [`DESIGN.md`](../../../DESIGN.md) — canonical "The Print" token/type/component spec.
- [`docs/2026-05-28-the-print-round3-findings.md`](../../2026-05-28-the-print-round3-findings.md) — round-3 audit; defines the R3-* finding IDs + U-1…U-5 operator adds the issues cite.
- [`docs/2026-05-27-the-print-round2-findings.md`](../../2026-05-27-the-print-round2-findings.md) + [`docs/2026-05-27-the-print-fidelity-findings.md`](../../2026-05-27-the-print-fidelity-findings.md) — prior rounds (re-judgement context).
- `docs/redesign-mockups/*.html` — pixel targets.

## Why no per-issue brainstorming

`brainstorming`'s HARD-GATE requires a **human-approved, spec-grade design before any implementation skill fires**. That deliverable **already exists upstream**: milestone #4's issues were authored from the round-3 findings doc + DESIGN.md + mockups (the same keystone-DAG design process milestone-3 used for R0a–R10), and each issue body carries Context + Scope + Done-when + `file:line`. We bypass the per-issue brainstorming *dialogue*, **not** the *gate* — the approved design artifact is the milestone + findings + DESIGN.md. (This is the legitimate justification; "one human can't run 8 brainstorms" is a throughput observation, not the reason.) The B+G+H multi-file collapse — the one decomposition call brainstorming would normally own — is made explicitly here (shared `IMAGE_PROCESSING_VERSION`, see DAG).

## Topology — one role-tagged persistent team

A **single** team (`print-r3`), not two. Reviewer≠author independence comes from distinct **agents**, not distinct **teams** — this avoids cross-team `SendMessage` routing risk and the dead-`agentId` name-reuse bug.

- **Orchestrator** = root session. Owns the shared task list (= the DAG via `addBlockedBy` edges). Runs the spec/plan gate. Bridges teammate↔human approval. Holds the orchestrator-side Monitor + nudge fallback (see experiment).
- **Implementer members** (`impl-252`, `impl-254`, …) — `general-purpose` agent type (need Edit/Write/Bash). One per issue or per file-cluster. **Distinct names, never reused.**
- **Reviewer members** (`rev-<PR#>`) — drive `review-pr`; invoke project subagent reviewers (`frontend-reviewer`, `himalaya-reviewer`, `saxs-physics-reviewer`, `queue-reviewer`) as their analysis engine. **Named per-PR, not reused across waves** (same name-routing reason).

## The dependency DAG (waves)

The 8+1 issues are **not** independent lanes. They cluster by shared file. Each wave is internally **file-disjoint** so parallel members never edit the same file; later waves branch fresh off `main` after earlier waves merge, so cross-wave overlaps resolve by rebase rather than by hand.

| Wave | Members → issues | Files | Ordering rationale |
|---|---|---|---|
| **1** | `impl-252` → #252 · `impl-254` → A/#254 · `impl-258` → E/#258 | #252: `lib/figure-export/presets.ts`, `lib/comparison/coloring.ts`, tests · A: `TraceViewer.tsx`, `PlotCard.tsx` · E: `GroupingModeToggle.tsx`, `SeriesBuilderRail.tsx`, `MemberHeatmapLayer.tsx`, `Series*Page.tsx`, `CrossTraceTrackingLayer.tsx` | Mutually file-disjoint. **Probe wave** for the wake/Monitor experiment. |
| **2** | `impl-256` → C/#256 · `impl-257` → D/#257 | C: `ContactSheetRow.tsx`, `LoupeSidebar.tsx`, `SamplesPage.tsx`, `DetectorImage.tsx` · D: `PlotCard.tsx`, `TraceViewer.tsx`, `PhasePanel.tsx`, `FocusNotesMargin.tsx`, `FocusReflectionsTable.tsx` | C ⊥ D. D inherits A's merged `PlotCard`/`TraceViewer` edits → no conflict. **`styles.css` caveat below.** |
| **3** | `impl-bgh` → B+G+H (#255+#260+#261), **one PR** | `image.jl`, `routes_exposures.jl`, `cli.jl`, `pipeline.jl`, `DetectorImage.tsx`, `bones/` | All three touch `image.jl`/`routes_exposures.jl` and **all bump `IMAGE_PROCESSING_VERSION`** → must be one body of work (the explicit one-responsibility exception). B's `DetectorImage` half inherits C's merge. |
| **4** | `impl-259` → F/#259, **last** | project-wide `frontend/src/` + new CI guard | Global token-name sweep + CI grep guard. Runs LAST so the sweep covers final state and the guard can't false-fail in-flight branches. |

**`styles.css` is a shared global file.** Wave-2's C (R3-S04, `--text-progress-numeral`) and D (R3-F09, `text-th-caps`) both *optionally* add scale-token roles to it. Both issues offer a "reuse an existing token" path (C → `text-headline-lg`; D → `text-meta`). **Members must take the reuse path** to keep wave 2 disjoint; if a new role is truly needed, that member owns the `styles.css` edit and flags it so the other rebases. This is a *soft* guardrail (member judgment, not structural disjointness), so the **orchestrator confirms at the plan gate (step 3) that both wave-2 plans took the reuse path** — caught there, not discovered at merge.

**Cluster cross-reference (why the grouping):**
- `image.jl` + `IMAGE_PROCESSING_VERSION`: **B, G, H** → one PR.
- `DetectorImage.tsx`: **B, C** → C (wave 2) before B (wave 3).
- `PlotCard.tsx`/`TraceViewer.tsx`: **A, D** → A (wave 1) before D (wave 2).
- Global class-name sweep + CI guard: **F** → strictly last; D/#257 R3-F03 also clears ~half of F's PlotCard token surface, so F rebases onto a smaller residue.

## Per-implementer runbook (plan-only)

The issue body **is** the spec. No re-brainstorm (see above).

1. **Workspace** — `using-git-worktrees` to create the worktree + branch (`r3-<letter>-<slug>`) off current `main`; then the project `worktree-setup` skill for `npm install` + `Pkg.instantiate` inside it.
2. **Plan** — `writing-plans`: issue body → implementation plan + bite-sized TDD task list, committed to the branch.
3. **Dual plan gate** — orchestrator runs `verify-before-review` (an explicit pass that **greps the plan's cited symbols / paths / line-numbers against `main`** — not a skim) + one project reviewer on the plan, then surfaces it to the **human** for approval. Member idles until the orchestrator sends "approved" (the gate is a discrete message, not an in-context pause — see experiment notes).
4. **Implement** — execute the plan via `executing-plans` (the outer loop: walk the task list, stop-when-blocked, review checkpoints) with `test-driven-development` as the inner loop per task (red → green → refactor, each step its own commit). For **visual-only** findings (token swap, accent color, checkbox state) the TDD artifact is a **`data-*`-attribute assertion** (per `frontend/test/AGENTS.md` — *never* a Tailwind class-string test) **plus** the operator visual check in the harness below; do not skip the test, and do not assert on class strings. Use **regression floors/ceilings, not hard counts**, for fixture-based assertions. (`executing-plans` would normally hand off to `finishing-a-development-branch` at task-list completion; here that terminal handoff is **replaced by step 6's `request-pr-review`** — do not invoke `finishing-a-development-branch`.)
5. **Verify** — `verification-before-completion`. Frontend gate **must** include `npm run build` (Vitest alone masks Vite prod-build breakage — milestone-3 lesson), plus Vitest. Backend (wave 3) runs the Julia suite **once, capture-to-file-and-grep** (it's 5–10 min and rebuilds fixtures every run — never re-run per grep). Visual Done-whens verified in the harness below.
6. **PR + review loop** — open PR → `request-pr-review` (author half; intentionally stands in for the `finishing-a-development-branch` handoff). Paired with a `rev-<PR#>` member running `review-pr` until convergence → merge → member self-resumes and claims its next unblocked task.

## Visual-acceptance harness (required for visual Done-whens)

Several round-3 "Done when"s are **visual** and unverifiable by build/test alone: #254 "≤2 terracotta marks per focus-indexed view," #255 "warm-toned detector image diff," #258 "terracotta checkbox check-state," and others. Milestone-3 hit this as an env gap; **this run the data is present**: `~/projects/himalaya-devdata/himalaya.db` exists and `/Volumes/data` (detector TIFFs) is mounted.

Harness (orchestrator stands up once, members/reviewers use):
- Backend `serve` against a **writable copy** of `~/projects/himalaya-devdata/himalaya.db` (never the original), `/Volumes/data` mounted for detector images. (No sysimage built — cold start is slower; build `make sysimage` once if serve cost dominates.)
- Each visual member runs its branch's Vite on a **unique port** pointed at the backend; captures the relevant surface via Playwright (mocked specs for structure, live backend for real-data visuals).
- The visual Done-when is checked against the rendered surface (e.g. count terracotta marks; compare detector warmth). A member may **not** mark a visual issue done on a green build alone. This visual check is an acceptance gate the **orchestrator owns and signs off**: `verification-before-completion` certifies only the build/test half, never the visual judgment — so "verification-before-completion passed" is necessary but not sufficient for a visual issue.
- Tag/curation-seeded surfaces (if any) need a private seeded backend on a separate port (milestone-3 used :8092 for R7's warm scoping path).

## Review insertion points

1. **Project subagent reviewers** (`frontend-reviewer`, `himalaya-reviewer`, `saxs-physics-reviewer`, `queue-reviewer`) — stateless one-shot diff reviewers. Used (a) by the orchestrator at the **plan gate** via `verify-before-review`, and (b) by a reviewer member as the **analysis engine** inside `review-pr` (which itself reads `AGENTS.md`/`CLAUDE.md` and folds them in).
2. **`request-pr-review` ↔ `review-pr`** — the multi-round converging loop (author/reviewer halves; they converge through the PR, never talk directly). **Distinct** from the superpowers `requesting-code-review`/`receiving-code-review` one-shot pair.

Issue → project reviewer: frontend (#252, A, C, D, E, F) → `frontend-reviewer`; backend image (B, G, H) → `himalaya-reviewer` **+ `saxs-physics-reviewer`** (B's LUT touches detector-image rendering).

## Wave-1 experiment & decision gate

Open question (`team-monitor-experiment` memory): *do persistent team members self-resume on their own armed Monitors without orchestrator nudging?* A research-based review argues **no** (idle teammates wake on messages/task-claims, not Monitor events); the operator's hands-on read is **yes** for *persistent* members (vs dead generic subagents). **Unresolved → test it cleanly.**

**Clean single-subject probe (per the memory's protocol):**
- `impl-252` is the **sole instrumented subject**: it arms its own Monitor via `request-pr-review`; the orchestrator runs **no** parallel watcher for its PR and does **not** nudge it. We observe whether it self-resumes when its `rev-252` reviewer posts.
- `impl-254` and `impl-258` (and all wave-2+ members) run on the **orchestrator-Monitor + `SendMessage` nudge** fallback **from the start** — so the wave stays productive regardless of the probe outcome, and nudging them doesn't contaminate the probe.

**Decision:**
- `impl-252` self-resumes cleanly → teammate-armed Monitors work; waves 2–4 drop the orchestrator watcher (orchestrator collapses to spawn + DAG-gating + plan/human-approval bridging + status collection). Update the memory: **confirmed**.
- It stalls → orchestrator-Monitor + nudge stays load-bearing for all PRs (milestone-3 model). Update the memory: **refuted**, and note `request-pr-review`/`review-pr` are root-session-only for now.

## Known risks & mitigations

| Risk | Mitigation |
|---|---|
| Teammate Monitor doesn't self-resume | The wave-1 single-subject experiment; orchestrator-Monitor + nudge fallback active for all non-probe PRs from the start. |
| Bare-name `SendMessage` routes to a dead `agentId` | Distinct member names, never reused (implementers AND per-PR reviewers); `@team` form if needed. |
| Visual Done-whens unverifiable by build/test | Live-screenshot harness (dev DB + `/Volumes/data`, per-branch Vite + Playwright) is a required acceptance step; visual-only TDD = `data-*` assertion + operator check. |
| Cross-wave file overlap (A↔D, B↔C) | Wave ordering + branch-off-`main`-after-merge → rebase, not hand-merge. |
| Within-wave `styles.css` overlap (C↔D) | Members take the reuse-existing-token path; else one member owns the additive edit and flags it. |
| `IMAGE_PROCESSING_VERSION` triple-bump | B+G+H collapsed into one member / one PR. |
| F's CI guard false-failing in-flight branches | F scheduled strictly last. |
| Stale `node_modules` on `main` after worktree installs (milestone-3 trap) | Pinned gate: **after each wave merges, re-run `npm install` + `npm run build` on the `main` checkout before opening the next wave.** |
| Slow-suite thrash across ~7 PRs | Capture-once + grep for Julia; `npm run build` required (not just Vitest); wave-3 member expects the 5–10 min Julia suite per review round. |
| Zombie/uncleaned idle teammates accumulating | Orchestrator tracks member liveness; shut down completed members via `shutdown_request`; don't reuse names. |
| One human can't run 8 interactive brainstorms | Plan-only: design exists upstream (milestone + findings + DESIGN.md). |

## Naming conventions

- Team: `print-r3`. Implementers: `impl-<issue#>` (or `impl-bgh`). Reviewers: `rev-<PR#>`.
- One worktree + one branch + one PR per member (except B+G+H = one of each).
- Branch names: `r3-<letter>-<slug>` (e.g. `r3-a-accent-rationing`).


## Wave-1 results & confirmed decisions (2026-05-28)

Wave 1 (#252 palette · #254 accent · #258 series) shipped end-to-end: 3 issues → 3 PRs (#262/#263/#264), all **round-1 approved**, merged to `main` @ `6d73d64`, build-green. The full pipeline ran: spawn → plan → dual gate → TDD → PR → review loop → visual acceptance → merge.

### The Monitor experiment — REFUTED

**Persistent team members do NOT self-resume on their own armed Monitors across idle.** `impl-252` armed a Monitor on PR #262 and went idle; `rev-262` posted its review; impl-252 stayed silent ~9.5 min. Confirmed: idle teammates wake on `SendMessage` / task-claims, not on their own Monitor stdout. (See the `team-monitor-experiment` memory.)

Sharper finding: **the nudge bridge is needed on BOTH sides** — the orchestrator had to nudge `rev-262` for round 2 as well as the author. So `request-pr-review` / `review-pr` is effectively **root-session-only inside a team**; the orchestrator drives every round transition and **cannot collapse to "spawn + collect."** Waves 2–4 run the orchestrator-watch + nudge model from the start (no per-wave probe).

### Confirmed decisions

- **Keep the full `request-pr-review` / `review-pr` loop on every PR, even when the review is trivial.** Rationale (operator): the loop's GitHub review thread is a first-class **audit-trail deliverable**, not just a quality gate. Do not downgrade low-risk PRs to gate-only.
- **Wave width stays 2 (no aggressive fan-out).** The 2-wide cadence keeps the orchestrator's per-round nudge load manageable; since the orchestrator is the serialization point for review rounds, width is bounded by bridge bandwidth, not by file-disjointness alone.
- **Seed the visual harness with one indexed series.** The dev DB has no indexed/ordered *series*, so series-heatmap Done-whens (e.g. #258's keyline/axis) fell back to `data-*` unit tests. Seeding one indexed series makes those pixel-verifiable for the remaining waves.

### Operational note — team-member worktree pinning

Spawned team members' `Write`/`Edit` tools (and the orchestrator's, once its session isolates) **pin to a single worktree** (`r3-e-series-finish` this run), not each member's own. Workaround (in use): author all edits via **Bash heredoc with absolute paths** in the correct worktree (`Bash` is not pinned). Consequences: (1) do **not** garbage-collect the pinned worktree mid-milestone — live members depend on it; (2) TDD edits are more fragile (no surgical `Edit`), so members must be careful. Not fixed this run (operator: stay with the workaround).

### What the gate caught (value evidence)

The plan gate (orchestrator verify-before-review grep + project reviewer + human) caught, **before any code**: #252's `[7]→285` swatch-collision (a bug in the *issue text itself*), #258's heatmap-keyline test-coupling (green unreachable as written) plus two non-existent test-helper names, and #256's "zero inline `text-[Npx]`" Done-when overshooting the issue's actual scope. Grepping cited values/paths against live source (not skimming) is what surfaced these.
