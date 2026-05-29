---
name: orchestrating-agent-teams
description: Use when you are the root session driving a GitHub milestone of several already-spec-grade issues, and want long-lived implementer + reviewer teammate agents to carry them to merge across dependency-ordered waves with a GitHub audit trail. Use when coordinating multiple persistent teammates (TeamCreate/Agent/SendMessage) rather than one-shot subagents. Not for a single ad-hoc task or work that still needs design.
---

# orchestrating-agent-teams

## Overview

You are the **orchestrator** (root session). You spawn `impl-<issue#>` implementers and `rev-<PR#>` reviewers and drive them to merged PRs.

**Core principle:** you hold the *bookkeeping layer* (nudging, watching PRs, running long verifies, sweeping for strays) — tedious for a human, cheap for you. **Keep your own judgment for the seams** (root-causing a regression, questioning an issue's premise, catching scope creep). The system works because a tireless executor is wrapped in *independent* checks (human judgment, a plan-gate reviewer, the push guard); surface seams to the human, don't sail past them.

## When to use

- A milestone whose **issue bodies are already spec-grade** (issue + linked findings = the approved design; no per-issue brainstorming).
- You want **parallel implementers + a per-PR review thread as an audit-trail deliverable**.
- You will **stay live** as the root session to drive it.

Not for: a single issue (just spawn one agent), or work needing design (brainstorm first).

## Load-bearing gotchas (read before spawning anyone)

These silently stall or corrupt — they are the whole reason naive "spawn + collect" fails.

1. **Idle teammates do NOT self-resume** — not on their own Monitors, not on their own background-job completion, and their background jobs get **reaped** on idle. They wake only on `SendMessage`/task-claims. So: **you nudge BOTH sides of every review round**, and **you run long verifies yourself** (a teammate that backgrounds a 5-min suite then idles loses the process). A root-session `gh pr list` poll Monitor DOES wake you — use one as a PR-open fallback.
2. **`Write`/`Edit` pin to ONE worktree** (and the root can pin once isolated). Author edits via **Bash heredoc/sed with absolute paths** (`Bash` is unpinned). Don't GC the pinned worktree mid-milestone.
3. **Teammates must run git from the worktree path explicitly**, or commits land on local `main`. Before any authorized `reset --hard main`, verify `git log main --not origin/main <branch>` is empty.
4. **Verify by invariant, not enumeration** — an acceptance grep keyed off a hand-listed set is blind to the leak it guards. Match the structural pattern across all forms.
5. **Discovery ≠ verification** — fan-out workflow to *scope* the unknown; deterministic grep+build+diff to *verify* the known. Don't re-run discovery to recount the known.

## Sequencing: the wave DAG and file collisions

The hard constraint: **two members editing the same file in parallel collide** — at merge the second PR conflicts, or (worse, with Bash-heredoc edits that can't surgically 3-way-merge) silently clobbers the first. So sequencing is about **file-disjointness**, not priority. Build the DAG before spawning anyone:

1. **Map issue → files.** For each issue, list the files its scope touches (grep the cited paths). Include the collision hotspots: shared primitives (`ui/` components), shared CSS/design tokens, and any file that gets a version bump — these are touched by many issues.
2. **Cluster by overlap.** Any two issues sharing a file **cannot be in the same wave.**
3. **Within a wave:** a maximally file-disjoint set. **Across waves:** dependency order — a wave whose surface builds on an earlier change comes later and **branches off `main` AFTER the earlier wave merges** (rebase onto the new `main`; never hand-merge between feature branches).
4. **Collapse hard-shared issues, don't split them.** If two issues *must* edit the same file (a shared component, one version bump), separate waves only defer the collision — put them in **one member / one PR**. (E.g. three issues all bumping a shared processing version ship as a single PR.)
5. **Cap width at ~2 regardless of disjointness.** Even five file-disjoint issues shouldn't run five-wide — you are the serialization point for review rounds, so width is bounded by your nudge bandwidth, not the DAG.
6. **Watch soft collisions on shared globals.** A wave can be file-disjoint yet still collide on a shared *global* file (`styles.css`, a design-token file, a shared config) that two members each *optionally* touch. This isn't structural — mitigate by having members take a reuse-an-existing path, and **confirm at the plan gate that they did** (caught there, not at merge).

```
Wave 1:  #A (files x,y)   #B (file z)            disjoint → parallel
         └─────── both merge to main ───────┘
Wave 2:  #C (files y,w)                          shares y with #A → must follow wave 1,
         rebased onto post-wave-1 main           branched after wave 1 merges
Wave 3:  #D+#E+#F  (all touch the same file)     hard-shared → ONE PR, not three
```

## Procedure

Create a TodoWrite item per step; work in order.

1. **Reconcile `main` with `origin`** and land any orchestration doc first, so every member inherits it.
2. **Build the wave DAG** (see *Sequencing* above) and confirm each wave is file-disjoint and width ~2 before spawning.
3. **One role-tagged team:** `impl-<issue#>`, `rev-<PR#>`. **Names never reused** (a bare-name `SendMessage` can hit a dead agentId). **Brief each member fully in its spawn prompt** — issue + worktree path + the edit/verify constraints (Bash-heredoc rule, "don't run the visual harness," "run git from the worktree"). A fresh teammate starts cold and inherits none of your context.
4. **Dual plan gate before ANY code** (members idle here until you relay approval): (a) you GREP every cited file:line/symbol/value against live source; (b) a project reviewer on the plan; (c) explicit human "approved".
5. **TDD implement.** Assert on `data-*` attributes, not class strings; regression floors, not exact counts; include the **prod build** in the verify gate.
6. **Full review loop on EVERY PR**, even trivial — the GitHub thread is the audit-trail deliverable. (Self-authored PRs can't be `--approve`d; reviewers `--comment` a verdict.)
7. **Orchestrator-owned visual sign-off** for visual issues: serve against a **writable copy** of the dev DB, per-branch Vite on unique ports. Probe hygiene: hold host background identical across branches; a uniform signed color offset means the harness, not the component.
8. **Between waves:** re-run `npm install` + build on the **`main` checkout** (worktree installs leave it stale).
9. **Merge on human go-ahead**, then **sweep the milestone for unclosed issues** — squash-merge does not auto-close issues lacking a per-issue `Closes #N`.
10. **Teardown:** `shutdown_request` members; GC a worktree only after confirming clean + merged + unlocked + not your cwd.

## More pitfalls (field-learned)

- **Fast-forward local `main` after every merge, before you build or validate on it.** `gh pr merge` updates `origin/main`, not your local checkout — validate the stale tree and you've proven nothing. (Bit twice: a build on pre-merge code; a teammate committing onto an un-synced local `main`.)
- **No CI? Run the full downstream suite at least once per milestone before declaring done.** Per-PR verification only covers touched files; a cross-package change (a core edit breaking a downstream test) slips straight through and surfaces waves later.
- **Mid-milestone discoveries → file a follow-up issue, don't scope-creep the in-flight PR.** Capture the new problem as its own issue and keep the wave scoped; folding unbounded work into a live PR breaks the file-disjoint plan (and the value-identity of a mechanical change).
- **Teardown: inspect before you GC.** A "dirty" worktree may be a regenerated artifact (a committed skeleton/boneyard capture), not real work — look before `--force`. Locks go stale after a milestone. In a shared repo with parallel jobs, a clean worktree may be another live session's — confirm it's yours/stale, and never remove your own cwd. (Squash-merge makes branch commits non-ancestors of `main`, so an `is-ancestor` test reads "unmerged" even when the work shipped — check PR state, not ancestry.)
- **The independent guards are part of the safety net — don't route around them.** The push-to-`main` guard and the permission classifier each caught orchestrator scope-creep this run (a direct-main push; an over-scoped delete). Route doc/skill commits through PRs, scope destructive ops to explicitly-named targets, and treat a denial as a signal to narrow — not to work around.

## Prefer self-revealing over guarded

When a defect class can re-enter, prefer a change that makes re-entry *visibly break* over a detector bolted onto a system that still accepts bad input (e.g. delete a duplicate shim so a stale reference fails to resolve, instead of a CI grep guard).

## References

- `docs/superpowers/specs/2026-05-28-print-r3-agent-team-orchestration.md` — full record + value evidence (what the gate caught before code).
- `project_agent_team_orchestration` and `team-monitor-experiment` memories — the idle/wake findings.
