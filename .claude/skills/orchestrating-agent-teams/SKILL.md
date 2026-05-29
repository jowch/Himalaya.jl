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

## Procedure

Create a TodoWrite item per step; work in order.

1. **Reconcile `main` with `origin`** and land any orchestration doc first, so every member inherits it.
2. **Build a cluster-aware, file-disjoint wave DAG.** Group by shared file; each wave is internally file-disjoint; later waves branch off `main` after earlier waves merge. **Keep wave width ~2** — you are the review-round serialization point; bound width by nudge bandwidth.
3. **One role-tagged team:** `impl-<issue#>`, `rev-<PR#>`. **Names never reused** (a bare-name `SendMessage` can hit a dead agentId).
4. **Dual plan gate before ANY code** (members idle here until you relay approval): (a) you GREP every cited file:line/symbol/value against live source; (b) a project reviewer on the plan; (c) explicit human "approved".
5. **TDD implement.** Assert on `data-*` attributes, not class strings; regression floors, not exact counts; include the **prod build** in the verify gate.
6. **Full review loop on EVERY PR**, even trivial — the GitHub thread is the audit-trail deliverable. (Self-authored PRs can't be `--approve`d; reviewers `--comment` a verdict.)
7. **Orchestrator-owned visual sign-off** for visual issues: serve against a **writable copy** of the dev DB, per-branch Vite on unique ports. Probe hygiene: hold host background identical across branches; a uniform signed color offset means the harness, not the component.
8. **Between waves:** re-run `npm install` + build on the **`main` checkout** (worktree installs leave it stale).
9. **Merge on human go-ahead**, then **sweep the milestone for unclosed issues** — squash-merge does not auto-close issues lacking a per-issue `Closes #N`.
10. **Teardown:** `shutdown_request` members; GC a worktree only after confirming clean + merged + unlocked + not your cwd.

## Prefer self-revealing over guarded

When a defect class can re-enter, prefer a change that makes re-entry *visibly break* over a detector bolted onto a system that still accepts bad input (e.g. delete a duplicate shim so a stale reference fails to resolve, instead of a CI grep guard).

## References

- `docs/superpowers/specs/2026-05-28-print-r3-agent-team-orchestration.md` — full record + value evidence (what the gate caught before code).
- `project_agent_team_orchestration` and `team-monitor-experiment` memories — the idle/wake findings.
