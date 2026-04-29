---
name: review-response
description: Triage received PR review comments, implement fixes, and post a structured rebuttal comment in one step. Handles the author side of the review loop — complements review-pr (which handles the reviewer side). Usage: /review-response [PR number]
---

# review-response

Reads review comments on a PR, triages each point, implements fixes, and posts a rebuttal comment explaining what was addressed and why, including explicit rationale for anything intentionally skipped.

## Steps

1. **Identify the PR.** Use the arg if provided; otherwise use the current branch:
   ```bash
   gh pr view --comments          # current branch
   gh pr view <number> --comments # named PR
   ```

2. **Triage every point** into one of four buckets:

   | Bucket | Criteria | Action |
   |--------|----------|--------|
   | **Bug** | Incorrect behaviour, broken edge case, security issue | Fix before posting comment |
   | **Doc gap** | Missing or stale documentation, CLAUDE.md note, README entry | Fix — these are quick and always worth it |
   | **Suggestion** | Non-blocking improvement, refactor, style | Assess value vs. complexity; implement or skip with rationale |
   | **Praise** | Compliments on good decisions | Acknowledge in the comment; no code change |

3. **Verify reviewer suggestions before implementing.** A reviewer's proposed fix can itself contain bugs — read it critically rather than applying it verbatim. If you diverge from the suggestion, document exactly why in the rebuttal comment. Example: a reviewer suggests a symlink-resolution loop that prepends the dirname unconditionally, but `readlink` on macOS returns an absolute path, so prepending produces a double-path — the fix needs an absolute-vs-relative branch.

4. **Implement bugs and doc gaps.** For suggestions, implement if the value clearly exceeds the complexity cost; otherwise prepare a written rationale. Common rationales for skipping:
   - "Adds complexity for marginal benefit" (e.g. stamp-file Make prerequisites for a 1s operation)
   - "Correct but out of scope for this PR" (defer to a follow-up)
   - "Reviewer's suggested approach has an issue" (explain and use an alternative)

5. **Commit and push fixes** following the repo's commit style:
   ```bash
   git add <files>
   git commit -m "fix: address PR review — <summary of what changed>"
   git push
   ```

6. **Post the rebuttal comment:**
   ```bash
   gh pr comment <number> --body "<comment>"
   ```

   Structure the comment as one paragraph per reviewed point:
   - **Fixed items**: what changed, and if you diverged from the reviewer's suggestion, why
   - **Skipped items**: explicit rationale — "not worth it because X", "correct but adds complexity Y", "out of scope, tracked in Z"
   - **Praise**: a brief acknowledgment (one sentence)

   Post the comment even when all suggestions are skipped — it closes the loop and demonstrates that each point was considered.

7. **Print the returned comment URL.**

## Args

- `<number>` — PR number to respond to (default: current branch's open PR)

## Comment template

```
Thanks for the review. Here's what I addressed:

**Fixed: <title>**
<What changed and why. If you diverged from the reviewer's suggested approach, explain: "The reviewer suggested X, but Y has a bug on macOS because Z — used W instead.">

**Fixed: <title>**
<...>

**Skipped: <title>**
<Explicit rationale: complexity cost, scope, or alternative approach preferred.>

**<Praise point>**
<One-sentence acknowledgment.>
```

Omit sections that don't apply. At minimum always include one sentence per reviewed point.

## Gotchas

- **Reviewer suggestions can have bugs.** The PR #11 symlink fix is the reference case: the reviewer's while-loop used `$(cd "$(dirname "$SCRIPT")" && pwd)/$(readlink "$SCRIPT")` which produces `/tmp//absolute/path` when `readlink` returns an absolute target on macOS. Always test proposed fixes in the actual environment before committing.
- **Skipping without rationale reads as dismissive.** Even a one-sentence explanation ("adds ~20 lines for a 1s gain") is enough. The goal is showing the reviewer their feedback was read and weighed.
- **Don't conflate "won't fix" with "wrong direction."** If the reviewer identified a real issue but you disagree with their proposed solution, fix the issue your way and explain the difference — don't skip it entirely.
- **Post the comment last**, after all commits are pushed, so the comment accurately reflects the final state of the branch.
