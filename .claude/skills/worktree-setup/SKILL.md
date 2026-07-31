---
name: worktree-setup
description: Complete first-time setup for a fresh git worktree — npm install for the frontend and Pkg.instantiate for the Julia projects. Run immediately after git worktree add.
---

# worktree-setup

Automates the steps every fresh worktree needs before Julia or npm will work correctly.

## Procedure

```bash
# 1. Confirm we're in a worktree (or the main repo — the procedure is idempotent)
root=$(git rev-parse --show-toplevel)
echo "Worktree root: $root"

# 2. Install frontend dependencies
(cd "$root/packages/HimalayaUI/frontend" && npm install)

# 3. Instantiate the Julia projects. HimalayaUI must `develop` the local core
#    FIRST: it declares `[compat] Himalaya = "0.6"` and the registry's newest
#    published version is older, so a bare instantiate fails to resolve. That
#    failure is deliberate — see "develop core Himalaya" below.
julia --project="$root" -e 'using Pkg; Pkg.instantiate()'
julia --project="$root/packages/HimalayaUI" -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# 4. Verify HimalayaUI loads
julia --project="$root/packages/HimalayaUI" -e 'using HimalayaUI; println("HimalayaUI ok")'
```

## Required: develop core Himalaya from the worktree

This used to be optional. It isn't any more.

`Manifest.toml` is gitignored and `HimalayaUI` has no `[sources]` (that needs Julia 1.11; the package declares `julia = "1.9"`), so a bare `Pkg.instantiate()` resolves `Himalaya` **from the registry** — the published copy of the physics, not the one in this checkout. Ratio positions renumber when the series changes, so a mismatch is wrong physics, not a load error: that is how a registry core still carrying the forbidden `√11` could back a working tree that had removed it (#304).

`[compat] Himalaya = "0.6"` catches that only **across versions**. Edit `phaseratios` in your worktree without bumping the core version and the registry copy resolves cleanly while disagreeing with your code — exactly the original failure, back again. Point at the local core instead:

```bash
julia --project="$root/packages/HimalayaUI" -e 'using Pkg; Pkg.develop(path="../..")'
```

(From the repo root the path is `.`; from the HimalayaUI project directory it is `../..`.)

This rewrites the HimalayaUI Manifest to pin `Himalaya` at the local path. `Manifest.toml` is gitignored, so it affects only your worktree — which is exactly why the bound has to live in `Project.toml`.

**Any change to `phaseratios` must bump the core minor version and this bound in lockstep.**
