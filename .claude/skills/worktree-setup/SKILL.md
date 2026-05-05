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

# 3. Instantiate the Julia projects (resolves against the registry; Himalaya v0.5+
#    is published, so no Manifest copy from main is required).
julia --project="$root" -e 'using Pkg; Pkg.instantiate()'
julia --project="$root/packages/HimalayaUI" -e 'using Pkg; Pkg.instantiate()'

# 4. Verify HimalayaUI loads
julia --project="$root/packages/HimalayaUI" -e 'using HimalayaUI; println("HimalayaUI ok")'
```

## Optional: develop core Himalaya from the worktree

By default, `Pkg.instantiate()` resolves `Himalaya` from the registry. That's fine for any work that only touches `HimalayaUI` (backend or frontend). If you need to edit core `Himalaya` in this worktree and have HimalayaUI pick up your changes, run once:

```bash
julia --project="$root/packages/HimalayaUI" -e 'using Pkg; Pkg.develop(path="../..")'
```

This rewrites the HimalayaUI Manifest to pin `Himalaya` at `../..`. `Manifest.toml` is gitignored, so this affects only your local worktree.
