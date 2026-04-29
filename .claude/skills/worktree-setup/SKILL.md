---
name: worktree-setup
description: Complete first-time setup for a fresh git worktree — npm install + copy both Manifest.toml files from main so Himalaya core resolves to local v0.5.0 instead of registry v0.4.5. Run immediately after git worktree add.
disable-model-invocation: true
---

# worktree-setup

Automates the three manual steps every fresh worktree needs before Julia or npm will work correctly.

## Why this matters

`Manifest.toml` is gitignored (Julia convention). A fresh worktree re-resolves dependencies against the registry and pulls the older published `Himalaya v0.4.5`, which has a different `findpeaks` signature from the local `v0.5.0`. Tests will fail silently with confusing errors until you copy the pinned Manifest from main. Same applies to `scripts/Manifest.toml` — needed before `make sysimage` will work.

## Procedure

```bash
# 1. Confirm we're in a worktree (not the main repo — no harm if run there, just redundant)
root=$(git rev-parse --show-toplevel)
echo "Worktree root: $root"

# 2. Install frontend dependencies
(cd "$root/packages/HimalayaUI/frontend" && npm install)

# 3. Copy packages/HimalayaUI/Manifest.toml from main
git show main:packages/HimalayaUI/Manifest.toml > "$root/packages/HimalayaUI/Manifest.toml"
echo "Copied packages/HimalayaUI/Manifest.toml from main"

# 4. Copy scripts/Manifest.toml from main (only if it exists there)
if git show main:scripts/Manifest.toml > /dev/null 2>&1; then
    git show main:scripts/Manifest.toml > "$root/scripts/Manifest.toml"
    echo "Copied scripts/Manifest.toml from main"
else
    echo "scripts/Manifest.toml not on main — run: julia --project=scripts -e 'using Pkg; Pkg.instantiate()' before make sysimage"
fi

# 5. Verify Julia loads HimalayaUI correctly
echo "Verifying HimalayaUI loads..."
julia --project="$root/packages/HimalayaUI" -e 'using HimalayaUI; println("HimalayaUI ok")'
```

If step 5 fails with a `findpeaks` signature error, the Manifest copy didn't take effect — check that `main` branch has the correct Manifest and retry.
