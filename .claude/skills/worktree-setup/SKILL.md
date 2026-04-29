---
name: worktree-setup
description: Complete first-time setup for a fresh git worktree — npm install + copy both Manifest.toml files from main so Himalaya core resolves to local v0.5.0 instead of registry v0.4.5. Run immediately after git worktree add.
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

# 3+4. Copy Manifest.toml files from the main worktree filesystem.
# Manifest.toml is gitignored (Julia convention) so `git show main:` won't work —
# we locate the main worktree via `git worktree list` instead.
main_root=$(git worktree list | awk '/\[main\]/{print $1}')
echo "Main worktree: $main_root"

src="$main_root/packages/HimalayaUI/Manifest.toml"
if [ -f "$src" ]; then
    cp "$src" "$root/packages/HimalayaUI/Manifest.toml"
    echo "Copied packages/HimalayaUI/Manifest.toml"
else
    echo "WARNING: $src not found — run Pkg.instantiate() in packages/HimalayaUI manually"
fi

src="$main_root/scripts/Manifest.toml"
if [ -f "$src" ]; then
    cp "$src" "$root/scripts/Manifest.toml"
    echo "Copied scripts/Manifest.toml"
else
    echo "scripts/Manifest.toml not in main worktree — run: julia --project=scripts -e 'using Pkg; Pkg.instantiate()' before make sysimage"
fi

# 5. Verify Julia loads HimalayaUI correctly
echo "Verifying HimalayaUI loads..."
julia --project="$root/packages/HimalayaUI" -e 'using HimalayaUI; println("HimalayaUI ok")'
```

If step 5 fails with a `findpeaks` signature error, the Manifest copy didn't take effect — check that `main` branch has the correct Manifest and retry.
