FRONTEND_DIR := packages/HimalayaUI/frontend
BUILD_DIR    := build

.PHONY: all frontend sysimage check-sysimage clean test-parallel setup demo-data demo demo-empty demo-clean

all: frontend sysimage

frontend:
	cd $(FRONTEND_DIR) && npm ci && npm run build

sysimage:
	julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
	julia --project=scripts scripts/build_sysimage.jl

check-sysimage:
	@stamp=$(BUILD_DIR)/julia_version; \
	if [ ! -f $$stamp ]; then \
		echo "No sysimage found. Run 'make sysimage'."; exit 1; \
	fi; \
	stored=$$(cat $$stamp); \
	current=$$(julia -e 'print(VERSION)'); \
	if [ "$$stored" = "$$current" ]; then \
		echo "Sysimage OK (Julia $$current)"; \
	else \
		echo "Julia version mismatch:"; \
		echo "  sysimage built with: $$stored"; \
		echo "  current Julia:       $$current"; \
		echo "Run 'make sysimage' to rebuild."; exit 1; \
	fi

clean:
	rm -rf $(BUILD_DIR) $(FRONTEND_DIR)/dist

# NOTE: must match the bucket names in packages/HimalayaUI/test/runtests.jl
# (the GROUPS const). This convenience runner is NOT cross-checked against it;
# the authoritative serial path (GROUP=All) is drift-guarded in runtests.jl. If
# you add/rename a bucket there, update this line too or it won't run in parallel.
GROUPS := db pipeline routes routes2 events wire migration
test-parallel:
	@mkdir -p build
	@echo "Running $(words $(GROUPS)) groups in parallel..."
	@pids=""; rc=0; \
	for g in $(GROUPS); do \
		( GROUP=$$g HIMALAYA_SUITE_PARALLEL=1 julia --project=packages/HimalayaUI \
			-e 'using Pkg; Pkg.test("HimalayaUI")' > build/test-$$g.log 2>&1 ) & \
		pids="$$pids $$!"; \
	done; \
	for p in $$pids; do wait $$p || rc=1; done; \
	for g in $(GROUPS); do \
		echo "== $$g =="; grep -E "Test Summary" build/test-$$g.log | tail -1 || true; \
	done; \
	exit $$rc

# ---------------------------------------------------------------------------
# Local / cloud dev setup and the synthetic demo (docs/demo.md)
# ---------------------------------------------------------------------------

# First-time setup, per CLAUDE.md. `npm ci` (not `npm install`) so the lockfile
# isn't rewritten. On Julia >= 1.11 `Pkg.develop` also writes a [sources] entry
# into HimalayaUI's Project.toml and drops its comments; the (gitignored)
# Manifest already records the path dep, so put back the file as it was before
# (a copy, not `git checkout`, so uncommitted edits to it survive).
UI_PROJECT := packages/HimalayaUI/Project.toml
setup:
	cd $(FRONTEND_DIR) && npm ci
	julia --project=. -e 'using Pkg; Pkg.instantiate()'
	cp $(UI_PROJECT) $(UI_PROJECT).setup-bak; \
	trap 'mv $(UI_PROJECT).setup-bak $(UI_PROJECT)' EXIT; \
	julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Synthetic SSRL-style experiment, generated on demand and NEVER committed
# (scratch/ is gitignored). Override the target dir with DEMO_DIR=/tmp/...
DEMO_DIR    ?= scratch/demo-experiment
DEMO_DB_DIR ?= scratch/demo
demo-data:
	julia --project=packages/HimalayaUI packages/HimalayaUI/scripts/make_demo_data.jl $(DEMO_DIR)

# Backend + Vite with the demo experiment ingested (through the same API calls
# the New-experiment funnel makes). Generates the data first if it is missing.
demo:
	DEMO_DIR=$(DEMO_DIR) DEMO_DB_DIR=$(DEMO_DB_DIR) packages/HimalayaUI/scripts/demo.sh

# Backend + Vite on a fresh empty DB, to walk the ingest funnel by hand.
demo-empty:
	DEMO_DIR=$(DEMO_DIR) DEMO_DB_DIR=$(DEMO_DB_DIR) packages/HimalayaUI/scripts/demo.sh --empty

# Deletes DEMO_DIR only when it is a demo tree (has TRUTH.csv), so a mistyped
# DEMO_DIR pointing at real data is left alone.
demo-clean:
	@if [ -f "$(DEMO_DIR)/TRUTH.csv" ]; then rm -rf "$(DEMO_DIR)"; \
	elif [ -e "$(DEMO_DIR)" ]; then echo "$(DEMO_DIR) has no TRUTH.csv; not a demo tree, leaving it" >&2; fi
	rm -rf "$(DEMO_DB_DIR)" "$(DEMO_DB_DIR)-empty"
