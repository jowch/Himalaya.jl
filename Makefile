FRONTEND_DIR := packages/HimalayaUI/frontend
BUILD_DIR    := build

.PHONY: all frontend sysimage check-sysimage deploy clean test-parallel setup dev-setup demo-data demo-large-data demo demo-large demo-empty demo-clean verify-truth lint

all: frontend sysimage

frontend:
	cd $(FRONTEND_DIR) && npm ci && npm run build

lint:
	cd $(FRONTEND_DIR) && npm run lint:design && npm run lint:ratchet && npm run lint

SHELL := /bin/bash

LOAD_ENV = set -a; for c in "$${HIMALAYA_ENV_FILE:-}" "$(CURDIR)/.env" "$(CURDIR)/packages/HimalayaUI/.env" /etc/himalaya/.env; do \
	           if [ -n "$$c" ] && [ -f "$$c" ]; then . "$$c"; break; fi; \
	         done; set +a;

JULIA = $(LOAD_ENV) julia

DEPOT_PROBE = 'd = first(DEPOT_PATH); print(ispath(d) ? realpath(d) : rstrip(normpath(d), collect("/")))'

sysimage:
	@$(JULIA) -e 'print("Building sysimage with depot ")' -e $(DEPOT_PROBE) -e 'println()'
	$(JULIA) --project=scripts -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=scripts scripts/build_sysimage.jl $(BUILD_DIR)

check-sysimage:
	@$(LOAD_ENV) \
	stamp=$(BUILD_DIR)/julia_version; \
	depot_stamp=$(BUILD_DIR)/depot_path; \
	if [ ! -f $(BUILD_DIR)/himalaya.so ] || [ ! -f $$stamp ]; then \
		echo "No sysimage found. Run 'make sysimage'."; exit 1; \
	fi; \
	stored=$$(cat $$stamp); \
	current=$$(julia -e 'print(VERSION)'); \
	if [ "$$stored" != "$$current" ]; then \
		echo "Julia version mismatch:"; \
		echo "  sysimage built with: $$stored"; \
		echo "  current Julia:       $$current"; \
		echo "Run 'make sysimage' to rebuild."; exit 1; \
	fi; \
	current_depot=$$(julia -e $(DEPOT_PROBE)); \
	if [ -f $$depot_stamp ]; then \
		stored_depot=$$(cat $$depot_stamp); \
		if [ "$$stored_depot" != "$$current_depot" ]; then \
			echo "Depot mismatch. This sysimage will fail to start here:"; \
			echo "  sysimage built with: $$stored_depot"; \
			echo "  runtime depot:       $$current_depot"; \
			echo "Run 'make sysimage' to rebuild against the runtime depot."; exit 1; \
		fi; \
	fi; \
	if ! julia --project=packages/HimalayaUI --sysimage=$(BUILD_DIR)/himalaya.so --startup-file=no -e 'exit()'; then \
		echo "The sysimage failed to start with depot $$current_depot."; \
		echo "Run 'make sysimage' to rebuild against the runtime depot."; exit 1; \
	fi; \
	echo "Sysimage OK (Julia $$current, depot $$current_depot)"

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
	@$(LOAD_ENV) pids=""; rc=0; \
	for g in $(GROUPS); do \
		( GROUP=$$g julia --project=packages/HimalayaUI \
			-e 'using Pkg; Pkg.test("HimalayaUI"; julia_args=["--threads=4,1"])' > build/test-$$g.log 2>&1 ) & \
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
UI_PREFS   := packages/HimalayaUI/LocalPreferences.toml
setup:
	cd $(FRONTEND_DIR) && npm ci
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'
	cp $(UI_PROJECT) $(UI_PROJECT).setup-bak; \
	trap 'mv $(UI_PROJECT).setup-bak $(UI_PROJECT)' EXIT; \
	$(JULIA) --project=packages/HimalayaUI -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
	$(JULIA) --project=packages/HimalayaUI -e 'using Pkg; Pkg.precompile()'
	$(JULIA) --project=packages/HimalayaPeakModel -e 'using Pkg; Pkg.instantiate()'
	@! grep -qs 'precompile_workload *= *false' $(UI_PREFS) || \
		echo "Note: $(UI_PREFS) turns the route warm-up off. Delete it on a server."

HIMALAYA_SERVICE ?= himalaya
HIMALAYA_RESTART ?= sudo systemctl restart $(HIMALAYA_SERVICE)
NEXT_BUILD_DIR   := $(BUILD_DIR)/next
RESOLVE_DIST = $(LOAD_ENV) dist="$${HIMALAYA_FRONTEND_DIST:-$(FRONTEND_DIR)/dist}"; \
	case "$$dist" in /*) ;; *) dist="$(CURDIR)/$$dist" ;; esac; dist="$${dist%/}";

deploy:
	@case '$(HIMALAYA_RESTART)' in sudo\ *) \
		sudo -n -l $(wordlist 2,99,$(HIMALAYA_RESTART)) >/dev/null 2>&1 || sudo -v ;; esac
	$(MAKE) setup
	@$(RESOLVE_DIST) rm -rf "$$dist.next"; \
	cd $(FRONTEND_DIR) && npm run build -- --outDir "$$dist.next"
	$(MAKE) sysimage BUILD_DIR=$(NEXT_BUILD_DIR)
	$(MAKE) check-sysimage BUILD_DIR=$(NEXT_BUILD_DIR)
	@$(RESOLVE_DIST) \
	if [ ! -f "$$dist.next/index.html" ]; then echo "No frontend build at $$dist.next; nothing swapped."; exit 1; fi; \
	rm -rf "$$dist.old"; \
	if [ -d "$$dist" ]; then mv "$$dist" "$$dist.old" || exit 1; fi; \
	if ! mv "$$dist.next" "$$dist"; then \
		[ ! -d "$$dist.old" ] || mv "$$dist.old" "$$dist"; \
		echo "Could not move the new frontend into $$dist; the old one is back."; exit 1; \
	fi; \
	if ! mv $(NEXT_BUILD_DIR)/himalaya.so $(BUILD_DIR)/himalaya.so; then \
		rm -rf "$$dist"; [ ! -d "$$dist.old" ] || mv "$$dist.old" "$$dist"; \
		echo "Could not move the new sysimage into $(BUILD_DIR); the old frontend is back."; exit 1; \
	fi; \
	mv $(NEXT_BUILD_DIR)/julia_version $(NEXT_BUILD_DIR)/depot_path $(BUILD_DIR)/; \
	rm -rf "$$dist.old"
	$(HIMALAYA_RESTART)

dev-setup:
	@test -f $(UI_PREFS) || printf '[HimalayaUI]\nprecompile_workload = false\n' > $(UI_PREFS)
	$(MAKE) setup

# Synthetic SSRL-style experiment, generated on demand and NEVER committed
# (scratch/ is gitignored). Override the target dir with DEMO_DIR=/tmp/...
DEMO_DIR    ?= scratch/demo-experiment
DEMO_DB_DIR ?= scratch/demo
demo-data:
	$(JULIA) --project=packages/HimalayaUI packages/HimalayaUI/scripts/make_demo_data.jl $(DEMO_DIR) $(DEMO_DATA_FLAGS)

# Realistic scale: DEMO_LARGE_DIR/large_01 ... large_10, 300 samples each, with
# images shared from a pool (docs/demo.md "Large mode"). Takes ~35 s, ~1.4 GB.
DEMO_LARGE_DIR    ?= scratch/demo-large
DEMO_LARGE_DB_DIR ?= scratch/demo-large-db
demo-large-data:
	$(JULIA) --project=packages/HimalayaUI packages/HimalayaUI/scripts/make_demo_data.jl $(DEMO_LARGE_DIR) --large

# Generate (if missing) and ingest every large experiment, then start the app.
demo-large:
	$(LOAD_ENV) DEMO_DIR=$(DEMO_LARGE_DIR) DEMO_DB_DIR=$(DEMO_LARGE_DB_DIR) packages/HimalayaUI/scripts/demo.sh --large

# Backend + Vite with the demo experiment ingested (through the same API calls
# the New-experiment funnel makes). Generates the data first if it is missing.
demo:
	$(LOAD_ENV) DEMO_DIR=$(DEMO_DIR) DEMO_DB_DIR=$(DEMO_DB_DIR) packages/HimalayaUI/scripts/demo.sh

# Backend + Vite on a fresh empty DB, to walk the ingest funnel by hand.
demo-empty:
	$(LOAD_ENV) DEMO_DIR=$(DEMO_DIR) DEMO_DB_DIR=$(DEMO_DB_DIR) packages/HimalayaUI/scripts/demo.sh --empty

verify-truth:
	$(LOAD_ENV) packages/HimalayaUI/scripts/verify_truth.sh

# Deletes DEMO_DIR only when it is a demo tree (has TRUTH.csv), and DEMO_LARGE_DIR
# only when it holds nothing but large_NN folders and the completion marker, so a
# mistyped path pointing at real data is left alone.
demo-clean:
	@if [ -f "$(DEMO_DIR)/TRUTH.csv" ]; then rm -rf "$(DEMO_DIR)"; \
	elif [ -e "$(DEMO_DIR)" ]; then echo "$(DEMO_DIR) has no TRUTH.csv; not a demo tree, leaving it" >&2; fi
	rm -rf "$(DEMO_DB_DIR)" "$(DEMO_DB_DIR)-empty"
	rm -rf "$(DEMO_LARGE_DB_DIR)"
	@if [ -e "$(DEMO_LARGE_DIR)" ] || [ -L "$(DEMO_LARGE_DIR)" ]; then \
		$(JULIA) --project=packages/HimalayaUI --startup-file=no packages/HimalayaUI/scripts/make_demo_data.jl --clean-large "$(DEMO_LARGE_DIR)"; fi
