FRONTEND_DIR := packages/HimalayaUI/frontend
BUILD_DIR    := build

.PHONY: all frontend sysimage check-sysimage clean test-parallel

all: frontend sysimage

frontend:
	cd $(FRONTEND_DIR) && npm ci && npm run build

# A sysimage bakes an ABSOLUTE artifact path for every JLL it embeds, so it must
# be built against the same depot the runtime will use. Build it under a
# different `JULIA_DEPOT_PATH` and the binary cannot start at all — the first JLL
# `__init__` throws `Artifact "…" was not found by looking in the path …`, before
# any of our code runs.
#
# `bin/himalaya` takes that variable from a .env file, so the runtime depot is
# whatever that file says. Resolve the same file the same way here and the two
# cannot diverge; without it, `make sysimage` silently inherits the invoking
# user's default depot while the service runs against the deploy's.
#
# NOTE: must match the search order in bin/himalaya. Not cross-checked (same
# caveat as GROUPS below) — if you change one, change the other.
# Absolute paths on purpose: recipes run under /bin/sh, whose `.` searches PATH
# (not the cwd) for an operand containing no slash, so a bare `.env` is "not
# found". bin/himalaya sources "$REPO/.env" and sidesteps this; $(CURDIR) is the
# same idea.
ENV_CANDIDATES = $(HIMALAYA_ENV_FILE) $(CURDIR)/.env $(CURDIR)/packages/HimalayaUI/.env /etc/himalaya/.env

LOAD_ENV = set -a; for c in $(ENV_CANDIDATES); do \
	           if [ -f "$$c" ]; then . "$$c"; break; fi; \
	         done; set +a;

sysimage:
	@$(LOAD_ENV) \
	echo "Building sysimage with JULIA_DEPOT_PATH=$${JULIA_DEPOT_PATH:-<julia default>}"; \
	julia --project=scripts -e 'using Pkg; Pkg.instantiate()' && \
	julia --project=scripts scripts/build_sysimage.jl

check-sysimage:
	@$(LOAD_ENV) \
	stamp=$(BUILD_DIR)/julia_version; \
	depot_stamp=$(BUILD_DIR)/depot_path; \
	if [ ! -f $$stamp ]; then \
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
	if [ ! -f $$depot_stamp ]; then \
		echo "Sysimage OK (Julia $$current), but it predates depot stamping"; \
		echo "so the depot it was built against cannot be verified. If it fails"; \
		echo "to start with an \"Artifact ... was not found\" error, rebuild."; \
		exit 0; \
	fi; \
	stored_depot=$$(cat $$depot_stamp); \
	current_depot=$${JULIA_DEPOT_PATH:-}; \
	if [ "$$stored_depot" != "$$current_depot" ]; then \
		echo "Depot mismatch — this sysimage will fail to start here:"; \
		echo "  sysimage built with: $${stored_depot:-<julia default>}"; \
		echo "  runtime depot:       $${current_depot:-<julia default>}"; \
		echo "The embedded JLL artifact paths point at the build depot."; \
		echo "Run 'make sysimage' to rebuild against the runtime depot."; exit 1; \
	fi; \
	echo "Sysimage OK (Julia $$current, depot $${current_depot:-<julia default>})"

clean:
	rm -rf $(BUILD_DIR) $(FRONTEND_DIR)/dist

# NOTE: must match the bucket names in packages/HimalayaUI/test/runtests.jl
# (the GROUPS const). This convenience runner is NOT cross-checked against it;
# the authoritative serial path (GROUP=All) is drift-guarded in runtests.jl. If
# you add/rename a bucket there, update this line too or it won't run in parallel.
GROUPS := db pipeline routes events wire migration
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
