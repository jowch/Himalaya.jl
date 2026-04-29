FRONTEND_DIR := packages/HimalayaUI/frontend
BUILD_DIR    := build

.PHONY: all frontend sysimage check-sysimage clean

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
