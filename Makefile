TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime

APP_SERVICE = app_service

WRAP_PYTHON_TOOL = wrap_python3
WRAP_PYTHON_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PYTHON3_TOOL).sh

SRC_SERVICE_PYTHON = $(wildcard scripts/*.py)
BIN_SERVICE_PYTHON = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PYTHON))))
DEPLOY_SERVICE_PYTHON = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PYTHON))))

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

SRC_R = $(wildcard scripts/*.R)
BIN_R = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_R))))
DEPLOY_R = $(addprefix $(TARGET)/bin/,$(basename $(notdir $(SRC_R))))

CLIENT_TESTS = $(wildcard t/client-tests/*.t)
SERVER_TESTS = $(wildcard t/server-tests/*.t)
PROD_TESTS = $(wildcard t/prod-tests/*.t)

TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT) --define kb_service_dir=$(SERVICE_DIR) \
	--define kb_sphinx_port=$(SPHINX_PORT) --define kb_sphinx_host=$(SPHINX_HOST) \
	--define kb_starman_workers=$(STARMAN_WORKERS) \
	--define kb_starman_max_requests=$(STARMAN_MAX_REQUESTS)

all: bin

local_tools: $(BIN_DIR)/diffdock_run $(BIN_DIR)/count-pdb-residues

# Simplified DiffDock wrapper (replaces run_local_docking)
$(BIN_DIR)/diffdock_run: bvbrc_docking/diffdock_run.py
	export KB_CONDA_BASE=$(BVDOCK_CONDA_BASE); \
	export KB_CONDA_ENV=$(BVDOCK_ENV); \
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

# Legacy wrapper (kept for backward compatibility)
$(BIN_DIR)/run_local_docking: bvbrc_docking/run_local_docking.py
	export KB_CONDA_BASE=$(BVDOCK_CONDA_BASE); \
	export KB_CONDA_ENV=$(BVDOCK_ENV); \
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/count-pdb-residues: bvbrc_docking/count-pdb-residues.py
	export KB_CONDA_BASE=$(BVDOCK_CONDA_BASE); \
	export KB_CONDA_ENV=$(BVDOCK_ENV); \
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@
$(BIN_DIR)/check_input_smile_strings: scripts/check_input_smile_strings.py
	export KB_CONDA_BASE=$(BVDOCK_CONDA_BASE); \
	export KB_CONDA_ENV=$(BVDOCK_ENV); \
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

deploy-local-tools:
	if [ "$(KB_OVERRIDE_TOP)" != "" ] ; then sbase=$(KB_OVERRIDE_TOP) ; else sbase=$(TARGET); fi; \
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PYTHON_PATH=$(TARGET)/lib ; \
	export KB_CONDA_BASE=$(BVDOCK_CONDA_BASE); \
	export KB_CONDA_ENV=$(BVDOCK_ENV); \
	for script in diffdock_run count-pdb-residues check_input_smile_strings ; do \
	    cp bvbrc_docking/$$script.py $(TARGET)/pybin; \
	    $(WRAP_PYTHON_SCRIPT) "$$sbase/pybin/$$script.py" $(TARGET)/bin/$$script; \
	done

bin: $(BIN_PERL) $(BIN_SERVICE_PERL) $(BIN_R) $(BIN_SERVICE_PYTHON) local_tools

deploy: deploy-all
deploy-all: deploy-client 
deploy-client: deploy-libs deploy-scripts deploy-docs deploy-local-tools

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs deploy-local-tools

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs:

clean:

include $(TOP_DIR)/tools/Makefile.common.rules

#
# This is a little ugly, but it works for now. Because lib/bvbrc_docking is a symlink, when we
# do the install we need to copy the link data. It's possible we should always do this in the default rule
#
# We place this under the include of Makefile.common.rules because we are overriding
# behavior defined there.
#

deploy-libs:
	rm -rf $(TARGET)/lib/bvbrc_docking
	rsync --copy-links --exclude '*.bak*' -arv lib/* $(TARGET)/lib

