.PHONY: all clean test download bam2fastq compress
.ONESHELL:

SHELL   = /bin/bash
WGET    = wget
CONDA   = conda
MKDIR_P = mkdir -p
RM_RF   = rm -rf

CONDA_BASE          = $(shell conda info --base)
CONDA_ENV           = $(addsuffix /env, $(realpath .))
CONDA_URL           = https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
CONDA_ACTIVATE      = source $(CONDA_BASE)/etc/profile.d/conda.sh ; conda activate ; conda activate

SYNAPSE_SESSION     = $(HOME)/.synapseSession
BEDTOOLS_URL        = https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary

DATA_DIR            = data
DATASETS_TSV        = $(DATA_DIR)/datasets.tsv
DATASETS_TEST_TSV   = $(DATA_DIR)/datasets.test.tsv

DEFAULT_DIR         = $(DATA_DIR) # Change this line if you want to change the default directory
#TODO: Dynamic set doesn't work correctly
# INSTALL_INPUT	:= $(realpath $(shell read -p "Enter path to installation directory [build]:"$$'\n' && echo $$REPLY))
INSTALL_DIR         ?= $(or $(INSTALL_INPUT), $(strip $(DEFAULT_DIR)))
TEST_DIR            := $(INSTALL_DIR)/test

URL_FILES             := $(shell awk 'NR<=1{next} $$3 !~ /null/{printf "$(INSTALL_DIR)/%s/%s.%s\n", $$1, $$2, "url"}' $(DATASETS_TSV))
URL_TEST_FILES        := $(shell awk 'NR<=1{next} $$3 !~ /null/{printf "$(TEST_DIR)/%s/%s.%s\n", $$1, $$2, "url"}' $(DATASETS_TEST_TSV))
SYNIDS_FILES          := $(shell awk 'NR<=1{next} $$4 !~ /null/{printf "$(INSTALL_DIR)/%s/%s.%s\n", $$1, $$2, "synid"}' $(DATASETS_TSV))
SYNIDS_TEST_FILES     := $(shell awk 'NR<=1{next} $$4 !~ /null/{printf "$(TEST_DIR)/%s/%s.%s\n", $$1, $$2, "synid"}' $(DATASETS_TEST_TSV))

DATASETS_FROM_URL_FILES         :=  $(URL_FILES:.url=)
DATASETS_FROM_URL_TEST_FILES    :=  $(URL_TEST_FILES:.url=)
DATASETS_FROM_SYNID_FILES       :=  $(SYNIDS_FILES:.synid=)
DATASETS_FROM_SYNID_TEST_FILES  :=  $(SYNIDS_TEST_FILES:.synid=)

DATASETS_FILES        := $(DATASETS_FROM_URL_FILES) $(DATASETS_FROM_SYNID_FILES)
DATASETS_TEST_FILES   := $(DATASETS_FROM_URL_TEST_FILES) $(DATASETS_FROM_SYNID_TEST_FILES)

DATASETS_DIRS       := $(shell echo "$(dir $(DATASETS_FILES))" | tr ' ' '\n' | uniq)
DATASETS_TEST_DIRS  := $(shell echo "$(dir $(DATASETS_TEST_FILES))" | tr ' ' '\n' | uniq)

BAM_FILES           := $(shell find $(INSTALL_DIR)/ -type f -name '*.bam')
BAM2FASTQ_FILES     := $(BAM_FILES:.bam=.r1.fastq) $(BAM_FILES:.bam=.r2.fastq)
FASTQ2GZ_FILES      := $(BAM2FASTQ_FILES:.fastq=.fastq.gz)

ifeq (, $(shell which $(CONDA)))
$(error "No $(CONDA) in $PATH, check miniconda website to install it ($(CONDA_URL))")
endif


###############################################################################
#                                  GOALS
###############################################################################

all: download

download: $(DATASETS_FILES)

test: $(DATASETS_TEST_FILES)

bam2fastq: $(BAM2FASTQ_FILES)

compress: $(FASTQ2GZ_FILES)

clean:
	$(RM_RF) $(URL_FILES) $(URL_TEST_FILES)
	$(RM_RF) $(DATASETS_TEST_FILES) $(DATASETS_FILES) 
	$(RM_RF) $(TEST_DIR) $(DATASETS_DIRS)


###############################################################################
#                                 CONDA
###############################################################################

# Create the environment, install and initialize all the tools required
$(CONDA_ENV): environment.yml
	$(CONDA) env create -p $@ --file $<


###############################################################################
#                                 DOWNLOAD
###############################################################################
.INTERMEDIATE: $(URL_TEST_FILES) $(URL_FILES) $(SYNIDS_FILES) $(SYNIDS_TEST_FILES)
# TODO: Check file hash

%/:
	$(MKDIR_P) $(@D)

%.synid:
	awk '/$(notdir $*)/ {print $$4}' $< > $@

%.url: 
	awk '/$(notdir $*)/ {print $$3}' $< > $@

$(URL_TEST_FILES) $(SYNIDS_TEST_FILES): $(DATASETS_TEST_TSV) | $(DATASETS_TEST_DIRS)
$(URL_FILES) $(SYNIDS_FILES): $(DATASETS_TSV) | $(DATASETS_DIRS)

# If synapse credentials not saved before, login and save them with the OS keyring tool
$(SYNAPSE_SESSION): | $(CONDA_ENV)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	synapse login --rememberMe

# Download file from synapse if there is an URL in the tsv file
$(DATASETS_FROM_URL_TEST_FILES) $(DATASETS_FROM_URL_FILES): $(INSTALL_DIR)/%: | $(INSTALL_DIR)/%.url
	$(WGET) -O $@ $(shell cat $(word 1, $|))

# Download file from synapse if there is a synapse ID in the tsv file
$(DATASETS_FROM_SYNID_FILES) $(DATASETS_FROM_SYNID_TEST_FILES): $(INSTALL_DIR)/%: | $(INSTALL_DIR)/%.synid $(SYNAPSE_SESSION)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	syn_file=$$(synapse get --downloadLocation $(@D) $(shell cat $(word 1, $|)) | awk 'BEGIN {FS=" "} /Downloaded file/ {print $$3}')
	mv $(@D)/$${syn_file} $@


###############################################################################
#                                 bam2fastq
###############################################################################

# bedtools bamtofastq -i $< -fq $(word 1, $@)
%.r1.fastq %.r2.fastq: %.bam | $(CONDA_ENV)
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	gatk SamToFastq --INPUT=$< --FASTQ=$*.r1.fastq --SECOND_END_FASTQ=$*.r2.fastq

%.fastq.gz: %.fastq
	$(CONDA_ACTIVATE) $(CONDA_ENV)
	pigz $<
