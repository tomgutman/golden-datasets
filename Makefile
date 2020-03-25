.PHONY: all clean test download bam2fastq

SHELL   = /bin/bash
WGET    = wget
MKDIR_P = mkdir -p
RM_RF   = rm -rf

BEDTOOLS_URL        = https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary
DATA_DIR            = data
DATASETS_TSV        = $(DATA_DIR)/datasets.tsv
DATASETS_TEST_TSV   = $(DATA_DIR)/datasets.test.tsv

DEFAULT_DIR         = $(DATA_DIR) # Change this line if you want to change the default directory
#TODO: Dynamic set doesn't work correctly
# INSTALL_INPUT	:= $(realpath $(shell read -p "Enter path to installation directory [build]:"$$'\n' && echo $$REPLY))
INSTALL_DIR         ?= $(or $(INSTALL_INPUT), $(strip $(DEFAULT_DIR)))
TEST_DIR            := $(INSTALL_DIR)/test

DATASETS_URLS       := $(shell awk 'NR>1 {print $$2}' $(DATASETS_TSV))
DATASETS_TEST_URLS  := $(shell awk 'NR>1 {print $$2}' $(DATASETS_TEST_TSV))
DATASETS_TAGS       := $(join $(addsuffix /,$(shell awk 'NR>1 {print $$1}' $(DATASETS_TSV))), $(notdir $(DATASETS_URLS)))
DATASETS_TEST_TAGS  := $(join $(addsuffix /,$(shell awk 'NR>1 {print $$1}' $(DATASETS_TEST_TSV))), $(notdir $(DATASETS_TEST_URLS)))
DATASETS_FILES      := $(addprefix $(INSTALL_DIR)/, $(DATASETS_TAGS))
DATASETS_TEST_FILES := $(addprefix $(TEST_DIR)/, $(DATASETS_TEST_TAGS))

DATASETS_DIRS       := $(shell echo "$(dir $(DATASETS_FILES))" | tr ' ' '\n' | uniq)
DATASETS_TEST_DIRS  := $(shell echo "$(dir $(DATASETS_TEST_FILES))" | tr ' ' '\n' | uniq)

URL_FILES           := $(addsuffix .url, $(DATASETS_FILES))
URL_TEST_FILES      := $(addsuffix .url, $(DATASETS_TEST_FILES))

BAM_FILES           := $(wildcard $(INSTALL_DIR)/*.bam)
BAM2FASTQ_FILES     := $(BAM_FILES:.bam=.fastq)

###############################################################################
#                                  GOALS
###############################################################################

all: download

download: $(DATASETS_FILES)

test: $(DATASETS_TEST_FILES)

bam2fastq: $(BAM2FASTQ_FILES)

clean:
	$(RM_RF) $(URL_FILES) $(URL_TEST_FILES)
	$(RM_RF) $(DATASETS_TEST_FILES) $(DATASETS_FILES) 
	$(RM_RF) $(TEST_DIR) $(DATASETS_DIRS)

###############################################################################
#                                 DOWNLOAD
###############################################################################
.INTERMEDIATE: $(URL_TEST_FILES) $(URL_FILES)


%/:
	$(MKDIR_P) $(@D)

%.url: 
	awk '/$(notdir $*)/ {print $$2}' $< > $@

$(URL_TEST_FILES): $(DATASETS_TEST_TSV) | $(DATASETS_TEST_DIRS)
$(URL_FILES): $(DATASETS_TSV) | $(DATASETS_DIRS)


$(INSTALL_DIR)/%.bam $(INSTALL_DIR)/%.bai $(INSTALL_DIR)/%.bz2 $(INSTALL_DIR)/%.fai $(INSTALL_DIR)/%.fastq:
	$(WGET) -O $@ $(shell cat $(word 1, $|))

$(DATASETS_TEST_FILES): $(INSTALL_DIR)/test/%: | $(INSTALL_DIR)/test/%.url
$(DATASETS_FILES): $(INSTALL_DIR)/%: | $(INSTALL_DIR)/%.url

# TODO: Check file hash



###############################################################################
#                                 bam2fastq
###############################################################################


bedtools:
	@echo "Install bedtools requirement"
	@wget $(BEDTOOLS_URL);
	@mv bedtools.static.binary bedtools;
	@chmod a+x bedtools


%.fastq: %.bam bedtools
	bedtools bamtofastq -i $< -fq $(word 1, $@)