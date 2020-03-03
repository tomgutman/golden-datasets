.PHONY: all clean test
.ONESHELL:

SHELL = /bin/bash
WGET = wget
MKDIR_P = mkdir -p

BEDTOOLS_URL = https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary

INSTALL_DIR := build
INSTALL_DIR ?= $(shell read -s -p "Enter path to installation directory:" DIR && echo $$DIR)

DATASETS_TSV = data/datasets.tsv
DATASETS_URLS = $(shell awk 'NR>1 {print $$2}' $(DATASETS_TSV))

DATASETS_FILES := $(addprefix $(INSTALL_DIR)/, $(notdir $(DATASETS_URLS)))

BAM_FILES := $(filter %.bam,$(DATASETS_FILES))
BAI_FILES := $(filter %.bai,$(DATASETS_FILES))
FASTQ_FILES := $(filter %.fastq,$(DATASETS_FILES))

BAM2FASTQ_FILES := $(BAM_FILES:.bam=.fastq)

INTURL_FILES := $(addsuffix .url, $(DATASETS_FILES))



# 1. DOWNLOAD FILES
$(INSTALL_DIR):
	@echo "Make installation directory"
	@$(MKDIR_P) $@


#TODO : Add subfolders for each dataset (field 1 in tsv)
$(INSTALL_DIR)/%.url: $(DATASETS_TSV) | $(INSTALL_DIR)
	awk '/$*/ {print $$2}' $^ > $@


# TODO: Check MD5
$(INSTALL_DIR)/%: | $(INSTALL_DIR)/%.url
	@$(WGET) -O $@ $(shell cat $|) 


# 2. BAM 2 FASTQ BAM FILES
bedtools:
	@echo "Install bedtools requirement"
	@wget $(BEDTOOLS_URL);
	@mv bedtools.static.binary bedtools;
	@chmod a+x bedtools


%.bam: %.fastq bedtools
	bedtools bamtofastq -i $@ -fq $(word 1, $<)



# GENERIC RULES
test: $(DATASETS_FILES) | $(INSTALL_DIR)


all:
	@echo "***"


clean:
	@echo "Clean commands"