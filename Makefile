.PHONY: all clean test
.ONESHELL:
.SUFFIXES:
.SUFFIXES: .bam .bai
# .SECONDEXPANSION:


SHELL	= /bin/bash
WGET	= wget
MKDIR_P	= mkdir -p
RM_RF	= rm -rf
EMPTY	=

BEDTOOLS_URL 	= https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary
DATA_DIR		= data
DATASETS_TSV	= $(DATA_DIR)/datasets.test.tsv

DEFAULT_DIR		= $(DATA_DIR) # Change this line if you want to change the default directory
#TODO: Dynamic set doesn't work correctly
# INSTALL_INPUT	:= $(realpath $(shell read -p "Enter path to installation directory [build]:"$$'\n' && echo $$REPLY))
INSTALL_DIR		= $(or $(INSTALL_INPUT), $(strip $(DEFAULT_DIR)))

DATASETS_URLS	:= $(shell awk 'NR>1 {print $$2}' $(DATASETS_TSV))
DATASETS_TAGS 	:= $(join $(addsuffix /,$(shell awk 'NR>1 {print $$1}' $(DATASETS_TSV))), $(notdir $(DATASETS_URLS)))
DATASETS_FILES	:= $(addprefix $(INSTALL_DIR)/, $(DATASETS_TAGS))

BAM_FILES		:= $(filter %.bam,$(DATASETS_FILES))
BAI_FILES		:= $(filter %.bai,$(DATASETS_FILES))
FASTQ_FILES		:= $(filter %.fastq,$(DATASETS_FILES))
BAM2FASTQ_FILES	:= $(BAM_FILES:.bam=.fastq)

URL_FILES	:= $(addsuffix .url, $(value DATASETS_FILES))



#TODO : Add subfolders for each dataset (field 1 in tsv)
$(URL_FILES): $(DATASETS_TSV)
	$(MKDIR_P) $(@D)
	awk '/$(basename $(@F))/ {print $$2}' $< > $@


$(DATASETS_FILES): | $(URL_FILES)
	$(WGET) -O $@ $(shell cat $(word 1, $|))


%.bam: | %.url
%.bam.bai: | %.url


# TODO: Check MD5


# 2. BAM 2 FASTQ BAM FILES
bedtools:
	@echo "Install bedtools requirement"
	@wget $(BEDTOOLS_URL);
	@mv bedtools.static.binary bedtools;
	@chmod a+x bedtools


%.bam: %.fastq | bedtools
	bedtools bamtofastq -i $@ -fq $(word 1, $<)


# GENERIC RULES
.PHONY: download
download: $(DATASETS_FILES)


test: $(DATASETS_FILES)


all:


clean:
	$(RM_RF) $(DATASETS_FILES)
	$(RM_RF) $(dir $(DATASETS_FILES))