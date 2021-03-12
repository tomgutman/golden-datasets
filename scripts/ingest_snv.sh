#!/bin/bash

#EUCANCAN SNV & INDEL vcf handling

while getopts "t:s:i:f:o:n:c:h" option; do
    case "${option}" in
        t) truth=${OPTARG};;
        s) snv=${OPTARG};;
        i) indel=${OPTARG};;
        f) HGREF=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        n) SAMPLE_NAME=${OPTARG};;
        c) CONFIG=${OPTARG};;
        h)  echo "Usage:"
            echo "bash tmb_dragon.sh -t truth_file.vcf"
            echo "                   -s snv.vcf"
            echo "                   -i indel.vcf"
            echo "                   -f ref_fasta.fa"
            echo "                   -o /OUTPUT_DIR/PATH"
            echo "                   -n Sample name"
            echo "                   -c config_file.txt"
            exit
    esac
done

echo " "
echo -e "[General Information]:\n"
echo "Truth File:" $truth
echo "SNV vcf file:" $snv
echo "INDEL vcf file:" $indel
echo "Reference fasta file:" $FASTA
echo "output path": $OUTPUT_DIR
echo "sample Name:" $SAMPLE_NAME
echo "config file:" $CONFIG
echo " "

#Load Config file:
source ${CONFIG}

# Load conda env:
source activate test_benchmark

# Sorting vcf files
echo -e "[Running Information]: sorting vcf files\n"
bcftools sort $snv -o `basename $snv .vcf`".sort.vcf.gz" -O z
bcftools sort $indel -o `basename $indel .vcf`".sort.vcf.gz" -O z

# indexing sorted vcf files
echo -e "[Running Information]: indexing vcf files\n"
bcftools index `basename $snv .vcf`".sort.vcf.gz"
bcftools index `basename $indel .vcf`".sort.vcf.gz"

# Merging SNV.vcf & INDEL.vcf:
echo -e "[Running Information]: concatenate vcf files\n"
bcftools concat -a `basename $snv .vcf`".sort.vcf.gz" `basename $indel .vcf`".sort.vcf.gz" -o $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf"

# Normalizing vcfs:
echo -e "[Running Information]: normalizing test file\n"

#export $FASTA
export HGREF=/data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

for merge_file in $OUTPUT_DIR/*"_merge.vcf"; do
echo $merge_file
pre.py $merge_file `basename $merge_file .vcf`".norm.vcf.gz" -L --decompose --somatic --pass-only
done

echo -e "\n[Running Information]: normalizing truth file\n"

pre.py $truth `basename $truth .vcf`".norm.vcf.gz" -L --decompose --somatic --pass-only

# Running ingestion script:
echo -e "[Running Information]: Running ingestion_snv.py script \n"

python $SCRIPT_DIR/ingest_snv.py -samplename "SAMPLE" -o `basename $merge_file .vcf` `basename $merge_file .vcf`".norm.vcf.gz"

python $SCRIPT_DIR/ingest_snv.py -samplename "SAMPLE" -o `basename $truth .vcf` `basename $truth .vcf`".norm.vcf.gz"

# Running som.py:
echo -e "[Running Information]: Running som.py evaluation script\n"

som.py `basename $truth .vcf`".filtered.vcf" `basename $merge_file .vcf`".filtered.vcf" -o `basename $merge_file .vcf`"eval.txt" --verbose -N

echo -e "[Running Information]: script ended successfully\n"



# Get the samplename:
#bcftools query -l test_merge.norm.vcf.gz

# example command:

#bash ../golden-datasets-fork-tom/golden-datasets/scripts/ingest_snv.sh -t test_truth.vcf -s results_dream/CURIE/insilico_1_snvs.vcf -i results_dream/CURIE/insilico_1_indels.vcf -o . -c test_ingestions_sh/config.txt  -n test -f /data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

#handle header problems with truth files:
# Add the contig info by hand with vim - with chr (example: chr1,chr2)
# awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' truth.snvs.synthetic.challenge.set2.vcf > test_truth.vcf
