#!/bin/bash

#EUCANCAN SNV & INDEL vcf handling

while getopts "t:s:i:v:f:o:n:c:h" option; do
    case "${option}" in
        t) truth=${OPTARG};;
        s) snv=${OPTARG};;
        i) indel=${OPTARG};;
        v) sv=${OPTARG};;
        u) truth_sv=${OPTARG};;
        f) HGREF=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        n) SAMPLE_NAME=${OPTARG};;
        c) CONFIG=${OPTARG};;
        h)  echo "Usage:"
            echo "bash tmb_dragon.sh -t truth_file.vcf"
            echo "                   -s snv.vcf"
            echo "                   -i indel.vcf"
            echo "                   -v sv.vcf"
            echo "                   -u truth_sv.vcf"
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
echo "SV vcf file:" $sv
echo "TRUTH SV file:" $truth_sv
echo "Reference fasta file:" $FASTA
echo "output path": $OUTPUT_DIR
echo "sample Name:" $SAMPLE_NAME
echo "config file:" $CONFIG
echo " "

#Load Config file:
source ${CONFIG}

# Load conda env:
#conda env create -n eucancan -f golden-datasets/scripts/environment.yml 

source activate eucancan

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
bcftools concat -a `basename $snv .vcf`".sort.vcf.gz" `basename $indel .vcf`".sort.vcf.gz" -O z -o $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz"

# Preparing for normalization
echo -e "[Running Information]: preparing for normalizing\n"

#export $FASTA
export HGREF=/data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

bcftools index $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz"
bcftools index $truth

multimerge $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz" -r $HGREF -f true -o $OUTPUT_DIR/$SAMPLE_NAME"_merge.prep.vcf" --process-full=1
multimerge $truth -r $HGREF -o `basename $truth .vcf.gz`".prep.vcf.gz" --process-full=1

# Normalizing vcfs:
echo -e "[Running Information]: normalizing test file\n"
merge_file=$OUTPUT_DIR/$SAMPLE_NAME"_merge.prep.vcf"

pre.py $merge_file `basename $merge_file _merge.prep.vcf`".norm.vcf.gz" -L --decompose --somatic --pass-only

echo -e "\n[Running Information]: normalizing truth file\n"

truth=`basename $truth .vcf.gz`".prep.vcf.gz"
pre.py $truth `basename $truth .vcf.gz`".norm.vcf.gz" -L --decompose --somatic --pass-only

# Running ingestion script:
echo -e "[Running Information]: Running ingestion_snv.py script \n"

python $SCRIPT_DIR/ingest_snv.py -samplename "SAMPLE" -o `basename $merge_file .vcf` `basename $merge_file _merge.prep.vcf`".norm.vcf.gz"

python $SCRIPT_DIR/ingest_snv.py -samplename "SAMPLE" -o `basename $truth .vcf.gz` `basename $truth .vcf.gz`".norm.vcf.gz"

# Running som.py:
echo -e "[Running Information]: Running som.py evaluation script\n"

som.py `basename $truth .vcf.gz`".filtered.vcf" `basename $merge_file .vcf`".filtered.vcf" -o `basename $merge_file .vcf`"eval.txt" --verbose -N

echo -e "[Running Information]: script ended successfully\n"

# Running SV ingestion script:

echo -e "[Running Information]: Running ingestion SV script\n"

# Cleaning
echo -e "[Running Information]: Cleaning files\n"
rm *{.tbi,csi,json}

#run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/Hartwig/5.19/COLO829v003T.purple.sv.vcf.gz" "COLO829v003T" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/hartwig.csv"

# Get the samplename:
#bcftools query -l test_merge.norm.vcf.gz

# example command:

#bash ../golden-datasets-fork-tom/golden-datasets/scripts/ingest_snv.sh -t test_truth.vcf -s results_dream/CURIE/insilico_1_snvs.vcf -i results_dream/CURIE/insilico_1_indels.vcf -o . -c test_ingestions_sh/config.txt  -n test -f /data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

# test with curie colo files:
#bash ../golden-datasets-fork-tom/golden-datasets/scripts/ingest_snv.sh -t ../Benchmark/test_truth.vcf -s curie_colo829_snps.sample.vcf -i curie_colo829_indels.sample.vcf -o . -c ../Benchmark/test_ingestions_sh/config.txt  -n COLO829T -f /data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

# bcftools norm -m -any ./COLO829T_merge.vcf > test.multi.vcf
# pre.py test.multi.vcf test.multi.norm.vcf -L --decompose --somatic --pass-only
# pre.py COLO829T_merge.vcf test_output.vcf -L --decompose --somatic --pass-only
#
# hap.py test.multi.vcf test.multi.norm.vcf -L --decompose --somatic --pass-only
# hap.py COLO829T_merge.vcf ../Benchmark/test_truth.vcf -o test.multi.norm.vcf -L --decompose --somatic --pass-only --unhappy
#
# /data/users/tgutman/.conda/envs/test_benchmark/bin/vcfcheck COLO829T_merge.vcf --check-bcf-errors 1
# /data/users/tgutman/.conda/envs/test_benchmark/bin/vcfcheck test.multi.vcf --check-bcf-errors 1

#handle header problems with truth files:
# Add the contig info by hand with vim - with chr (example: chr1,chr2)
# awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' truth.snvs.synthetic.challenge.set2.vcf > test_truth.vcf

#handle mulisample files:

# bcftools view -c1 -Oz -s COLO829T -o curie_colo829_indels.sample.vcf curie_colo829_indels.vcf
# bcftools view -c1 -Oz -s COLO829T -o curie_colo829_snps.sample.vcf curie_colo829_snps.vcf
#
# gunzip $OUTPUT_DIR/$dataset/$sample_vcf"_multisample.hg19_multianno.rec.vcf.gz"
#
# ${GATK_DIR}/gatk --java-options "-Xmx20g" SelectVariants \
#     -O test.multi.vcf \
#     -V COLO829T_merge.vcf \
#     --select-type-to-include SNP --select-type-to-include INDEL --exclude-filtered --select-type-to-exclude MNP
