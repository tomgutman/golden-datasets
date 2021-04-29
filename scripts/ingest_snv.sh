#!/bin/bash
set -e
#EUCANCAN SNV & INDEL vcf handling
KEEP=false

while [ $# -gt 0 ] ; do
  case $1 in
      -t | --truth) truth="$2";;
      -s | --snv) snv="$2";;
      -i | --indel) indel="$2";;
      -m | --snvindel) snvindel="$2";;
      -v | --sv) sv="$2";;
      -u | --truth_sv) truth_sv="$2";;
      -f | --fasta) FASTA="$2";;
      -d | --outdir) OUTPUT_DIR="$2";;
      -o | --outname) OUT_NAME="$2";;
      -n | --sname) SAMPLE_NAME="$2";;
      -k | --keep) KEEP=true;;
      -h | --help)  echo "Usage:"
          echo "bash ingest_snv.sh --truth truth_file.vcf"
          echo "                   --snv snv.vcf"
          echo "                   --indel indel.vcf"
          echo "                   --snvindel indels_and_vcfs.vcf"
          echo "                   --sv sv.vcf"
          echo "                   --truth_sv truth_sv.vcf"
          echo "                   --fasta ref_fasta.fa"
          echo "                   --outdir /OUTPUT_DIR/PATH"
          echo "                   --outname output file name"
          echo "                   --sname vcf sample name"
          echo "                   --keep (to keep intermediates files)"
          exit

  esac
  shift
done

echo " "
echo -e "[General Information]:\n"
echo "Truth File:" $truth
echo "SNV vcf file:" $snv
echo "INDEL vcf file:" $indel
echo "SNV + INDEL vcf file:" $snvindel
echo "SV vcf file:" $sv
echo "TRUTH SV file:" $truth_sv
echo "Reference fasta file:" $FASTA
echo "output path:" $OUTPUT_DIR
echo "output file Name:" $OUT_NAME
echo "vcf sample name:" $SAMPLE_NAME
echo "keep intermediate files ?:" $KEEP
echo " "

# Create output dir:

mkdir -p $OUTPUT_DIR/$OUT_NAME
OUTPUT_DIR=$OUTPUT_DIR/$OUT_NAME

# Load conda env:
#conda env create -n eucancan -f golden-datasets/scripts/environment_snv.yml

#source activate eucancan

# If SNV and INDEL in two files

if [[ ! -z "$snv" && ! -z "$indel" ]]; then
    echo "snv and indel not empty"
    if [[ $snv == *.vcf ]]; then
        bgzip -c $snv > $OUTPUT_DIR/$snv".gz"
        snv=$OUTPUT_DIR/$snv".gz"
    fi
    if [[ $indel == *.vcf ]]; then
        bgzip -c $indel > $OUTPUT_DIR/$indel".gz"
        indel=$OUTPUT_DIR/$indel".gz"
    fi

    echo -e "[Running Information]: indexing vcf files\n"
    bcftools index -f -o $snv".csi" $snv
    bcftools index -f -o $indel".csi" $indel

    echo -e "[Running Information]: concatenate vcf files\n"
    bcftools concat -a $snv $indel -O z -o $OUTPUT_DIR/"snv_indel_temp.vcf.gz"
    snvindel=$OUTPUT_DIR/"snv_indel_temp.vcf.gz"
fi

echo $snvindel
#bgzip
# index
# concat



# If SNV_INDEL:

## if snvindel vcf
if [[ $snvindel == *.vcf ]]; then
    cat $snvindel > $OUTPUT_DIR/snv_indel_temp.vcf
    gzip $OUTPUT_DIR/snv_indel_temp.vcf
    snvindel=$OUTPUT_DIR/snv_indel_temp.vcf.gz
fi

## if truth vcf
if [[ $truth == *.vcf ]]; then
    cat $truth > $OUTPUT_DIR/truth_temp.vcf
    gzip $OUTPUT_DIR/truth_temp.vcf
    truth=$OUTPUT_DIR/truth_temp.vcf.gz
fi

# Check if vcf is mono sample:

echo -e "[Running Information]: checking if vcf is multisample\n"

if [[ `bcftools query -l $snvindel |wc -l` -gt 1  && -z "$SAMPLE_NAME" ]]; then
    echo $snvindel "is a multisample"
    echo "sample name must be specified in -n parameter"
    exit
elif [[ `bcftools query -l $snvindel |wc -l` -gt 1  && ! -z "$SAMPLE_NAME" ]]; then
    echo $SAMPLE_NAME
    bcftools view -c1 -O z -s $SAMPLE_NAME -o $OUTPUT_DIR/snv_indel_temp.sample.vcf.gz $snvindel
    snvindel=$OUTPUT_DIR/snv_indel_temp.sample.vcf.gz
fi

## Replace the 'chr' with '' in the VCFs

echo -e "[Running Information]: replacing "" by "chr"\n"

zcat $snvindel | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | awk '{gsub(/contig=\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > $OUTPUT_DIR/snv_indel_temp.vcf
zcat $truth | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | awk '{gsub(/contig=\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > $OUTPUT_DIR/truth_temp.vcf

## Filtering PASS variants:

echo -e "[Running Information]: Filtering PASS variant\n"

grep "PASS\|#" $OUTPUT_DIR/snv_indel_temp.vcf > $OUTPUT_DIR/"snv_indel.pass.vcf"

snvindel=$OUTPUT_DIR/"snv_indel.pass.vcf"

# Sorting vcf files
echo -e "[Running Information]: sorting vcf files\n"

bcftools sort $snvindel -o $OUTPUT_DIR/"snv_indel.pass.sort.vcf.gz" -O z
bcftools sort $truth -o $OUTPUT_DIR/truth_temp.sort.vcf.gz -O z

snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.vcf.gz"
truth=$OUTPUT_DIR/"truth_temp.sort.vcf.gz"

# indexing sorted vcf files
echo -e "[Running Information]: indexing vcf files\n"

bcftools index -f -o $snvindel".csi" $snvindel
bcftools index -f -o $truth".csi" $truth

# Preparing for normalization
echo -e "[Running Information]: preparing for normalizing\n"
echo -e "[Running Information]: multimerge...]"

export HGREF=$FASTA

multimerge $snvindel -r $HGREF -o $OUTPUT_DIR/"snv_indel.pass.sort.prep.vcf.gz" --calls-only=1
multimerge $truth -r $HGREF -o $OUTPUT_DIR/truth_temp.sort.prep.vcf.gz --calls-only=1

snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.prep.vcf.gz"
truth=$OUTPUT_DIR/"truth_temp.sort.prep.vcf.gz"

# Normalizing vcfs:
echo -e "[Running Information]: normalizing test file\n"

pre.py $snvindel $OUTPUT_DIR/"snv_indel.pass.sort.prep.norm.vcf.gz" -L --decompose --somatic

echo -e "\n[Running Information]: normalizing truth file\n"

pre.py $truth $OUTPUT_DIR/truth_temp.sort.prep.norm.vcf.gz -L --decompose --somatic

snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.prep.norm.vcf.gz"
truth=$OUTPUT_DIR/"truth_temp.sort.prep.norm.vcf.gz"

# Running ingestion script:
echo -e "[Running Information]: Running ingestion_snv.py script \n"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

python $DIR/ingest_snv.py -samplename "SAMPLE" -o $OUTPUT_DIR/"snv_indel.pass.sort.prep.norm" $snvindel

python $DIR/ingest_snv.py -samplename "SAMPLE" -o $OUTPUT_DIR/"truth_temp.sort.prep.norm" $truth

snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.prep.norm.filtered.vcf"
truth=$OUTPUT_DIR/"truth_temp.sort.prep.norm.filtered.vcf"

# Running som.py:
echo -e "[Running Information]: Running som.py evaluation script\n"

som.py $truth $snvindel -o $OUTPUT_DIR/$OUT_NAME --verbose -N

echo -e "[Running Information]: script ended successfully\n"

# Running SV ingestion script:

#echo -e "[Running Information]: Running ingestion SV script\n"

# Cleaning
echo -e "[Running Information]: Cleaning files\n"

rm $OUTPUT_DIR/*{tbi,csi,json}

if [[ "$KEEP" = false ]];then
    rm $OUTPUT_DIR/*{vcf,vcf.gz}
fi


# Example commands
# SCRIPT_DIR=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/golden-datasets-fork-tom/golden-datasets/scripts
# OUT_DIR=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/test_ingestions_sh
# TRUTH=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/truth_dream/truth.snvs.synthetic.challenge.set1.chr.vcf.gz
# SNV=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/curie_colo829_snps.sample.vcf.gz
# INDEL=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/curie_colo829_indels.sample.vcf.gz
# SNV_INDEL=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/head_filtered_Mutect2_ERR2752450_vs_ERR2752449_snpEff.ann.vcf.gz
#
# ## with mixed indels & snvs
# bash $SCRIPT_DIR/temp_ingest_snv.sh --truth $TRUTH --snvindel $SNV_INDEL --outdir $OUT_DIR --outname test_multisample --sname COLO829T  --keep -f /data/annotations/pipelines/Human/hg19/genome/hg19.fa
#
# # with split snvs & indels
# bash $SCRIPT_DIR/temp_ingest_snv.sh --truth $TRUTH --snv $SNV --indel $INDEL --outdir $OUT_DIR --outname test_multisample --keep -f /data/annotations/pipelines/Human/hg19/genome/hg19.fa
