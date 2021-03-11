repo_path=/Users/daphnevanbeek/IdeaProjects/golden-datasets/scripts/
source ${repo_path}/../venv/bin/activate


function run {
    repo=$1
    SV_VCF=$2
    samplename=$3
    out_sv_df=$4
    
    
    
    python ${repo}/ingest.py ${SV_VCF} -samplename ${samplename} -outputfile ${out_sv_df}
    
    #python ${repo}/ingest_snv.py
    
    
    

}


# Hartwig
run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/Hartwig/5.19/COLO829v003T.purple.sv.vcf.gz" "COLO829v003T" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/hartwig.csv"

# BSC
run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/BSC_SV_COLO.tsv" "none" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/bsc.csv"

# Curie
run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/Curie/COLO_curie_SV.vcf" "COLO829T" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/curie.csv"

# Truth file
run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/truth_file/COLO829_truth.tsv" "none" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/truth.csv"
