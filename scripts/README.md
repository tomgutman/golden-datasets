# Scripts

Folder with development scripts.

To ensure execution use: 

```shell
conda env create --name ingestion python=3.6 -f environment.yml
```

```shell
conda activate ingestion
```

Add any new requirement for your code to [environment.yml](https://github.com/EUCANCan/golden-datasets/blob/devel/scripts/environment.yml) (if it is compatible)


## 1) ingestion

Creates dataframes from VCF/TSV files. This step is required both for test files and truth files.
Flag -samplename is required for VCF files and must contain the Tumor sample name

```shell
python ingest.py VCF_FILE -outputfile DATAFRAME_FILE
```


## 2) Benchmarking

Computes metrics for SV calls. Optional: `-metrics` flag to print result to a file

```shell
python compare_node_to_truth.py DATAFRAME_TEST DATAFRAME_TRUTH -metrics OUTPUT_METRICS_FILE
```
