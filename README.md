# golden-datasets

Download, extract and format golden datasets to benchmark pipelines

## Usage


```shell
make COMMAND [OPTIONS]
```

> To process in parallel, the `-j N` option can be used where `N` is the is the maximum number of jobs that `make` will run in parallel.

| Command | Description |
| --- | --- |
|`download`| Download datasets in `data/datasets.tsv` |
|`test`| Download test data in `data/datasets.test.tsv` |
|`bam2fastq`| Convert downloaded `bam` files into `fastq` files |
|`clean`| Remove downloaded data |

## Download

To download datasets listed in `data/datasets.tsv` use the `make` command described below. 

```shell
make INSTALL_DIR=$HOME/golden-datasets
```

The installation directory is set from the environment with the variable `INSTALL_DIR` or directly through the command line as a positional arguments of the `make` command.