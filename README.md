# golden-datasets

Download, extract and format golden datasets to benchmark pipelines

## Usage

To download datasets listed in `data/datasets.tsv` use the `make` command described below. 

```shell
make INSTALL_DIR=/tmp
```

The installation directory is set from the environment with the variable `INSTALL_DIR` or directly through the command line as a positional arguments of the `make` command.

| Command | Description |
| --- | --- |
|`download`| Download datasets in `data/datasets.tsv` |
|`test`| Download test data in `data/datasets.test.tsv` |
|`clean`| Remove downloaded data |