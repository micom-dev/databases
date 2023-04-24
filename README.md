<img src="https://github.com/micom-dev/micom/raw/main/docs/source/micom.png" width="60%">

[![Test and release](https://github.com/micom-dev/databases/actions/workflows/main.yml/badge.svg)](https://github.com/micom-dev/databases/actions/workflows/main.yml)
[![DOI for the current release](https://zenodo.org/badge/DOI/10.5281/zenodo.7739096.svg)](https://doi.org/10.5281/zenodo.7739096)
[![Follow MICOM on Mastodon](https://img.shields.io/mastodon/follow/109960852316221526?domain=https%3A%2F%2Fmstdn.science&style=social)](https://mstdn.science/@micom)

# MICOM model databases

Workflows and input data for the construction of model databases for use in MICOM.

See the releases for the Zenodo links to download databases for a specific version.

## Local setup

Install Qiime2 following the normal conda instructions. Let's assume your Qiime2 environments
is named `qiime2-2023.2`.

Add the additional requirements.

```bash
conda env update -n qiime2-2023.2 -f conda.yml

conda activate qiime2-2023.2
```

You will now be able to build the databases by running the Nextflow workflows in each
directory. If there is a download workflow this has to be run first. For instance, to
build AGORA2:

```bash
cd recipes/agora2
nextflow run download_agora2.nf -resume
nextflow run agora2.nf -resume
```

This will create the databases in `recipes/agora2/databases`.


## Organization

```
root
|---|release_notes
|---|recipes
      |----|[SOURCE]
           |---|[SOURCE].nf         # workflow to generate the DBs
               |data                # required input data
               |manifests           # lists of all models in a DB
               |databases           # databases that were built
```
