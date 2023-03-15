<img src="https://github.com/micom-dev/q2-micom/raw/main/docs/assets/logo.png" width="75%">

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
