nextflow.enable.dsl=2

params.repo_version = "260d0f1"
params.url = "https://github.com/cdanielmachado/embl_gems/archive/refs/heads/master.zip"
params.refseq = "216"
params.gtdb = "207"
params.database_version = "1"

process download_models {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3

    output:
    tuple path("model_list.tsv"), path("models")

    """
    wget ${params.url}
    unzip master.zip
    mv embl_gems-master/models .
    mv embl_gems-master/model_list.tsv .
    rm -rf embl_gems-master
    """
}

process assign_taxonomy {
    cpus 1

    input:
    tuple path(manifest), path(models)

    output:
    tuple path(manifest), path("lineage.txt")

    """
    awk -F"\\t" 'NR>1 {print \$3}' ${manifest} | \
        taxonkit reformat -I 1 -f "{k};{p};{c};{o};{f};{g};{s};{t}" -r NA > lineage.txt
    """
}

process convert_to_refseq {
    cpus 1
    publishDir "${baseDir}/manifests", mode: "copy", overwrite: true

    input:
    tuple val(level), path(manifest), path(lineage)

    output:
    tuple val(level), val("refseq"), val("${params.refseq}"), path("carveme_refseq_${level}.tsv")
    """
    #!/usr/bin/env python

    import pandas as pd
    from os import path

    manifest = pd.read_csv("${manifest}", sep="\\t")
    manifest.rename(columns={"file_path": "file", "assembly_accession": "id"}, inplace=True)
    lineage = pd.read_csv("${lineage}", sep="\\t")
    tax = lineage.iloc[:, 1].str.split(";", expand=True)
    tax.columns = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
    merged = pd.concat([manifest, tax], axis=1)
    merged = merged[merged["${level}"] != "NA"]
    merged.to_csv("carveme_refseq_${level}.tsv", index=False, sep="\\t")
    """
}

process download_gtdb_tables {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3

    output:
    path("*.tsv")

    """
    wget --no-check-certificate -O gtdb_bac.tar.gz https://data.gtdb.ecogenomic.org/releases/release${params.gtdb}/${params.gtdb}.0/bac120_metadata_r${params.gtdb}.tar.gz
    wget --no-check-certificate -O gtdb_ar.tar.gz https://data.gtdb.ecogenomic.org/releases/release${params.gtdb}/${params.gtdb}.0/ar53_metadata_r${params.gtdb}.tar.gz
    tar -xf gtdb_bac.tar.gz
    tar -xf gtdb_ar.tar.gz
    """
}

process convert_to_gtdb {
    cpus 1
    publishDir "${baseDir}/manifests", mode: "copy", overwrite: true

    input:
    tuple val(level), val(db), val(ver), path(manifest)
    path(tables)

    output:
    tuple val(level), val("gtdb"), val("${params.gtdb}"), path("carveme_gtdb_${level}.tsv")

    """
    #!/usr/bin/env python

    import pandas as pd

    gtdb_rank_names = pd.Index(["kingdom", "phylum", "class", "order", "family", "genus", "species"])

    files = "${tables}".split()
    meta = pd.concat(
        pd.read_csv(fi, sep="\\t", usecols=["gtdb_taxonomy", "ncbi_taxid"])
        for fi in files
    )

    rank_idx = gtdb_rank_names.get_loc("${level}")
    tax = meta.gtdb_taxonomy.str.split(";", expand=True)
    tax.columns = gtdb_rank_names
    meta = pd.concat([meta, tax], 1).drop(["gtdb_taxonomy"], axis=1).drop_duplicates()

    manifest = pd.read_csv("${manifest}", sep="\\t").drop(
        ["strain", "species", "genus" , "family","order", "class", "phylum", "kingdom"],
        axis=1, errors="ignore"
    )
    merged = pd.merge(manifest, meta, left_on="taxid", right_on="ncbi_taxid")
    rank_mappings = merged.groupby("file")["${level}"].nunique()
    valid = rank_mappings.index[rank_mappings == 1]
    merged = merged[merged.file.isin(valid)].drop_duplicates(subset=["file", "${level}"])

    print(f"rank: ${level} - matched: {merged.shape[0]}/{manifest.shape[0]}")

    merged.to_csv("carveme_gtdb_${level}.tsv", index=False, sep="\\t")
    """

}

process build_db {
    publishDir "${baseDir}/databases", mode: "copy", overwrite: true
    cpus 12

    input:
    tuple val(level), val(db), val(ver), path(manifest), path(models)

    output:
    path("*.qza")

    """
    qiime micom db --m-meta-file ${manifest} \
        --p-rank ${level} \
        --p-threads ${task.cpus} \
        --verbose \
        --o-metabolic-models carveme${params.repo_version}_${db}${ver}_${level}_${params.database_version}.qza
    """
}

workflow {
    levels = Channel.from(["genus", "species"])

    download_models | assign_taxonomy

    download_gtdb_tables()

    refseq_manifests = convert_to_refseq(levels.combine(assign_taxonomy.out))
    gtdb_manifests = convert_to_gtdb(refseq_manifests, download_gtdb_tables.out)

    refseq_manifests.concat(gtdb_manifests).view()

    build_db(refseq_manifests.concat(gtdb_manifests).combine(download_models.out.map{it -> it[1]}))
}
