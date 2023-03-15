nextflow.enable.dsl=2

params.manifest = "${baseDir}/data/agora201.csv"
params.agora_version = "201"
params.refseq = "216"
params.gtdb = "207"
params.database_version = "1"
params.models = "${baseDir}/models"

process convert_to_refseq {
    cpus 1
    publishDir "${baseDir}/manifests", mode: "copy", overwrite: true

    input:
    val(level)

    output:
    tuple val(level), val("refseq"), val("${params.refseq}"), path("agora2_refseq_${level}.tsv")

    """
    #!/usr/bin/env python

    import pandas as pd
    from os import path

    manifest = pd.read_csv("${params.manifest}")
    del manifest["sbml_model"]
    manifest["file"] = [path.join("${params.models}", id + ".xml.gz") for id in manifest["MicrobeID"]]
    ranks = pd.Index(["Strain", "Species", "Genus" , "Family","Order", "Class", "Phylum", "Kingdom"][::-1])
    manifest.rename(columns=dict(zip(ranks, ranks.str.lower())), inplace=True)
    manifest.rename(columns={"MicrobeID": "id"}, inplace=True)
    print(manifest)
    manifest.to_csv("agora2_refseq_${level}.tsv", index=False, sep="\\t")
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
    tuple val(level), val("gtdb"), val("${params.gtdb}"), path("agora2_gtdb_${level}.tsv")

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

    agora = pd.read_csv("${manifest}", sep="\\t").drop(
        ["strain", "species", "genus" , "family","order", "class", "phylum", "kingdom"],
        axis=1, errors="ignore"
    )
    merged = pd.merge(agora, meta, left_on="NCBI Taxonomy ID", right_on="ncbi_taxid")
    rank_mappings = merged.groupby("file")["${level}"].nunique()
    valid = rank_mappings.index[rank_mappings == 1]
    merged = merged[merged.file.isin(valid)].drop_duplicates(subset=["file", "${level}"])

    print(f"rank: ${level} - matched: {merged.shape[0]}/{agora.shape[0]}")

    merged.to_csv("agora2_gtdb_${level}.tsv", index=False, sep="\\t")
    """

}

process build_db {
    publishDir "${baseDir}/databases", mode: "copy", overwrite: true
    cpus 12

    input:
    tuple val(level), val(db), val(ver), path(manifest)

    output:
    path("*.qza")

    """
    qiime micom db --m-meta-file ${manifest} \
        --p-rank ${level} \
        --p-threads ${task.cpus} \
        --verbose \
        --o-metabolic-models agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.qza
    """
}

workflow {
    levels = Channel.from(["genus", "species"])

    download_gtdb_tables()

    refseq_manifests = convert_to_refseq(levels)
    gtdb_manifests = convert_to_gtdb(refseq_manifests, download_gtdb_tables.out)

    refseq_manifests.concat(gtdb_manifests).view()

    build_db(refseq_manifests.concat(gtdb_manifests))
}
