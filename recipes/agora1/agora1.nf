nextflow.enable.dsl=2

params.out = "${launchDir}"
params.manifest = "${params.out}/data/agora103_new_taxa.tsv"
params.agora_version = "103"
params.refseq = "216"
params.gtdb = "207"
params.database_version = "2"
params.models = "${params.out}/models"


workflow {
    levels = channel.from(["genus", "species"])

    DownloadGtdbTables()

    refseq_manifests = ConvertToRefseq(levels)
    gtdb_manifests = ConvertToGtdb(refseq_manifests, DownloadGtdbTables.out)

    refseq_manifests.concat(gtdb_manifests).view()

    BuildDb(refseq_manifests.concat(gtdb_manifests)) | GetManifest
}


// Process definitions

process ConvertToRefseq {
    cpus 1
    memory 4.GB
    time 1.h

    input:
    val(level)

    output:
    tuple val(level), val("refseq"), val("${params.refseq}"), path("agora1_refseq_${level}.tsv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from os import path

    manifest = pd.read_csv("${params.manifest}", sep="\\t")
    manifest["file"] = manifest.file.str.replace("agora/", "${params.models}/")
    ranks = pd.Index(["strain", "species", "genus" , "family","order", "class", "phylum", "kingdom"][::-1])
    manifest.to_csv("agora1_refseq_${level}.tsv", index=False, sep="\\t")
    """
}

process DownloadGtdbTables {
    cpus 1
    memory 8.GB
    time 8.h
    errorStrategy 'retry'
    maxRetries 3

    output:
    path("*.tsv")

    script:
    """
    wget --no-check-certificate -O gtdb_bac.tar.gz https://data.gtdb.ecogenomic.org/releases/release${params.gtdb}/${params.gtdb}.0/bac120_metadata_r${params.gtdb}.tar.gz
    wget --no-check-certificate -O gtdb_ar.tar.gz https://data.gtdb.ecogenomic.org/releases/release${params.gtdb}/${params.gtdb}.0/ar53_metadata_r${params.gtdb}.tar.gz
    tar -xf gtdb_bac.tar.gz
    tar -xf gtdb_ar.tar.gz
    """
}

process ConvertToGtdb {
    cpus 1
    memory 4.GB
    time 1.h

    input:
    tuple val(level), val(db), val(ver), path(manifest)
    path(tables)

    output:
    tuple val(level), val("gtdb"), val("${params.gtdb}"), path("agora1_gtdb_${level}.tsv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    gtdb_rank_names = pd.Index(["domain", "phylum", "class", "order", "family", "genus", "species"])

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
    merged = pd.merge(agora, meta, left_on="ncbiid", right_on="ncbi_taxid")
    rank_mappings = merged.groupby("file")["${level}"].nunique()
    valid = rank_mappings.index[rank_mappings == 1]
    merged = merged[merged.file.isin(valid)].drop_duplicates(subset=["file", "${level}"])

    print(f"rank: ${level} - matched: {merged.shape[0]}/{agora.shape[0]}")

    merged.to_csv("agora1_gtdb_${level}.tsv", index=False, sep="\\t")
    """

}

process BuildDb {
    publishDir "${params.out}/databases", mode: "copy", overwrite: true
    cpus 12
    memory 16.GB
    time 8.h

    input:
    tuple val(level), val(db), val(ver), path(manifest)

    output:
    tuple val(level), val(db), val(ver), path("*.qza")

    script:
    """
    qiime micom db --m-meta-file ${manifest} \
        --p-rank ${level} \
        --p-threads ${task.cpus} \
        --verbose \
        --o-metabolic-models agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.qza
    """
}

process GetManifest {
    cpus 1
    memory 4.GB
    time 1.h
    publishDir "${params.out}/manifests", mode: "copy", overwrite: true

    input:
    tuple val(level), val(db), val(ver), path(arti)

    output:
    path("agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.tsv")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import qiime2 as q2

    db = q2.Artifact.load("${arti}")
    manifest = db.manifest.view(pd.DataFrame)
    manifest.to_csv(
        "agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.tsv",
        sep="\\t",
        index=False
    )
    """
}
