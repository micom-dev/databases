nextflow.enable.dsl=2

params.manifest = "data/agora201.csv"
params.url = "https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/individual_reconstructions/"

workflow {
    channel.fromPath("${baseDir}/${params.manifest}")
        .splitCsv(header: true, quote: '"')
        .map{row -> row.MicrobeID}
        .set{ids}

    ids | download | convert_to_sbml
}


// Process definitions

process download {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(id)

    output:
    tuple val(id), path("${id}.mat")

    script:
    """
    wget --no-check-certificate -O ${id}.mat ${params.url}/${id}.mat
    """
}

process convert_to_sbml {
    publishDir "${baseDir}/models/", mode: 'copy', overwrite: false
    cpus 1

    input:
    tuple val(id), path(mat)

    output:
    tuple val(id), path("${id}.xml.gz")

    script:
    """
    #!/usr/bin/env python

    from cobra.io import load_matlab_model, write_sbml_model

    mod = load_matlab_model("${mat}")
    assert mod.slim_optimize() > 1e-3
    write_sbml_model(mod, "${id}.xml.gz")
    """
}
