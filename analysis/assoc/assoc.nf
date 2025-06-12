params.gtPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc'

process pca_prepare {
    

    input:
    tuple val(sampleId), path(gtFile)

    output:
    tuple val(sampleId), path("*.pca")

    script:
    """
    Rscript ${params.scriptDir}/pca_prepare.R ${gtFile} ${sampleId}
    """
}