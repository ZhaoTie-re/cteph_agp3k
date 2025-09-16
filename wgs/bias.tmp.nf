params.vcfFilePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/06.variant_filter'
params.plinkFilePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter'
params.infoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'

Channel
    .from((1..22).collect { "chr${it}" })
    .set { chr_ch }

process BiasCalculation {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outDir}/21.geno_bias_calc", mode: 'symlink'

    input:
    val(chr) from chr_ch
    val(infoPath) from params.infoPath
    val(vcfFilePath) from params.vcfFilePath
    val(plinkFilePath) from params.plinkFilePath

    output:
    file("*.log")
    tuple val(chr), file("${chr}.coverage_transitions.tsv"), file("${chr}.stat_summary.tsv") into bias_result_ch

    script:
    bed_prefix = "${plinkFilePath}/cteph_agp3k.lowfreq_common"
    vcf = "${vcfFilePath}/${chr}.pass.mac1.vfilter.vcf.gz"
    info = "${infoPath}/cteph_agp3k_jhrpv4.xlsx"
    """
    source activate compute_env
    python ${params.scriptDir}/geno_miss_bias_chr_main.py \
        --plink-prefix ${bed_prefix} \
        --vcf ${vcf} \
        --info-xls ${info} \
        --chr ${chr} \
        --n-chunks 12 \
        --max-parallel 8 \
        --plink-threads 8 \
        --bcftools-threads 8 \
        --keep-temp
    """
}

bias_result_ch
    .map { the_chr, cov, sum -> [
        chr: the_chr.toString(),
        coverage_transitions: cov.toString(),
        stat_summary: sum.toString()
    ] }
    .collect()
    .map { entries -> groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([bias_results: entries])) }
    .set { bias_results_json_text_ch }

process BiasResultsToJson {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'bias_json'

    publishDir "${params.outDir}/21.geno_bias_calc", mode: 'copy'

    input:
    val(json_text) from bias_results_json_text_ch

    output:
    file('bias_results.json') into bias_results_json_ch

    script:
    """
    echo '${json_text}' > bias_results.json
    """
}


process MergeBiasResultsJson {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag 'merge_bias_json'

    publishDir "${params.outDir}/22.geno_bias_corr", mode: 'symlink'

    input:
    file(json_file) from bias_results_json_ch
    val(scriptDir) from params.scriptDir

    output:
    file('all.coverage_transitions.tsv')
    file('all.stat_summary.tsv')
    file('bias_results.merged.json') into merged_bias_results_json_ch

    script:
    """
    source activate compute_env
    export JSON_PATH="${json_file}"
    export SCRIPT_DIR="${scriptDir}"
    python - <<'PY'
import os, sys, importlib
# 将脚本目录加入 sys.path（来自环境变量，避免 Nextflow 字符串插值）
script_path = os.path.abspath(os.environ['SCRIPT_DIR'])
if script_path not in sys.path:
    sys.path.append(script_path)

import geno_miss_bias_tools
importlib.reload(geno_miss_bias_tools)
from geno_miss_bias_tools import merge_bias_results_json

json_path = os.environ['JSON_PATH']
# 在当前工作目录输出合并结果（函数内部已固定输出到 CWD）
ret = merge_bias_results_json(json_path, out_prefix="all")
print("[INFO] merged_coverage_transitions=", ret.get("merged_coverage_transitions"))
print("[INFO] merged_stat_summary=", ret.get("merged_stat_summary"))
print("[INFO] merged_json=", ret.get("merged_json"))
PY
    """
}


