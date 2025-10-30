import sys, importlib
sys.path.insert(0, "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/scripts")

import array_check_prepare_tools
importlib.reload(array_check_prepare_tools)

import matplotlib
matplotlib.use('Agg')

from array_check_prepare_tools import plot_wgs_array_beta_comparison

# 使用和notebook中相同的参数
plot = plot_wgs_array_beta_comparison(
    merged_sumstat_file='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute/cteph_agp3k.array_wgs_compration.merged_sumstat.tsv',
    output_prefix="cteph_agp3k.array_wgs_compration",
    save_pdf=False,
)

print("绘图完成！")
print(plot.get('summary_text', ''))
