import pandas as pd
import numpy as np
import gwaslab as gl
import logging
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt

# 设置全局日志级别变量
logging.getLogger("matplotlib").setLevel(logging.ERROR)
plt.style.use('default')
plt.ioff()  # 关闭交互模式

# 设置更保守的图形参数
plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['agg.path.chunksize'] = 10000

sig_level = 5e-8

array_add_path = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute/04.array_assoc_result/cteph_agp3k.array.sex.10pc.additive.PHENO1.glm.logistic"

print("="*80)
print("Loading GWAS summary statistics...")
print("="*80)

# 加载数据
sumstats_add = gl.Sumstats(
    array_add_path, 
    fmt="plink2", 
    build='38', 
    ea='A1', 
    nea='OMITTED', 
    OR_95L='L95', 
    OR_95U='U95'
)

print(f"\n数据加载完成！")
print(f"总变异数: {len(sumstats_add.data)}")
print(f"数据列: {sumstats_add.data.columns.tolist()}")

# 检查 NaN 值
print(f"\nNaN 值统计:")
print(sumstats_add.data.isnull().sum())

print("\n" + "="*80)
print("Generating Manhattan and QQ plots...")
print("="*80)

# 绘制 Manhattan 和 QQ 图
try:
    fig = sumstats_add.plot_mqq(
        mode='mqq',
        dpi=300,
        stratified=False,
        verbose=True,
        save='array_mqq_plot.png',
        save_args={'dpi': 300, 'bbox_inches': 'tight'}
    )
    print("\n✓ Manhattan and QQ plots saved successfully!")
except Exception as e:
    print(f"\n✗ Error plotting mqq: {e}")
    print("\nTrying to plot separately...")
    
    # 尝试分别绘制
    try:
        print("  - Plotting Manhattan...")
        fig_m = sumstats_add.plot_mqq(
            mode='m',
            dpi=300,
            stratified=False,
            verbose=True,
            save='array_manhattan.png',
            save_args={'dpi': 300, 'bbox_inches': 'tight'}
        )
        plt.close('all')
        print("  ✓ Manhattan plot saved!")
    except Exception as e2:
        print(f"  ✗ Manhattan plot failed: {e2}")
    
    try:
        print("  - Plotting QQ...")
        fig_qq = sumstats_add.plot_mqq(
            mode='qq',
            dpi=300,
            stratified=False,
            verbose=True,
            save='array_qq.png',
            save_args={'dpi': 300, 'bbox_inches': 'tight'}
        )
        plt.close('all')
        print("  ✓ QQ plot saved!")
    except Exception as e3:
        print(f"  ✗ QQ plot failed: {e3}")

finally:
    plt.close('all')

print("\n" + "="*80)
print("Extracting lead SNPs...")
print("="*80)

# 获取 lead SNPs
try:
    lead_snps = sumstats_add.get_lead(sig_level=5e-8, verbose=True)
    print(f"\n✓ Found {len(lead_snps)} lead SNPs at p < 5e-8")
    
    # 保存 lead SNPs
    lead_snps.to_csv('array_lead_snps.csv', index=False)
    print("✓ Lead SNPs saved to: array_lead_snps.csv")
    
    # 显示前10个
    if len(lead_snps) > 0:
        print("\nTop 10 lead SNPs:")
        print(lead_snps.head(10).to_string())
except Exception as e:
    print(f"✗ Error getting lead SNPs: {e}")

print("\n" + "="*80)
print("Checking specific SNP: chr17:13528059:G:A")
print("="*80)

# 检查特定 SNP
try:
    specific_snp = sumstats_add.data[sumstats_add.data['SNPID'] == 'chr17:13528059:G:A']
    if len(specific_snp) > 0:
        print("\n✓ SNP found:")
        print(specific_snp.to_string())
        specific_snp.to_csv('array_specific_snp_chr17_13528059.csv', index=False)
        print("✓ Saved to: array_specific_snp_chr17_13528059.csv")
    else:
        print("✗ SNP not found in dataset")
except Exception as e:
    print(f"✗ Error checking specific SNP: {e}")

print("\n" + "="*80)
print("Analysis completed!")
print("="*80)

# 保存完整的summary statistics
print("\nSaving full summary statistics...")
sumstats_add.data.to_csv('array_full_sumstats.csv.gz', index=False, compression='gzip')
print("✓ Full summary statistics saved to: array_full_sumstats.csv.gz")

print("\nAll outputs saved in: /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute/")
