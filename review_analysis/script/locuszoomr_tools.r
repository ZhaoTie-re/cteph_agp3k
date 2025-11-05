# ==============================================================================
# 区域关联图绘制工具 (locuszoomr)
# ==============================================================================
# 描述: 生成 LocusZoom 风格的区域关联图，包含 LD 信息和重组率
# 作者: ZHAO TIE
# 日期: 2025-11-05
# 
# 环境要求:
#   - Conda 环境: r_work
#   - 激活命令: conda activate r_work
# ==============================================================================

# %%
# 加载必需的R包
# ------------------------------------------------------------------------------
library(EnsDb.Hsapiens.v86)  # 人类基因组注释 (Ensembl v86)
library(locuszoomr)           # 区域关联图绘制
library(rtracklayer)          # 基因组数据格式导入

# %%
# 配置参数
# ------------------------------------------------------------------------------
# 输入文件路径
sum_path <- '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/ld_summary/chr3_154069965_A_G.merged_ld_sum_stat.tsv'

# LD 数据源选择
# 可选项: "ld_r2_with_lead" (AGP3K), "tommo_r2_with_lead" (ToMMo), "eas_r2_with_lead" (1000G-EAS)
ld_column <- "tommo_r2_with_lead"

# 重组率文件 (hg38)
recomb_file <- "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/info/recomb1000GAvg.bw"

# 输出文件夹
output_dir <- "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/result_vis"

# 目标基因
target_gene <- "ARHGEF26-AS1"

# 目标 SNP
target_snp <- c("chr3:154069965:A:G")

# 图形参数设置
plot_width <- 8         # 图形宽度（英寸）
plot_height <- 5        # 图形高度（英寸）
png_res <- 400          # PNG 分辨率（DPI）
cex_main <- 1.2         # 主标题字体大小
cex_lab <- 1.3          # 轴标签字体大小
cex_axis <- 1.1         # 轴刻度字体大小
cex_gene <- 1.0         # 基因标签字体大小

# %%
# 加载和准备数据
# ------------------------------------------------------------------------------
sum_df <- read.table(sum_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# head(sum_df)

# %%
# 创建 locus 对象
# ------------------------------------------------------------------------------
# 提取基因组区域信息
chrom_num <- unique(sum_df$CHR)[1]
start_pos <- min(sum_df$POS)
end_pos <- max(sum_df$POS)

# 使用 GWAS 汇总统计初始化 locus 对象
loc <- locus(data = sum_df,
             chrom = "CHR",           # 染色体列名
             pos = "POS",             # 位置列名
             p = "P",                 # P值列名
             labs = "SNPID",          # SNP ID列名
             LD = ld_column,          # LD r²列名
             seqname = chrom_num,     # 染色体编号
             xrange = c(start_pos, end_pos),  # 基因组范围
             ens_db = "EnsDb.Hsapiens.v86")   # 基因注释数据库

summary(loc)

# %%
# 添加重组率信息
# ------------------------------------------------------------------------------
recomb.hg38 <- import.bw(recomb_file)
loc <- link_recomb(loc, recomb = recomb.hg38)

# %%
# 生成区域关联图
# ------------------------------------------------------------------------------
# 确定 LD 数据来源标签
ld_source <- switch(ld_column,
                    "ld_r2_with_lead" = "AGP3K",
                    "tommo_r2_with_lead" = "ToMMo",
                    "eas_r2_with_lead" = "1000G-EAS",
                    "Unknown")

plot_title <- paste0("LD R² source: ", ld_source)

# 自定义面板，用于添加标题
panel_title <- quote({
  title(main = plot_title, line = 0.2, cex.main = cex_main, font.main = 2)
})

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 生成输出文件名（基于输入文件名和LD来源）
base_name <- sub(".merged_ld_sum_stat.tsv$", "", basename(sum_path))
output_prefix <- file.path(output_dir, paste0(base_name, "_", ld_source))

# 定义统一的绘图函数
plot_locus <- function() {
  par(cex.lab = cex_lab, cex.axis = cex_axis, mar = c(5,5,2,2))
  locus_plot(loc, 
            #  filter_gene_biotype = c("protein_coding"),
             filter_gene_biotype = c("protein_coding", "processed_transcript"),
             highlight = target_gene,
             labels = target_snp,
            #  label_x = c(4, -5),
             cex = cex_gene,
            #  gene_col = 'grey', exon_col = 'orange', exon_border = 'darkgrey',
             panel.last = panel_title)
}

# 保存为高清 PDF
pdf(paste0(output_prefix, ".pdf"), width = plot_width, height = plot_height)
plot_locus()
dev.off()

# 保存为高清 PNG
png(paste0(output_prefix, ".png"), 
    width = plot_width * png_res, 
    height = plot_height * png_res, 
    res = png_res)
plot_locus()
dev.off()

# 在屏幕显示图形
plot_locus()

# %%
