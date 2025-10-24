#!/usr/bin/env Rscript
# 作者: ZHAO TIE
# =================================================================
# SuSiE 精细定位分析命令行工具 (SuSiE Fine-mapping Analysis CLI Tool)
# =================================================================
# 
# 描述 (Description): 使用SuSiE (Sum of Single Effects) 方法进行精细定位分析
# 
# SuSiE是一种用于基因组关联研究(GWAS)精细定位的统计方法，能够：
# 1. 识别独立的因果变异位点
# 2. 计算每个变异位点的后验包含概率(PIP)
# 3. 构建可信集合(Credible Sets)
# 4. 评估信号的纯度(Purity)
#
# 输入文件要求：
# - {lead_id}.sum_stat.tsv: 汇总统计数据，包含BETA, SE等列
# - {lead_id}.r.ld: LD相关系数矩阵(无表头的数值矩阵)
# - {lead_id}.ld_r.tsv: 带变异ID的LD相关系数矩阵
# - {lead_id}.ld_r2.tsv: 带变异ID的LD决定系数矩阵
# =================================================================

# 加载必需的R包
suppressPackageStartupMessages({
  library(susieR)    # SuSiE精细定位算法
  library(jsonlite)  # JSON数据格式处理
  library(optparse)  # 命令行参数解析
})

# =================================================================
# 命令行参数定义 (Command Line Options Definition)
# =================================================================
option_list <- list(
  # 样本量 (必需参数)
  make_option(c("-s", "--sample_size"), 
              type="integer", 
              default=NULL,
              help="分析所用的样本量 [必需] (Sample size for analysis [required])", 
              metavar="整数"),
  
  # 主导变异ID (必需参数) 
  make_option(c("-l", "--lead_id"), 
              type="character", 
              default=NULL,
              help="主导变异ID，格式为 chr:pos:ref:alt [必需] (Lead variant ID [required])", 
              metavar="字符串"),
  
  # 输入文件基础路径 (必需参数)
  make_option(c("-b", "--base_path"), 
              type="character", 
              default=NULL,
              help="包含输入文件的基础路径 [必需] (Base path for input files [required])", 
              metavar="路径"),
  
  # 输出目录 (必需参数)
  make_option(c("-o", "--out_dir"), 
              type="character", 
              default=NULL,
              help="结果输出目录 [必需] (Output directory [required])", 
              metavar="路径"),
  
  # PIP阈值 (可选参数，默认0.001)
  make_option(c("-p", "--pip_threshold"), 
              type="double", 
              default=0.001,
              help="变异位点过滤的PIP阈值 [默认: 0.001] (PIP threshold for filtering [default: 0.001])", 
              metavar="小数"),
  
  # SuSiE最大组件数 (可选参数，默认10)
  make_option(c("-L", "--max_components"), 
              type="integer", 
              default=10,
              help="SuSiE算法的最大组件数 [默认: 10] (Max components for SuSiE [default: 10])", 
              metavar="整数"),
  
  # 可信集合覆盖度 (可选参数，默认0.95)
  make_option(c("-c", "--coverage"), 
              type="double", 
              default=0.95,
              help="可信集合的覆盖度 [默认: 0.95] (Credible set coverage [default: 0.95])", 
              metavar="小数"),
  
  # 随机种子 (可选参数，默认123)
  make_option(c("--seed"), 
              type="integer", 
              default=123,
              help="随机种子，确保结果可重现 [默认: 123] (Random seed [default: 123])", 
              metavar="整数"),
  
  # 详细输出模式 (可选参数)
  make_option(c("-v", "--verbose"), 
              action="store_true", 
              default=FALSE,
              help="启用详细输出模式 (Enable verbose output)"),
  
  # 帮助信息
  make_option(c("-h", "--help"), 
              action="store_true", 
              default=FALSE,
              help="显示帮助信息并退出 (Show help message and exit)")
)

# =================================================================
# 命令行参数解析器配置 (Command Line Parser Configuration)
# =================================================================
opt_parser <- OptionParser(
  option_list=option_list,
  usage="用法: %prog [选项] (Usage: %prog [options])",
  add_help_option=FALSE,
  description=paste(
    "=================================================================",
    "SuSiE 精细定位分析工具 (SuSiE Fine-mapping Analysis Tool)",
    "=================================================================",
    "",
    "本工具使用SuSiE方法进行GWAS精细定位分析。",
    "This tool performs fine-mapping analysis using the SuSiE method.",
    "",
    "需要在base_path目录中准备以下输入文件:",
    "Required input files in base_path:",
    "  - {lead_id}.sum_stat.tsv: 汇总统计数据 (Summary statistics)",
    "  - {lead_id}.r.ld: LD相关系数矩阵 (LD correlation matrix)",
    "  - {lead_id}.ld_r.tsv: 带ID的LD相关系数矩阵 (LD matrix with IDs)",
    "  - {lead_id}.ld_r2.tsv: 带ID的LD决定系数矩阵 (LD r-squared matrix with IDs)",
    "",
    "使用示例 (Example):",
    "  Rscript susie_cli.r -s 2644 -l 'chr17:13528059:G:A' \\",
    "    -b /path/to/input -o /path/to/output -p 0.001",
    "",
    "输出文件说明 (Output files):",
    "  - {lead_id_formatted}.susie_results.json: SuSiE分析结果",
    "  - {lead_id_formatted}.susie_analysis.log: 详细日志文件",
    "=================================================================",
    sep="\n"
  )
)

# 解析命令行参数
opt <- parse_args(opt_parser)

# 检查帮助标志
if (!is.null(opt$help) && isTRUE(opt$help)) {
  print_help(opt_parser)
  quit(status=0)
}

# =================================================================
# 参数验证 (Parameter Validation)
# =================================================================
# 验证必需参数是否提供
required_args <- c("sample_size", "lead_id", "base_path", "out_dir")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

if (length(missing_args) > 0) {
  cat("错误: 缺少必需参数 (Error: Missing required arguments):", 
      paste(missing_args, collapse=", "), "\n\n")
  print_help(opt_parser)
  quit(status=1)
}

# =================================================================
# 参数提取和初始化 (Parameter Extraction and Initialization)
# =================================================================
# 从命令行参数中提取分析参数
sample_size <- opt$sample_size      # 样本量
lead_id <- opt$lead_id              # 主导变异ID
base_path <- opt$base_path          # 输入文件基础路径
out_dir <- opt$out_dir              # 输出目录
pip_threshold <- opt$pip_threshold  # PIP阈值
max_components <- opt$max_components # SuSiE最大组件数
coverage <- opt$coverage            # 可信集合覆盖度
seed <- opt$seed                    # 随机种子
verbose <- opt$verbose              # 详细输出模式

# =================================================================
# 日志系统设置 (Logging System Setup)
# =================================================================
# 生成日志文件路径（将冒号替换为下划线以避免文件名问题）
log_file <- file.path(out_dir, paste0(gsub(":", "_", lead_id), ".susie_analysis.log"))

# 创建输出目录（如果不存在）
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 初始化日志文件连接
log_conn <- file(log_file, "w")

# 日志写入函数
write_log <- function(message, level = "信息") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- sprintf("[%s] %s: %s", timestamp, level, message)
  writeLines(log_message, log_conn)
  flush(log_conn)  # 立即写入文件
  if (verbose || level == "错误") {
    cat(log_message, "\n")  # 如果是详细模式或错误，同时输出到控制台
  }
}

# =================================================================
# 主要分析流程 - 错误处理包装器 (Main Analysis with Error Handling)
# =================================================================
tryCatch({
  
  # 记录分析开始
  write_log("=== SuSiE 精细定位分析开始 ===")
  write_log(sprintf("主导变异: %s", lead_id))
  write_log(sprintf("样本量: %d", sample_size))
  write_log(sprintf("输入路径: %s", base_path))
  write_log(sprintf("输出目录: %s", out_dir))
  write_log(sprintf("PIP阈值: %.3f", pip_threshold))
  write_log(sprintf("最大组件数: %d", max_components))
  write_log(sprintf("覆盖度: %.2f", coverage))
  write_log(sprintf("随机种子: %d", seed))
  
  # =================================================================
  # 输入文件路径定义 (Input File Path Definition)
  # =================================================================
  # 根据主导变异ID构建输入文件路径
  sumstat_path <- file.path(base_path, paste0(lead_id, ".sum_stat.tsv"))
  ld_r_path <- file.path(base_path, paste0(lead_id, ".r.ld"))
  ld_r_ids_tsv_path <- file.path(base_path, paste0(lead_id, ".ld_r.tsv"))
  ld_r2_ids_tsv_path <- file.path(base_path, paste0(lead_id, ".ld_r2.tsv"))
  
  write_log("=== 输入文件检查 ===")
  write_log(sprintf("汇总统计文件: %s", sumstat_path))
  write_log(sprintf("LD矩阵文件: %s", ld_r_path))
  write_log(sprintf("带ID的LD相关系数矩阵: %s", ld_r_ids_tsv_path))
  write_log(sprintf("带ID的LD决定系数矩阵: %s", ld_r2_ids_tsv_path))
  
  # 检查输入文件是否存在
  input_files <- c(sumstat_path, ld_r_path, ld_r_ids_tsv_path, ld_r2_ids_tsv_path)
  missing_files <- input_files[!file.exists(input_files)]
  
  if (length(missing_files) > 0) {
    write_log(sprintf("缺少输入文件: %s", paste(missing_files, collapse=", ")), "错误")
    quit(status=1)
  }
  
  write_log("所有输入文件均存在")
  
  # =================================================================
  # 数据加载 (Data Loading)
  # =================================================================
  write_log("=== 数据加载 ===")
  
  # 加载汇总统计数据
  sumstats <- read.csv(sumstat_path, header = TRUE, sep = "\t")
  write_log(sprintf("汇总统计数据已加载: %d 个变异位点", nrow(sumstats)))
  
  # 加载LD矩阵（无表头的数值矩阵，用于SuSiE计算）
  ld_mat <- as.matrix(read.csv(ld_r_path, header = FALSE, sep = "\t"))
  write_log(sprintf("LD矩阵已加载: %d x %d", nrow(ld_mat), ncol(ld_mat)))
  
  # 加载带变异ID的LD矩阵（用于结果解释和LD值提取）
  ld_r_with_ids <- read.csv(ld_r_ids_tsv_path, header = TRUE, sep = "\t", row.names = 1)
  ld_r2_with_ids <- read.csv(ld_r2_ids_tsv_path, header = TRUE, sep = "\t", row.names = 1)
  write_log(sprintf("带ID的LD矩阵已加载: %d x %d", nrow(ld_r_with_ids), ncol(ld_r_with_ids)))
  
  # =================================================================
  # 数据验证 (Data Validation)
  # =================================================================
  write_log("=== 数据质量检查 ===")
  
  # 检查LD矩阵对称性（SuSiE要求LD矩阵为对称矩阵）
  max_diff <- max(abs(ld_mat - t(ld_mat)), na.rm = TRUE)
  is_symmetric <- isSymmetric(ld_mat)
  write_log(sprintf("LD矩阵对称性: %s (最大差异: %.2e)", 
                   ifelse(is_symmetric, "对称", "非对称"), max_diff))
  
  # 检查对角线元素（应该接近1，代表变异与自身的LD）
  diag_stats <- summary(diag(ld_mat))
  write_log(sprintf("对角线元素统计 - 最小值: %.3f, 平均值: %.3f, 最大值: %.3f", 
                   diag_stats["Min."], diag_stats["Mean"], diag_stats["Max."]))
  
  # =================================================================
  # 主导变异检测 (Lead Variant Detection)
  # =================================================================
  write_log("=== 主导变异位点检测 ===")
  
  # 初始化主导变异相关变量
  lead_variant_found <- FALSE    # 是否找到主导变异
  lead_r_values <- NULL          # 主导变异的LD相关系数向量
  lead_r2_values <- NULL         # 主导变异的LD决定系数向量
  
  # 尝试直接匹配主导变异ID
  if (lead_id %in% colnames(ld_r_with_ids)) {
    lead_variant_found <- TRUE
    lead_r_values <- ld_r_with_ids[, lead_id]
    lead_r2_values <- ld_r2_with_ids[, lead_id]
    write_log(sprintf("主导变异找到 (直接匹配): %s", lead_id))
  } else {
    # 尝试格式转换：冒号转点号（某些数据格式可能使用点号）
    lead_id_dots <- gsub(":", ".", lead_id, fixed=TRUE)
    if (lead_id_dots %in% colnames(ld_r_with_ids)) {
      lead_variant_found <- TRUE
      lead_r_values <- ld_r_with_ids[, lead_id_dots]
      lead_r2_values <- ld_r2_with_ids[, lead_id_dots]
      write_log(sprintf("主导变异找到 (格式转换): %s -> %s", lead_id, lead_id_dots))
    } else {
      write_log(sprintf("未找到主导变异: %s", lead_id), "警告")
      write_log(sprintf("格式转换后也未找到: %s", lead_id_dots), "警告")
      
      # 搜索相似的ID（可能存在命名不一致）
      similar_ids <- grep(gsub(":", ".", lead_id, fixed=TRUE), colnames(ld_r_with_ids), value=TRUE)
      if (length(similar_ids) > 0) {
        write_log(sprintf("发现相似的ID: %s", paste(head(similar_ids, 5), collapse=", ")), "警告")
      }
    }
  }
  
  # 如果找到主导变异，记录其LD值信息
  if (lead_variant_found) {
    write_log(sprintf("主导变异LD值数量: %d 个", length(lead_r_values)))
  }
  
  # =================================================================
  # SuSiE 算法执行 (SuSiE Algorithm Execution)
  # =================================================================
  write_log("=== 运行 SuSiE 分析 ===")
  
  # 设置随机种子确保结果可重现
  set.seed(seed)
  write_log(sprintf("设置随机种子: %d", seed))
  
  # 执行SuSiE RSS (Reference panel Summary Statistics)分析
  # bhat: 效应量估计值 (BETA)
  # shat: 效应量标准误 (SE) 
  # R: LD相关系数矩阵
  # n: 样本量
  # L: 最大独立信号数量
  susie_res <- susie_rss(
    bhat = sumstats$BETA,    # 变异效应量
    shat = sumstats$SE,      # 效应量标准误
    R = ld_mat,              # LD矩阵
    n = sample_size,         # 样本量
    L = max_components       # 最大组件数
  )
  
  write_log(sprintf("SuSiE分析完成 - 收敛状态: %s, 迭代次数: %d", 
                   ifelse(susie_res$converged, "已收敛", "未收敛"), susie_res$niter))
  
  # =================================================================
  # 结果提取 (Results Extraction)
  # =================================================================
  write_log("=== 提取分析结果 ===")
  
  # 提取后验包含概率 (Posterior Inclusion Probability, PIP)
  # PIP表示每个变异位点是因果变异的概率
  pip_values <- susie_get_pip(susie_res)
  
  # 提取可信集合 (Credible Sets)
  # 可信集合包含真实因果变异的概率达到指定覆盖度
  credible_sets <- susie_get_cs(susie_res, coverage = coverage, X = ld_mat)
  
  # 计算关键统计量
  n_credible_sets <- length(credible_sets$cs %||% list())  # 可信集合数量
  max_pip <- max(pip_values)                               # 最大PIP值
  n_variants_high_pip <- sum(pip_values > 0.1)           # 高PIP变异数量
  n_variants_very_high_pip <- sum(pip_values > 0.5)      # 极高PIP变异数量
  
  write_log(sprintf("识别的可信集合数量: %d", n_credible_sets))
  write_log(sprintf("最大PIP值: %.4f", max_pip))
  write_log(sprintf("PIP > 0.1 的变异数量: %d", n_variants_high_pip))
  write_log(sprintf("PIP > 0.5 的变异数量: %d", n_variants_very_high_pip))
  
  # 检查纯度矩阵是否可用（用于评估可信集合质量）
  if (!is.null(credible_sets$purity)) {
    write_log(sprintf("纯度矩阵可用: %d x %d", 
                     nrow(credible_sets$purity), ncol(credible_sets$purity)))
  }
  
  # =================================================================
  # 结果处理和整理 (Results Processing and Organization)
  # =================================================================
  write_log("=== 处理输出结果 ===")
  
  # 初始化结果列表结构
  results <- list(
    lead_id = lead_id,                    # 主导变异ID
    sample_size = sample_size,            # 样本量
    n_variants = nrow(sumstats),          # 变异位点总数
    parameters = list(                    # 分析参数
      pip_threshold = pip_threshold,      # PIP阈值
      max_components = max_components,    # 最大组件数
      coverage = coverage,                # 覆盖度
      seed = seed                         # 随机种子
    ),
    pip = list(),                         # PIP结果列表
    credible_sets = list()                # 可信集合列表
  )
  
  # =================================================================
  # PIP信息处理 (PIP Information Processing)
  # =================================================================
  # 添加超过阈值的变异位点PIP信息
  variants_above_threshold <- 0
  for (i in 1:length(pip_values)) {
    if (pip_values[i] >= pip_threshold) {  # 只处理超过PIP阈值的变异
      variants_above_threshold <- variants_above_threshold + 1
      # 获取变异ID（优先使用SNPID列，否则生成通用ID）
      variant_id <- if("SNPID" %in% colnames(sumstats)) sumstats$SNPID[i] else paste0("variant_", i)
      
      # 提取当前变异与主导变异的LD信息
      ld_r_with_lead <- NA     # LD相关系数
      ld_r2_with_lead <- NA    # LD决定系数
      if (lead_variant_found) {
        # 尝试直接从lead_r_values中获取（通过变异名称索引）
        if (variant_id %in% names(lead_r_values)) {
          ld_r_with_lead <- lead_r_values[variant_id]
          ld_r2_with_lead <- lead_r2_values[variant_id]
        } else {
          # 通过矩阵行索引方式获取LD值
          row_idx <- which(rownames(ld_r_with_ids) == variant_id)
          if (length(row_idx) > 0) {
            col_name <- if (lead_id %in% colnames(ld_r_with_ids)) lead_id else gsub(":", ".", lead_id, fixed=TRUE)
            ld_r_with_lead <- ld_r_with_ids[row_idx, col_name]
            ld_r2_with_lead <- ld_r2_with_ids[row_idx, col_name]
          }
        }
      }
      
      # 构建变异信息记录
      variant_info <- list(
        index = i,                           # 变异在数据中的索引
        variant_id = variant_id,             # 变异标识符
        pip = pip_values[i],                 # 后验包含概率
        beta = sumstats$BETA[i],             # 效应量
        se = sumstats$SE[i],                 # 标准误
        pvalue = if("P" %in% colnames(sumstats)) sumstats$P[i] else NA,  # P值
        ld_r_with_lead = ld_r_with_lead,     # 与主导变异的LD相关系数
        ld_r2_with_lead = ld_r2_with_lead    # 与主导变异的LD决定系数
      )
      results$pip[[length(results$pip) + 1]] <- variant_info
    }
  }
  
  write_log(sprintf("超过PIP阈值 (%.3f) 的变异数量: %d", pip_threshold, variants_above_threshold))
  
  # =================================================================
  # 可信集合处理 (Credible Sets Processing)  
  # =================================================================
  # 处理每个识别出的可信集合
  if (!is.null(credible_sets$cs)) {
    for (cs_idx in 1:length(credible_sets$cs)) {
      cs_variants <- credible_sets$cs[[cs_idx]]  # 当前可信集合中的变异索引
      
      # 提取纯度信息（衡量可信集合中变异间的相关性）
      purity_info <- NULL
      if (!is.null(credible_sets$purity) && nrow(credible_sets$purity) >= cs_idx) {
        purity_row <- credible_sets$purity[cs_idx, ]
        purity_info <- list(
          min_abs_corr = purity_row["min.abs.corr"],       # 最小绝对相关系数
          mean_abs_corr = purity_row["mean.abs.corr"],     # 平均绝对相关系数
          median_abs_corr = purity_row["median.abs.corr"]  # 中位绝对相关系数
        )
      }
      
      # 构建可信集合信息记录
      cs_info <- list(
        cs_index = cs_idx,                            # 可信集合索引
        coverage = credible_sets$coverage[cs_idx],    # 覆盖度
        purity = purity_info,                         # 纯度信息
        size = length(cs_variants),                   # 可信集合大小
        variants = list()                             # 变异列表
      )
      
      # 处理可信集合中的每个变异
      for (var_idx in cs_variants) {
        # 获取变异ID
        variant_id <- if("SNPID" %in% colnames(sumstats)) sumstats$SNPID[var_idx] else paste0("variant_", var_idx)
        
        # 提取与主导变异的LD信息
        ld_r_with_lead <- NA
        ld_r2_with_lead <- NA
        if (lead_variant_found) {
          if (variant_id %in% names(lead_r_values)) {
            ld_r_with_lead <- lead_r_values[variant_id]
            ld_r2_with_lead <- lead_r2_values[variant_id]
          } else {
            # 通过矩阵索引获取LD值
            row_idx <- which(rownames(ld_r_with_ids) == variant_id)
            if (length(row_idx) > 0) {
              col_name <- if (lead_id %in% colnames(ld_r_with_ids)) lead_id else gsub(":", ".", lead_id, fixed=TRUE)
              ld_r_with_lead <- ld_r_with_ids[row_idx, col_name]
              ld_r2_with_lead <- ld_r2_with_ids[row_idx, col_name]
            }
          }
        }
        
        # 构建可信集合中变异的详细信息
        variant_cs_info <- list(
          index = var_idx,                         # 变异索引
          variant_id = variant_id,                 # 变异标识符
          pip = pip_values[var_idx],               # 后验包含概率
          beta = sumstats$BETA[var_idx],           # 效应量
          se = sumstats$SE[var_idx],               # 标准误
          pvalue = if("P" %in% colnames(sumstats)) sumstats$P[var_idx] else NA,  # P值
          ld_r_with_lead = ld_r_with_lead,         # 与主导变异的LD相关系数
          ld_r2_with_lead = ld_r2_with_lead        # 与主导变异的LD决定系数
        )
        cs_info$variants[[length(cs_info$variants) + 1]] <- variant_cs_info
      }
      
      # 将可信集合信息添加到结果中
      results$credible_sets[[length(results$credible_sets) + 1]] <- cs_info
      
      write_log(sprintf("可信集合 %d: %d 个变异, 覆盖度=%.3f", 
                       cs_idx, length(cs_variants), credible_sets$coverage[cs_idx]))
    }
  }
  
  # =================================================================
  # 添加汇总信息 (Add Summary Information)
  # =================================================================
  results$summary <- list(
    converged = susie_res$converged,                # 是否收敛
    niter = susie_res$niter,                        # 迭代次数
    elbo = tail(susie_res$elbo, 1),                # 证据下界(ELBO)
    n_credible_sets = n_credible_sets,              # 可信集合数量
    max_pip = max_pip,                              # 最大PIP值
    n_variants_pip_gt_0.1 = n_variants_high_pip,   # PIP>0.1的变异数
    n_variants_pip_gt_0.5 = n_variants_very_high_pip # PIP>0.5的变异数
  )
  
  # =================================================================
  # 保存结果 (Save Results)
  # =================================================================
  write_log("=== 保存结果 ===")
  
  # 构建输出文件名（将冒号替换为下划线）
  output_filename <- paste0(gsub(":", "_", lead_id), ".susie_results.json")
  output_path <- file.path(out_dir, output_filename)
  
  # 将结果以JSON格式写入文件
  writeLines(toJSON(results, pretty = TRUE, auto_unbox = TRUE), output_path)
  write_log(sprintf("结果已保存至: %s", output_path))
  
  # =================================================================
  # 最终总结 (Final Summary)
  # =================================================================
  write_log("=== 分析总结 ===")
  write_log(sprintf("主导变异 %s 的分析已成功完成", lead_id))
  write_log(sprintf("收敛状态: %s", ifelse(susie_res$converged, "已收敛", "未收敛")))
  write_log(sprintf("迭代次数: %d", susie_res$niter))
  write_log(sprintf("证据下界(ELBO): %.4f", tail(susie_res$elbo, 1)))
  write_log(sprintf("识别的可信集合数量: %d", n_credible_sets))
  write_log(sprintf("最大后验包含概率(PIP): %.4f", max_pip))
  write_log(sprintf("高置信度变异 (PIP > 0.1): %d 个", n_variants_high_pip))
  write_log(sprintf("极高置信度变异 (PIP > 0.5): %d 个", n_variants_very_high_pip))
  write_log("=== 分析完成 ===")
  
}, error = function(e) {
  # 错误处理
  write_log(sprintf("分析过程中发生错误: %s", e$message), "错误")
  quit(status=1)
}, finally = {
  # 无论成功或失败都要关闭日志文件连接
  close(log_conn)
})

# =================================================================
# 程序成功退出 (Successful Exit)
# =================================================================
quit(status=0)
