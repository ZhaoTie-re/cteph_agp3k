import json
import pandas as pd
import numpy as np
import os
import subprocess
import shutil
import tempfile
from io import StringIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import logging


def setup_logging():
    """设置日志配置"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
        ]
    )
    return logging.getLogger(__name__)


def process_single_lead_variant(lead_variant, lead_data, output_dir):
    """
    处理单个lead variant的LD和统计数据
    
    Args:
        lead_variant (str): lead variant的ID
        lead_data (dict): 包含该lead variant相关文件路径的字典
        output_dir (str): 输出目录
    
    Returns:
        dict: 处理结果信息
    """
    logger = logging.getLogger(__name__)
    
    try:
        logger.info(f"开始处理lead variant: {lead_variant}")
        
        # 读取统计数据文件
        sum_stat_file = lead_data['sum_stat_tsv']
        if not os.path.exists(sum_stat_file):
            logger.error(f"统计数据文件不存在: {sum_stat_file}")
            return None
            
        logger.info(f"读取统计数据文件: {sum_stat_file}")
        sum_stat_df = pd.read_csv(sum_stat_file, sep='\t')
        
        # 检查必要的列
        if 'SNPID' not in sum_stat_df.columns:
            logger.error(f"统计数据文件缺少SNPID列: {sum_stat_file}")
            return None
        
        # 添加is_lead_variant列
        sum_stat_df['is_lead_variant'] = sum_stat_df['SNPID'] == lead_variant
        logger.info(f"添加is_lead_variant列，lead variant: {lead_variant}")
        
        # 读取LD r值数据
        ld_r_file = lead_data['ld_r_tsv']
        ld_r2_file = lead_data['ld_r2_tsv']
        
        # 初始化LD列
        sum_stat_df['ld_r_with_lead'] = np.nan
        sum_stat_df['ld_r2_with_lead'] = np.nan
        
        # 处理LD r值
        if os.path.exists(ld_r_file):
            logger.info(f"读取LD r值文件: {ld_r_file}")
            ld_r_df = pd.read_csv(ld_r_file, sep='\t', index_col=0)
            
            # 检查lead variant是否在LD矩阵中
            if lead_variant in ld_r_df.index and lead_variant in ld_r_df.columns:
                # 提取与lead variant的LD r值
                ld_r_values = ld_r_df.loc[:, lead_variant]
                
                # 将LD r值匹配到统计数据中
                for idx, row in sum_stat_df.iterrows():
                    snpid = row['SNPID']
                    if snpid in ld_r_values.index:
                        sum_stat_df.at[idx, 'ld_r_with_lead'] = ld_r_values[snpid]
                
                logger.info(f"成功添加LD r值，匹配到 {sum_stat_df['ld_r_with_lead'].notna().sum()} 个变体")
            else:
                logger.warning(f"Lead variant {lead_variant} 不在LD r矩阵中")
        else:
            logger.warning(f"LD r值文件不存在: {ld_r_file}")
        
        # 处理LD r2值
        if os.path.exists(ld_r2_file):
            logger.info(f"读取LD r2值文件: {ld_r2_file}")
            ld_r2_df = pd.read_csv(ld_r2_file, sep='\t', index_col=0)
            
            # 检查lead variant是否在LD矩阵中
            if lead_variant in ld_r2_df.index and lead_variant in ld_r2_df.columns:
                # 提取与lead variant的LD r2值
                ld_r2_values = ld_r2_df.loc[:, lead_variant]
                
                # 将LD r2值匹配到统计数据中
                for idx, row in sum_stat_df.iterrows():
                    snpid = row['SNPID']
                    if snpid in ld_r2_values.index:
                        sum_stat_df.at[idx, 'ld_r2_with_lead'] = ld_r2_values[snpid]
                
                logger.info(f"成功添加LD r2值，匹配到 {sum_stat_df['ld_r2_with_lead'].notna().sum()} 个变体")
            else:
                logger.warning(f"Lead variant {lead_variant} 不在LD r2矩阵中")
        else:
            logger.warning(f"LD r2值文件不存在: {ld_r2_file}")
        
        # 为lead variant本身设置LD值为1.0
        lead_mask = sum_stat_df['is_lead_variant']
        sum_stat_df.loc[lead_mask, 'ld_r_with_lead'] = 1.0
        sum_stat_df.loc[lead_mask, 'ld_r2_with_lead'] = 1.0
        
        # 生成输出文件名
        safe_variant_name = lead_variant.replace(':', '_').replace('>', '_').replace('<', '_')
        output_file = os.path.join(output_dir, f"{safe_variant_name}.ld_sum_stat.tsv")
        
        # 保存结果
        sum_stat_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"保存结果到: {output_file}")
        
        # 返回处理结果信息
        result = {
            'lead_variant': lead_variant,
            'output_file': output_file,
            'n_variants': int(len(sum_stat_df)),
            'n_variants_with_ld_r': int(sum_stat_df['ld_r_with_lead'].notna().sum()),
            'n_variants_with_ld_r2': int(sum_stat_df['ld_r2_with_lead'].notna().sum()),
            'processed_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        logger.info(f"完成处理lead variant: {lead_variant}")
        return result
        
    except Exception as e:
        logger.error(f"处理lead variant {lead_variant} 时发生错误: {str(e)}")
        return None


def process_ld_matrices_with_sum_stats(ld_json_file, output_dir=None, n_jobs=4):
    """
    处理LD矩阵和统计数据，为每个lead variant生成包含LD信息的统计文件
    
    Args:
        ld_json_file (str): LD矩阵摘要JSON文件路径
        output_dir (str, optional): 输出目录，默认为当前目录下的self_ld文件夹
        n_jobs (int): 并行处理的进程数，默认为4
    
    Returns:
        str: 输出摘要JSON文件路径
    """
    logger = setup_logging()
    
    # 设置输出目录
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'self_ld')
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"输出目录: {output_dir}")
    
    # 读取LD JSON文件
    logger.info(f"读取LD JSON文件: {ld_json_file}")
    try:
        with open(ld_json_file, 'r') as f:
            ld_data = json.load(f)
    except Exception as e:
        logger.error(f"读取LD JSON文件失败: {str(e)}")
        return None
    
    # 获取per_lead_outputs
    if 'per_lead_outputs' not in ld_data:
        logger.error("LD JSON文件中缺少per_lead_outputs字段")
        return None
    
    per_lead_outputs = ld_data['per_lead_outputs']
    lead_variants = list(per_lead_outputs.keys())
    logger.info(f"发现 {len(lead_variants)} 个lead variants: {lead_variants}")
    
    # 准备并行处理的参数
    process_args = [
        (lead_variant, lead_data, output_dir) 
        for lead_variant, lead_data in per_lead_outputs.items()
    ]
    
    # 并行处理
    results = []
    logger.info(f"开始并行处理，使用 {n_jobs} 个进程")
    
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # 提交任务
        future_to_variant = {
            executor.submit(process_single_lead_variant, *args): args[0] 
            for args in process_args
        }
        
        # 收集结果
        for future in as_completed(future_to_variant):
            lead_variant = future_to_variant[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    logger.info(f"✓ 成功处理: {lead_variant}")
                else:
                    logger.error(f"✗ 处理失败: {lead_variant}")
            except Exception as e:
                logger.error(f"✗ 处理 {lead_variant} 时发生异常: {str(e)}")
    
    # 生成输出摘要
    output_summary = {
        'created_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'source_ld_json': ld_json_file,
        'output_directory': output_dir,
        'n_lead_variants_total': int(len(lead_variants)),
        'n_lead_variants_processed': int(len(results)),
        'processing_parameters': {
            'n_jobs': int(n_jobs)
        },
        'per_lead_results': {result['lead_variant']: result for result in results}
    }
    
    # 保存输出摘要
    base_name = os.path.splitext(os.path.basename(ld_json_file))[0]
    output_summary_file = os.path.join(output_dir, f"{base_name}.ld_sum_stat_summary.json")
    
    with open(output_summary_file, 'w') as f:
        json.dump(output_summary, f, indent=2, ensure_ascii=False)
    
    logger.info(f"处理完成！")
    logger.info(f"- 总lead variants: {int(len(lead_variants))}")
    logger.info(f"- 成功处理: {int(len(results))}")
    logger.info(f"- 输出摘要文件: {output_summary_file}")
    
    return output_summary_file


def process_single_lead_variant_with_tommo(lead_variant, lead_data, output_dir, tommo_dir):
    """
    处理单个lead variant的统计数据并添加ToMMo LD信息
    
    Args:
        lead_variant (str): lead variant的ID (格式: chr:pos:ref:alt)
        lead_data (dict): 包含该lead variant相关文件路径的字典
        output_dir (str): 输出目录
        tommo_dir (str): ToMMo数据目录路径
    
    Returns:
        dict: 处理结果信息
    """
    logger = logging.getLogger(__name__)
    
    try:
        logger.info(f"开始处理lead variant (ToMMo): {lead_variant}")
        
        # 解析lead variant信息
        try:
            chr_name, pos, ref, alt = lead_variant.split(':')
            chr_num = chr_name.replace('chr', '')
        except ValueError:
            logger.error(f"Lead variant格式错误: {lead_variant}")
            return None
        
        # 读取统计数据文件
        sum_stat_file = lead_data['sum_stat_tsv']
        if not os.path.exists(sum_stat_file):
            logger.error(f"统计数据文件不存在: {sum_stat_file}")
            return None
            
        logger.info(f"读取统计数据文件: {sum_stat_file}")
        sum_stat_df = pd.read_csv(sum_stat_file, sep='\t')
        
        # 检查必要的列
        if 'SNPID' not in sum_stat_df.columns:
            logger.error(f"统计数据文件缺少SNPID列: {sum_stat_file}")
            return None
        
        # 添加is_lead_variant列
        sum_stat_df['is_lead_variant'] = sum_stat_df['SNPID'] == lead_variant
        logger.info(f"添加is_lead_variant列，lead variant: {lead_variant}")
        
        # 初始化ToMMo LD列
        sum_stat_df['tommo_r2_with_lead'] = np.nan
        
        # 初始化ToMMo记录计数
        tommo_lead_records = 0
        
        # 构建ToMMo文件路径
        tommo_file = os.path.join(tommo_dir, f"tommo-54kjpn-20230828-GRCh38-autosome-chr{chr_num}-plink-r2.tsv.gz")
        
        if not os.path.exists(tommo_file):
            logger.warning(f"ToMMo文件不存在: {tommo_file}")
        else:
            logger.info(f"查询ToMMo LD数据: {tommo_file}")
            
            # 使用tabix查询特定位点
            tabix_cmd = f"~/htslib-1.9/tabix {tommo_file} {chr_name}:{pos}-{pos}"
            
            try:
                result = subprocess.run(tabix_cmd, shell=True, capture_output=True, text=True)
                
                if result.returncode == 0 and result.stdout.strip():
                    lines = result.stdout.strip().split('\n')
                    logger.info(f"成功查询到ToMMo数据，记录数: {len(lines)}")
                    
                    # 解析tabix输出
                    tommo_records = []
                    
                    for line in lines:
                        if line.strip():
                            fields = line.split('\t')
                            if len(fields) >= 13:
                                # 构建variation1_id和variation2_id
                                var1_id = f"{fields[1]}:{fields[2]}:{fields[3]}:{fields[4]}"
                                var2_id = f"{fields[7]}:{fields[8]}:{fields[9]}:{fields[10]}"
                                r2_value = float(fields[12])
                                
                                tommo_records.append({
                                    'variation1_id': var1_id,
                                    'variation2_id': var2_id,
                                    'r2': r2_value
                                })
                    
                    logger.info(f"解析得到 {len(tommo_records)} 条ToMMo记录")
                    
                    # 查找与lead variant匹配的记录
                    matched_records = 0
                    
                    for record in tommo_records:
                        if record['variation1_id'] == lead_variant:
                            tommo_lead_records += 1  # 计数ToMMo中关于lead variant的记录
                            # 查找对应的SNP并填充r2值
                            mask = sum_stat_df['SNPID'] == record['variation2_id']
                            if mask.any():
                                sum_stat_df.loc[mask, 'tommo_r2_with_lead'] = record['r2']
                                matched_records += 1
                    
                    logger.info(f"ToMMo中关于lead variant的记录数: {tommo_lead_records}")
                    logger.info(f"匹配并填充了 {matched_records} 个变体的ToMMo r2值")
                    
                else:
                    logger.warning(f"tabix查询无结果或失败: {result.stderr}")
                    
            except Exception as e:
                logger.error(f"执行tabix查询时出错: {str(e)}")
        
        # 为lead variant本身设置r2值为1.0
        lead_mask = sum_stat_df['is_lead_variant']
        sum_stat_df.loc[lead_mask, 'tommo_r2_with_lead'] = 1.0
        
        # 生成输出文件名
        safe_variant_name = lead_variant.replace(':', '_').replace('>', '_').replace('<', '_')
        output_file = os.path.join(output_dir, f"{safe_variant_name}.tommo_ld_sum_stat.tsv")
        
        # 保存结果
        sum_stat_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"保存结果到: {output_file}")
        
        # 返回处理结果信息
        result = {
            'lead_variant': lead_variant,
            'output_file': output_file,
            'n_variants': int(len(sum_stat_df)),
            'n_variants_with_tommo_r2': int(sum_stat_df['tommo_r2_with_lead'].notna().sum()),
            'n_tommo_lead_records': int(tommo_lead_records),  # ToMMo中关于该lead variant的记录数
            'processed_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        logger.info(f"完成处理lead variant (ToMMo): {lead_variant}")
        return result
        
    except Exception as e:
        logger.error(f"处理lead variant {lead_variant} 时发生错误: {str(e)}")
        return None


def process_ld_matrices_with_tommo(ld_json_file, tommo_dir=None, output_dir=None, n_jobs=4):
    """
    处理LD矩阵和统计数据，为每个lead variant添加ToMMo LD信息
    
    Args:
        ld_json_file (str): LD矩阵摘要JSON文件路径
        tommo_dir (str, optional): ToMMo数据目录，默认为ToMMo_60KJPN/co-occurrence
        output_dir (str, optional): 输出目录，默认为当前目录下的tommo_ld文件夹
        n_jobs (int): 并行处理的进程数，默认为4
    
    Returns:
        str: 输出摘要JSON文件路径
    """
    logger = setup_logging()
    
    # 设置默认路径
    if tommo_dir is None:
        tommo_dir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/co-occurrence"
    
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'tommo_ld')
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"ToMMo数据目录: {tommo_dir}")
    
    # 检查ToMMo目录是否存在
    if not os.path.exists(tommo_dir):
        logger.error(f"ToMMo数据目录不存在: {tommo_dir}")
        return None
    
    # 读取LD JSON文件
    logger.info(f"读取LD JSON文件: {ld_json_file}")
    try:
        with open(ld_json_file, 'r') as f:
            ld_data = json.load(f)
    except Exception as e:
        logger.error(f"读取LD JSON文件失败: {str(e)}")
        return None
    
    # 获取per_lead_outputs
    if 'per_lead_outputs' not in ld_data:
        logger.error("LD JSON文件中缺少per_lead_outputs字段")
        return None
    
    per_lead_outputs = ld_data['per_lead_outputs']
    lead_variants = list(per_lead_outputs.keys())
    logger.info(f"发现 {len(lead_variants)} 个lead variants: {lead_variants}")
    
    # 准备并行处理的参数
    process_args = [
        (lead_variant, lead_data, output_dir, tommo_dir) 
        for lead_variant, lead_data in per_lead_outputs.items()
    ]
    
    # 并行处理
    results = []
    logger.info(f"开始并行处理，使用 {n_jobs} 个进程")
    
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # 提交任务
        future_to_variant = {
            executor.submit(process_single_lead_variant_with_tommo, *args): args[0] 
            for args in process_args
        }
        
        # 收集结果
        for future in as_completed(future_to_variant):
            lead_variant = future_to_variant[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    logger.info(f"✓ 成功处理(ToMMo): {lead_variant}")
                else:
                    logger.error(f"✗ 处理失败(ToMMo): {lead_variant}")
            except Exception as e:
                logger.error(f"✗ 处理 {lead_variant} 时发生异常: {str(e)}")
    
    # 生成输出摘要
    total_tommo_records = sum(result.get('n_tommo_lead_records', 0) for result in results)
    total_matched_variants = sum(result.get('n_variants_with_tommo_r2', 0) for result in results)
    
    output_summary = {
        'created_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'source_ld_json': ld_json_file,
        'tommo_data_directory': tommo_dir,
        'output_directory': output_dir,
        'n_lead_variants_total': int(len(lead_variants)),
        'n_lead_variants_processed': int(len(results)),
        'n_tommo_records_total': int(total_tommo_records),  # 所有lead variants在ToMMo中的总记录数
        'n_variants_with_tommo_r2_total': int(total_matched_variants),  # 总匹配的变体数
        'processing_parameters': {
            'n_jobs': int(n_jobs)
        },
        'per_lead_results': {result['lead_variant']: result for result in results}
    }
    
    # 保存输出摘要
    base_name = os.path.splitext(os.path.basename(ld_json_file))[0]
    output_summary_file = os.path.join(output_dir, f"{base_name}.tommo_ld_sum_stat_summary.json")
    
    with open(output_summary_file, 'w') as f:
        json.dump(output_summary, f, indent=2, ensure_ascii=False)
    
    logger.info(f"ToMMo LD处理完成！")
    logger.info(f"- 总lead variants: {int(len(lead_variants))}")
    logger.info(f"- 成功处理: {int(len(results))}")
    
    # 计算ToMMo记录统计
    total_tommo_records = sum(result.get('n_tommo_lead_records', 0) for result in results)
    logger.info(f"- ToMMo中总记录数: {total_tommo_records}")
    
    logger.info(f"- 输出摘要文件: {output_summary_file}")
    
    return output_summary_file

# new function achivment to add EAS LD info

def process_single_lead_variant_with_eas(lead_variant, lead_data, output_dir, eas_bed_prefix, intersection_summary):
    """
    处理单个lead variant的统计数据并添加EAS 1000G LD信息
    
    Args:
        lead_variant (str): lead variant的ID (格式: chr:pos:ref:alt)
        lead_data (dict): 包含该lead variant相关文件路径的字典
        output_dir (str): 输出目录
        eas_bed_prefix (str): EAS plink bed文件的前缀路径
        intersection_summary (dict): 变体交集摘要信息
    
    Returns:
        dict: 处理结果信息
    """
    logger = logging.getLogger(__name__)
    
    try:
        logger.info(f"开始处理lead variant (EAS): {lead_variant}")
        
        # 读取统计数据文件
        sum_stat_file = lead_data['sum_stat_tsv']
        if not os.path.exists(sum_stat_file):
            logger.error(f"统计数据文件不存在: {sum_stat_file}")
            return None
            
        logger.info(f"读取统计数据文件: {sum_stat_file}")
        sum_stat_df = pd.read_csv(sum_stat_file, sep='\t')
        
        # 检查必要的列
        if 'SNPID' not in sum_stat_df.columns:
            logger.error(f"统计数据文件缺少SNPID列: {sum_stat_file}")
            return None
        
        # 添加is_lead_variant列
        sum_stat_df['is_lead_variant'] = sum_stat_df['SNPID'] == lead_variant
        logger.info(f"添加is_lead_variant列，lead variant: {lead_variant}")
        
        # 初始化EAS LD列
        sum_stat_df['eas_r2_with_lead'] = np.nan
        
        # 检查lead variant是否在交集中
        lead_intersect_info = intersection_summary['per_lead_intersections'].get(lead_variant)
        
        if lead_intersect_info is None:
            logger.warning(f"未找到lead variant {lead_variant} 的交集信息，保持eas_r2_with_lead为NaN")
        elif not lead_intersect_info.get('lead_variant_in_intersection', False):
            logger.warning(f"Lead variant {lead_variant} 不在EAS交集中，保持eas_r2_with_lead为NaN")
        else:
            # Lead variant在交集中，为其设置r2值为1.0
            lead_mask = sum_stat_df['is_lead_variant']
            sum_stat_df.loc[lead_mask, 'eas_r2_with_lead'] = 1.0
            
            # 读取预计算的交集文件
            intersect_file = lead_intersect_info['intersect_file']
            if not os.path.exists(intersect_file):
                logger.error(f"交集文件不存在: {intersect_file}")
            else:
                logger.info(f"读取交集文件: {intersect_file}")
                with open(intersect_file, 'r') as f:
                    common_variants = set(line.strip() for line in f if line.strip())
                
                logger.info(f"交集文件中有 {len(common_variants)} 个共同变体")
                logger.info(f"交集文件中有 {len(common_variants)} 个共同变体")
                
                if len(common_variants) == 0:
                    logger.warning("交集为空，跳过LD计算")
                else:
                    # 创建临时目录用于此lead variant的计算
                    safe_variant_name = lead_variant.replace(':', '_').replace('>', '_').replace('<', '_')
                    variant_tmp_dir = os.path.join(os.getcwd(), 'eas_ld_tmp', safe_variant_name)
                    os.makedirs(variant_tmp_dir, exist_ok=True)
                    
                    logger.info(f"临时目录: {variant_tmp_dir}")
                    
                    try:
                        # 使用plink2提取子集，plink1.9计算LD
                        subset_prefix = os.path.join(variant_tmp_dir, "eas_subset")
                        plink2_path = "/home/b/b37974/plink2"
                        plink19_path = "/home/b/b37974/plink"
                        
                        # 步骤1: 使用plink2提取变体子集
                        extract_cmd = [
                            plink2_path,
                            "--bfile", eas_bed_prefix,
                            "--extract", intersect_file,
                            "--make-bed",
                            "--out", subset_prefix
                        ]
                        
                        logger.info("使用plink2提取变体子集...")
                        result = subprocess.run(extract_cmd, capture_output=True, text=True)
                        
                        if result.returncode != 0:
                            logger.error(f"plink2提取变体失败: {result.stderr}")
                        elif not os.path.exists(f"{subset_prefix}.bed"):
                            logger.error(f"子集bed文件未生成: {subset_prefix}.bed")
                        else:
                            # 步骤2: 使用plink1.9计算LD矩阵
                            ld_prefix = os.path.join(variant_tmp_dir, "ld_matrix")
                            ld_cmd = [
                                plink19_path,
                                "--bfile", subset_prefix,
                                "--r2", "square",
                                "--keep-allele-order",
                                "--out", ld_prefix
                            ]
                            
                            logger.info("使用plink1.9计算LD矩阵...")
                            result = subprocess.run(ld_cmd, capture_output=True, text=True)
                            
                            if result.returncode != 0:
                                logger.error(f"plink1.9计算LD失败: {result.stderr}")
                            else:
                                # 读取LD矩阵结果
                                ld_file = f"{ld_prefix}.ld"
                                if not os.path.exists(ld_file):
                                    logger.error(f"LD矩阵文件未生成: {ld_file}")
                                else:
                                    logger.info(f"LD矩阵文件已保存: {ld_file}")
                                    
                                    # 读取子集的bim文件以获取变体顺序
                                    subset_bim = f"{subset_prefix}.bim"
                                    if not os.path.exists(subset_bim):
                                        logger.error(f"子集bim文件不存在: {subset_bim}")
                                    else:
                                        # 读取bim文件获取变体ID列表
                                        bim_df = pd.read_csv(subset_bim, sep='\t', header=None,
                                                            names=['chr', 'snp_id', 'cm', 'pos', 'a1', 'a2'])
                                        variant_ids = bim_df['snp_id'].tolist()
                                        logger.info(f"从bim文件读取到 {len(variant_ids)} 个变体ID")
                                        
                                        # 读取LD矩阵（无表头的方阵）
                                        ld_matrix = pd.read_csv(ld_file, sep=r'\s+', header=None)
                                        logger.info(f"LD矩阵维度: {ld_matrix.shape}")
                                        
                                        # 验证矩阵维度与变体数量匹配
                                        if ld_matrix.shape[0] != len(variant_ids) or ld_matrix.shape[1] != len(variant_ids):
                                            logger.error(f"LD矩阵维度 {ld_matrix.shape} 与变体数量 {len(variant_ids)} 不匹配")
                                        else:
                                            # 设置行名和列名
                                            ld_matrix = ld_matrix.set_axis(variant_ids, axis=0)
                                            ld_matrix = ld_matrix.set_axis(variant_ids, axis=1)
                                            logger.info(f"为LD矩阵添加了行名和列名")
                                            
                                            # 保存带标签的LD矩阵用于调试
                                            labeled_ld_file = os.path.join(variant_tmp_dir, "ld_matrix_labeled.tsv")
                                            ld_matrix.to_csv(labeled_ld_file, sep='\t')
                                            logger.info(f"带标签的LD矩阵已保存: {labeled_ld_file}")
                                        
                                            # 检查lead variant是否在LD矩阵中
                                            if lead_variant not in ld_matrix.index:
                                                logger.warning(f"Lead variant {lead_variant} 不在LD矩阵中")
                                            else:
                                                # 提取lead variant的LD值
                                                lead_ld_values = ld_matrix.loc[lead_variant, :]
                                                
                                                # 匹配LD值到统计数据
                                                matched_count = 0
                                                for snpid, r2_value in lead_ld_values.items():
                                                    mask = sum_stat_df['SNPID'] == snpid
                                                    if mask.any():
                                                        sum_stat_df.loc[mask, 'eas_r2_with_lead'] = float(r2_value)
                                                        matched_count += 1
                                                
                                                logger.info(f"成功匹配 {matched_count} 个变体的EAS LD值")
                        
                    except Exception as e:
                        logger.error(f"处理EAS LD时发生错误: {str(e)}")
                        import traceback
                        logger.error(traceback.format_exc())
                    
                    # 注意：不删除临时目录，保留用于调试
                    logger.info(f"临时文件保留在: {variant_tmp_dir}")
        
        # 保存结果
        safe_variant_name = lead_variant.replace(':', '_').replace('>', '_').replace('<', '_')
        output_file = os.path.join(output_dir, f"{safe_variant_name}.eas_ld_sum_stat.tsv")
        sum_stat_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"保存结果到: {output_file}")
        
        # 返回处理结果信息
        result = {
            'lead_variant': lead_variant,
            'output_file': output_file,
            'n_variants': int(len(sum_stat_df)),
            'n_variants_with_eas_r2': int(sum_stat_df['eas_r2_with_lead'].notna().sum()),
            'processed_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        logger.info(f"完成处理lead variant (EAS): {lead_variant}")
        return result
        
    except Exception as e:
        logger.error(f"处理lead variant {lead_variant} 时发生错误: {str(e)}")
        return None


def process_ld_matrices_with_eas(ld_json_file, eas_bed_prefix=None, intersection_json=None, output_dir=None, n_jobs=4):
    """
    处理LD矩阵和统计数据，为每个lead variant添加EAS 1000G LD信息
    
    Args:
        ld_json_file (str): LD矩阵摘要JSON文件路径
        eas_bed_prefix (str, optional): EAS plink bed文件的前缀路径（不含.bed/.bim/.fam后缀）
        intersection_json (str, optional): 变体交集摘要JSON文件路径
        output_dir (str, optional): 输出目录，默认为当前目录下的eas_ld文件夹
        n_jobs (int): 并行处理的进程数，默认为4
    
    Returns:
        str: 输出摘要JSON文件路径
    """
    logger = setup_logging()
    
    # 设置默认路径
    if eas_bed_prefix is None:
        eas_bed_prefix = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all"
    
    if intersection_json is None:
        intersection_json = os.path.join(os.getcwd(), 'eas_ld_tmp', 'variant_intersection_summary.json')
    
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'eas_ld')
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"EAS bed前缀: {eas_bed_prefix}")
    logger.info(f"交集摘要JSON: {intersection_json}")
    
    # 检查必要文件和工具是否存在
    if not os.path.exists(f"{eas_bed_prefix}.bed"):
        logger.error(f"EAS bed文件不存在: {eas_bed_prefix}.bed")
        return None
    
    if not os.path.exists(f"{eas_bed_prefix}.bim"):
        logger.error(f"EAS bim文件不存在: {eas_bed_prefix}.bim")
        return None
    
    if not os.path.exists(f"{eas_bed_prefix}.fam"):
        logger.error(f"EAS fam文件不存在: {eas_bed_prefix}.fam")
        return None
    
    plink2_path = "/home/b/b37974/plink2"
    if not os.path.exists(plink2_path):
        logger.error(f"plink2不存在: {plink2_path}")
        return None
    
    plink19_path = "/home/b/b37974/plink"
    if not os.path.exists(plink19_path):
        logger.error(f"plink1.9不存在: {plink19_path}")
        return None
    
    # 读取交集摘要JSON
    if not os.path.exists(intersection_json):
        logger.error(f"交集摘要JSON文件不存在: {intersection_json}")
        return None
    
    logger.info(f"读取交集摘要JSON文件: {intersection_json}")
    try:
        with open(intersection_json, 'r') as f:
            intersection_summary = json.load(f)
    except Exception as e:
        logger.error(f"读取交集摘要JSON失败: {str(e)}")
        return None
    
    # 读取LD JSON文件
    logger.info(f"读取LD JSON文件: {ld_json_file}")
    try:
        with open(ld_json_file, 'r') as f:
            ld_data = json.load(f)
    except Exception as e:
        logger.error(f"读取LD JSON文件失败: {str(e)}")
        return None
    
    # 获取per_lead_outputs
    if 'per_lead_outputs' not in ld_data:
        logger.error("LD JSON文件中缺少per_lead_outputs字段")
        return None
    
    per_lead_outputs = ld_data['per_lead_outputs']
    lead_variants = list(per_lead_outputs.keys())
    logger.info(f"发现 {len(lead_variants)} 个lead variants: {lead_variants}")
    
    # 逐个处理lead variants（避免内存问题）
    results = []
    logger.info(f"开始顺序处理lead variants")
    
    for lead_variant, lead_data in per_lead_outputs.items():
        try:
            logger.info(f"处理 {lead_variant}...")
            result = process_single_lead_variant_with_eas(
                lead_variant, lead_data, output_dir, eas_bed_prefix, intersection_summary
            )
            if result is not None:
                results.append(result)
                logger.info(f"✓ 成功处理(EAS): {lead_variant}")
            else:
                logger.error(f"✗ 处理失败(EAS): {lead_variant}")
        except Exception as e:
            logger.error(f"✗ 处理 {lead_variant} 时发生异常: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
    
    # 生成输出摘要
    total_matched_variants = sum(result.get('n_variants_with_eas_r2', 0) for result in results)
    
    output_summary = {
        'created_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'source_ld_json': ld_json_file,
        'intersection_json': intersection_json,
        'eas_bed_prefix': eas_bed_prefix,
        'output_directory': output_dir,
        'n_lead_variants_total': int(len(lead_variants)),
        'n_lead_variants_processed': int(len(results)),
        'n_variants_with_eas_r2_total': int(total_matched_variants),
        'processing_parameters': {
            'n_jobs': int(n_jobs),
            'plink2_path': plink2_path,
            'plink19_path': plink19_path
        },
        'per_lead_results': {result['lead_variant']: result for result in results}
    }
    
    # 保存输出摘要
    base_name = os.path.splitext(os.path.basename(ld_json_file))[0]
    output_summary_file = os.path.join(output_dir, f"{base_name}.eas_ld_sum_stat_summary.json")
    
    with open(output_summary_file, 'w') as f:
        json.dump(output_summary, f, indent=2, ensure_ascii=False)
    
    logger.info(f"EAS LD处理完成！")
    logger.info(f"- 总lead variants: {int(len(lead_variants))}")
    logger.info(f"- 成功处理: {int(len(results))}")
    logger.info(f"- 总匹配变体数: {total_matched_variants}")
    logger.info(f"- 输出摘要文件: {output_summary_file}")
    
    return output_summary_file


def merge_ld_data_sources(self_ld_summary, tommo_ld_summary, eas_ld_summary, output_dir=None):
    """
    合并自己的LD数据、ToMMo LD数据和EAS 1000G LD数据
    
    Args:
        self_ld_summary (str): 自己LD数据的摘要JSON文件路径
        tommo_ld_summary (str): ToMMo LD数据的摘要JSON文件路径
        eas_ld_summary (str): EAS LD数据的摘要JSON文件路径
        output_dir (str, optional): 输出目录，默认为当前目录下的ld_summary文件夹
    
    Returns:
        str: 输出摘要JSON文件路径
    """
    logger = setup_logging()
    
    # 设置默认输出目录
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'ld_summary')
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"输出目录: {output_dir}")
    
    # 读取三个摘要JSON文件
    logger.info("读取LD摘要JSON文件...")
    
    try:
        with open(self_ld_summary, 'r') as f:
            self_data = json.load(f)
        logger.info(f"✓ 读取自己的LD数据摘要: {self_ld_summary}")
    except Exception as e:
        logger.error(f"读取自己的LD数据摘要失败: {str(e)}")
        return None
    
    try:
        with open(tommo_ld_summary, 'r') as f:
            tommo_data = json.load(f)
        logger.info(f"✓ 读取ToMMo LD数据摘要: {tommo_ld_summary}")
    except Exception as e:
        logger.error(f"读取ToMMo LD数据摘要失败: {str(e)}")
        return None
    
    try:
        with open(eas_ld_summary, 'r') as f:
            eas_data = json.load(f)
        logger.info(f"✓ 读取EAS LD数据摘要: {eas_ld_summary}")
    except Exception as e:
        logger.error(f"读取EAS LD数据摘要失败: {str(e)}")
        return None
    
    # 获取所有lead variants
    self_results = self_data.get('per_lead_results', {})
    tommo_results = tommo_data.get('per_lead_results', {})
    eas_results = eas_data.get('per_lead_results', {})
    
    lead_variants = list(self_results.keys())
    logger.info(f"发现 {len(lead_variants)} 个lead variants")
    
    # 处理每个lead variant
    merged_results = []
    
    for lead_variant in lead_variants:
        try:
            logger.info(f"处理lead variant: {lead_variant}")
            
            # 获取自己的LD数据文件（作为基础）
            self_result = self_results.get(lead_variant)
            if self_result is None:
                logger.warning(f"未找到lead variant {lead_variant} 的自己LD数据")
                continue
            
            self_file = self_result.get('output_file')
            if not os.path.exists(self_file):
                logger.error(f"自己的LD数据文件不存在: {self_file}")
                continue
            
            # 读取自己的LD数据（作为基础表）
            logger.info(f"  读取自己的LD数据: {self_file}")
            base_df = pd.read_csv(self_file, sep='\t')
            
            # 检查SNPID列
            if 'SNPID' not in base_df.columns:
                logger.error(f"自己的LD数据缺少SNPID列: {self_file}")
                continue
            
            logger.info(f"  基础数据: {len(base_df)} 个变体")
            
            # 添加ToMMo LD数据
            tommo_result = tommo_results.get(lead_variant)
            if tommo_result is not None:
                tommo_file = tommo_result.get('output_file')
                if os.path.exists(tommo_file):
                    logger.info(f"  读取ToMMo LD数据: {tommo_file}")
                    tommo_df = pd.read_csv(tommo_file, sep='\t')
                    
                    # 只保留SNPID和tommo_r2_with_lead列
                    if 'SNPID' in tommo_df.columns and 'tommo_r2_with_lead' in tommo_df.columns:
                        tommo_subset = tommo_df[['SNPID', 'tommo_r2_with_lead']].copy()
                        
                        # 合并到基础表
                        base_df = base_df.merge(tommo_subset, on='SNPID', how='left')
                        logger.info(f"  ✓ 添加ToMMo LD数据")
                    else:
                        logger.warning(f"  ToMMo数据缺少必要的列")
                        base_df['tommo_r2_with_lead'] = np.nan
                else:
                    logger.warning(f"  ToMMo LD数据文件不存在: {tommo_file}")
                    base_df['tommo_r2_with_lead'] = np.nan
            else:
                logger.warning(f"  未找到lead variant {lead_variant} 的ToMMo LD结果")
                base_df['tommo_r2_with_lead'] = np.nan
            
            # 添加EAS LD数据
            eas_result = eas_results.get(lead_variant)
            if eas_result is not None:
                eas_file = eas_result.get('output_file')
                if os.path.exists(eas_file):
                    logger.info(f"  读取EAS LD数据: {eas_file}")
                    eas_df = pd.read_csv(eas_file, sep='\t')
                    
                    # 只保留SNPID和eas_r2_with_lead列
                    if 'SNPID' in eas_df.columns and 'eas_r2_with_lead' in eas_df.columns:
                        eas_subset = eas_df[['SNPID', 'eas_r2_with_lead']].copy()
                        
                        # 合并到基础表
                        base_df = base_df.merge(eas_subset, on='SNPID', how='left')
                        logger.info(f"  ✓ 添加EAS LD数据")
                    else:
                        logger.warning(f"  EAS数据缺少必要的列")
                        base_df['eas_r2_with_lead'] = np.nan
                else:
                    logger.warning(f"  EAS LD数据文件不存在: {eas_file}")
                    base_df['eas_r2_with_lead'] = np.nan
            else:
                logger.warning(f"  未找到lead variant {lead_variant} 的EAS LD结果")
                base_df['eas_r2_with_lead'] = np.nan
            
            # 保存合并后的数据
            safe_variant_name = lead_variant.replace(':', '_').replace('>', '_').replace('<', '_')
            output_file = os.path.join(output_dir, f"{safe_variant_name}.merged_ld_sum_stat.tsv")
            base_df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"  保存合并数据到: {output_file}")
            
            # 统计信息
            result = {
                'lead_variant': lead_variant,
                'output_file': output_file,
                'n_variants': int(len(base_df)),
                'n_variants_with_self_ld_r2': int(base_df['ld_r2_with_lead'].notna().sum()) if 'ld_r2_with_lead' in base_df.columns else 0,
                'n_variants_with_tommo_r2': int(base_df['tommo_r2_with_lead'].notna().sum()) if 'tommo_r2_with_lead' in base_df.columns else 0,
                'n_variants_with_eas_r2': int(base_df['eas_r2_with_lead'].notna().sum()) if 'eas_r2_with_lead' in base_df.columns else 0,
                'processed_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            }
            
            merged_results.append(result)
            logger.info(f"✓ 成功合并lead variant {lead_variant} 的LD数据")
            
        except Exception as e:
            logger.error(f"处理lead variant {lead_variant} 时发生错误: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            continue
    
    # 生成输出摘要
    output_summary = {
        'created_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'source_summaries': {
            'self_ld_summary': self_ld_summary,
            'tommo_ld_summary': tommo_ld_summary,
            'eas_ld_summary': eas_ld_summary
        },
        'output_directory': output_dir,
        'n_lead_variants_total': int(len(lead_variants)),
        'n_lead_variants_merged': int(len(merged_results)),
        'per_lead_results': {result['lead_variant']: result for result in merged_results}
    }
    
    # 保存输出摘要
    output_summary_file = os.path.join(output_dir, 'merged_ld_sum_stat_summary.json')
    
    with open(output_summary_file, 'w') as f:
        json.dump(output_summary, f, indent=2, ensure_ascii=False)
    
    logger.info(f"\n合并LD数据完成！")
    logger.info(f"- 总lead variants: {int(len(lead_variants))}")
    logger.info(f"- 成功合并: {int(len(merged_results))}")
    logger.info(f"- 输出摘要文件: {output_summary_file}")
    
    # 打印详细统计
    if len(merged_results) > 0:
        total_self_ld = sum(r.get('n_variants_with_self_ld_r2', 0) for r in merged_results)
        total_tommo = sum(r.get('n_variants_with_tommo_r2', 0) for r in merged_results)
        total_eas = sum(r.get('n_variants_with_eas_r2', 0) for r in merged_results)
        
        logger.info(f"\n总体统计:")
        logger.info(f"- 包含自己LD r²数据的变体总数: {total_self_ld}")
        logger.info(f"- 包含ToMMo r²数据的变体总数: {total_tommo}")
        logger.info(f"- 包含EAS r²数据的变体总数: {total_eas}")
    
    return output_summary_file


# new funtion

