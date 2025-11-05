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

