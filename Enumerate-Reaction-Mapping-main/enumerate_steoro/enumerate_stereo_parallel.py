#!/usr/bin/env python3
"""
为反应SMILES枚举所有可能的立体异构体 - 多进程并行版本

对反应物和产物的立体异构体进行组合枚举
保留原子映射号信息
"""

import pandas as pd
import argparse
import os
import multiprocessing as mp
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers


def get_atom_map_number(atom):
    """获取原子的映射号"""
    if atom.HasProp('molAtomMapNumber'):
        return int(atom.GetProp('molAtomMapNumber'))
    return None


def set_atom_map_number(atom, map_num):
    """设置原子的映射号"""
    if map_num is not None:
        atom.SetProp('molAtomMapNumber', str(map_num))


def clear_atom_map_numbers(mol):
    """
    清除分子的原子映射号，返回 idx->map_num 的映射字典
    """
    idx_to_map = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        map_num = get_atom_map_number(atom)
        if map_num is not None:
            idx_to_map[idx] = map_num
        atom.ClearProp('molAtomMapNumber')
    return idx_to_map


def restore_atom_map_numbers(mol, idx_to_map):
    """根据 idx->map_num 映射字典恢复原子映射号"""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in idx_to_map:
            set_atom_map_number(atom, idx_to_map[idx])


def clear_non_tetrahedral_chiral_tags(mol):
    """
    清除非四配位原子上的手性标签。

    用于生成去重键，避免 RDKit 将非四配位原子（如部分 N）枚举成额外立体中心。
    """
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        if atom.GetDegree() != 4:
            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
            if atom.HasProp('_CIPCode'):
                atom.ClearProp('_CIPCode')


def enumerate_all_stereoisomers(smi, max_isomers=1000):
    """
    枚举SMILES的所有立体异构体

    Args:
        smi: SMILES字符串
        max_isomers: 最大立体异构体数量限制

    Returns:
        立体异构体SMILES列表
        - 如果枚举成功：只返回枚举出的立体异构体（不包含原始版本）
        - 如果枚举失败：返回原始SMILES作为保底
    """
    try:
        # 尝试解析SMILES，保留显式氢
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        mol = Chem.MolFromSmiles(smi, ps)

        if mol is None:
            return [smi]

        # 保存原始原子映射号
        idx_to_map = clear_atom_map_numbers(mol)

        # 枚举立体异构体
        try:
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(
                mol,
                options=EnumerateStereoisomers.StereoEnumerationOptions(
                    unique=True,
                    maxIsomers=max_isomers,
                    rand=False,
                    onlyUnassigned=False,
                    tryEmbedding=True,
                    onlyStereoGroups=False
                )
            ))
        except Exception as e:
            # 如果枚举失败，返回原始SMILES
            restore_atom_map_numbers(mol, idx_to_map)
            original_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            return [original_smi if original_smi else smi]

        # 如果没有额外的立体异构体，返回原始SMILES
        if not isomers:
            restore_atom_map_numbers(mol, idx_to_map)
            original_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            return [original_smi if original_smi else smi]

        # 枚举成功：清理非四配位原子的手性标签后，再去重并返回规范化结果
        seen = set()
        unique_results = []
        for isomer in isomers:
            restore_atom_map_numbers(isomer, idx_to_map)
            normalized_isomer = Chem.Mol(isomer)
            clear_non_tetrahedral_chiral_tags(normalized_isomer)
            normalized_smi = Chem.MolToSmiles(normalized_isomer, isomericSmiles=True)
            if not normalized_smi:
                normalized_smi = Chem.MolToSmiles(isomer, isomericSmiles=True)
            if not normalized_smi:
                continue

            if normalized_smi in seen:
                continue

            seen.add(normalized_smi)
            unique_results.append(normalized_smi)

        return unique_results

    except Exception as e:
        return [smi]


def enumerate_stereochemistry_for_reaction(reaction_smi, max_isomers=1000):
    """
    为反应SMILES枚举所有可能的立体异构体组合

    Args:
        reaction_smi: 反应SMILES字符串 (格式: reactants>>products)
        max_isomers: 每个分子的最大立体异构体数量限制

    Returns:
        所有立体异构体组合的反应SMILES列表
    """
    try:
        # 分离反应物和生成物
        parts = reaction_smi.split('>>')
        if len(parts) != 2:
            return [reaction_smi]

        reactants_smi, products_smi = parts

        # 处理反应物 - 枚举每个反应物分子的立体异构体
        reactants = reactants_smi.split('.')
        reactant_isomers_list = []
        for r_smi in reactants:
            isomers = enumerate_all_stereoisomers(r_smi, max_isomers)
            reactant_isomers_list.append(isomers)

        # 处理生成物 - 枚举每个生成物分子的立体异构体
        products = products_smi.split('.')
        product_isomers_list = []
        for p_smi in products:
            isomers = enumerate_all_stereoisomers(p_smi, max_isomers)
            product_isomers_list.append(isomers)

        # 生成所有可能的组合
        from itertools import product

        results = []

        # 反应物立体异构体组合
        for reactant_combo in product(*reactant_isomers_list):
            new_reactants = '.'.join(reactant_combo)

            # 生成物立体异构体组合
            for product_combo in product(*product_isomers_list):
                new_products = '.'.join(product_combo)
                results.append(new_reactants + '>>' + new_products)

        return results

    except Exception as e:
        return [reaction_smi]


def process_chunk(chunk_data, max_isomers, worker_id):
    """
    处理数据chunk的函数，用于多进程

    Args:
        chunk_data: 包含rows和start_idx的元组
        max_isomers: 每个分子的最大立体异构体数量限制
        worker_id: 工作进程ID

    Returns:
        处理结果列表
    """
    rows, start_idx = chunk_data
    results = []

    for row in rows:
        origin_idx = row['origin_idx']
        origin_rsmi = row['origin_rsmi']
        aug_idx = row['aug_idx']

        # 枚举立体异构体组合
        stereo_reactions = enumerate_stereochemistry_for_reaction(origin_rsmi, max_isomers)

        # 每个反应的立体异构体独立编号（从1开始）
        for stereo_idx, stereo_rsmi in enumerate(stereo_reactions, start=1):
            results.append({
                'origin_idx': origin_idx,
                'origin_rsmi': origin_rsmi,
                'aug_idx': aug_idx,
                'stereo_idx': stereo_idx,
                'stereo_rsmi': stereo_rsmi
            })

    return results, worker_id


def process_csv_parallel(input_csv, output_csv, max_isomers=1000, num_workers=None, debug=False):
    """
    多进程并行处理CSV文件

    Args:
        input_csv: 输入CSV文件路径
        output_csv: 输出CSV文件路径
        max_isomers: 每个分子的最大立体异构体数量限制
        num_workers: 工作进程数，None则使用所有CPU核心
        debug: 是否输出调试信息
    """
    if num_workers is None:
        num_workers = mp.cpu_count()

    df = pd.read_csv(input_csv)

    print(f"处理 {len(df)} 条反应...")
    print(f"每个分子最大立体异构体数量: {max_isomers}")
    print(f"使用 {num_workers} 个进程并行处理")

    # 将数据分成多个chunk
    chunk_size = (len(df) + num_workers - 1) // num_workers
    chunks = []

    for i in range(num_workers):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(df))
        chunk_rows = df.iloc[start_idx:end_idx].to_dict('records')
        chunks.append((chunk_rows, i * chunk_size + 1))

    # 使用进程池处理
    with mp.Pool(processes=num_workers) as pool:
        # 使用starmap传递多个参数
        args = [(chunk, max_isomers, i) for i, chunk in enumerate(chunks)]
        async_results = [pool.apply_async(process_chunk, arg) for arg in args]

        # 收集结果
        all_results = []
        completed = 0
        total = len(async_results)

        for async_result in async_results:
            chunk_results, worker_id = async_result.get()
            all_results.extend(chunk_results)
            completed += 1

            if debug:
                print(f"[Worker {worker_id}] 完成 ({completed}/{total})")
            elif completed % 5 == 0 or completed == total:
                print(f"进度: {completed}/{num_workers} 个chunk已完成")

    # 创建结果DataFrame并保存
    result_df = pd.DataFrame(all_results)
    result_df.to_csv(output_csv, index=False)
    print(f"\n已保存到: {output_csv}")
    print(f"总输出行数: {len(result_df)}")

    return result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='为反应SMILES枚举所有可能的立体异构体 (多进程版本)')
    parser.add_argument('input_csv', help='输入CSV文件路径')
    parser.add_argument('output_csv', help='输出CSV文件路径')
    parser.add_argument('--max-isomers', type=int, default=128,
                        help='每个分子的最大立体异构体数量限制 (默认: 128)')
    parser.add_argument('--num-workers', type=int, default=None,
                        help='工作进程数 (默认: 所有CPU核心)')
    parser.add_argument('--debug', action='store_true',
                        help='输出调试信息')

    args = parser.parse_args()

    process_csv_parallel(args.input_csv, args.output_csv, args.max_isomers, args.num_workers, args.debug)
