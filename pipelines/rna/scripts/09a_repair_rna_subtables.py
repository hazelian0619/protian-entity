#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from pathlib import Path

EXPECTED = [
    'rna_id','rna_type','rna_name','sequence','sequence_len','taxon_id',
    'symbol','source','fetch_date','source_version','hgnc_id',
    'ensembl_gene_id','ncbi_gene_id','mirbase_id','ensembl_transcript_id',
    'rfam_id','secondary_structure','pdb_ids'
]

def split_line(line: str):
    # 优先按 tab；如果没有 tab，就按任意空白
    if '\t' in line:
        return line.rstrip('\n').split('\t')
    return re.split(r'\s+', line.strip())

def repair_symbol_source(token: str):
    # 形如 MIR660RNAcentral;miRBase -> (MIR660, RNAcentral;miRBase)
    m = re.match(r'^(.*?)(RNAcentral;.*)$', token)
    if m:
        return m.group(1), m.group(2)
    return None

def repair_ncbi_sequence(token: str):
    # 形如 381.0CUGC... -> (381.0, CUGC...)
    m = re.match(r'^(\d+(?:\.\d+)?)([ACGUN]+)$', token)
    if m:
        return m.group(1), m.group(2)
    return None

def repair_sequence_len(token: str):
    # 形如 ACGU...1032 -> (ACGU..., 1032)
    m = re.match(r'^([ACGUN]+)(\d+)$', token)
    if m:
        return m.group(1), m.group(2)
    return None

def repair_file(inp: Path, outp: Path):
    with open(inp, 'r') as f:
        header_raw = f.readline().rstrip('\n')
        hdr = split_line(header_raw)

        # 如果 header 不是 18 列（你现在就是这种情况），直接用标准 header
        if len(hdr) != 18 or ('symbolsource' in hdr):
            header = EXPECTED
        else:
            header = hdr

        fixed = []
        bad = 0

        for line in f:
            if not line.strip():
                continue
            parts = split_line(line)

            # 常见问题1：header里/数据里把 symbol 和 source 粘成一个字段（导致少1列）
            if len(parts) == 17:
                # 尝试在第7列附近找包含 RNAcentral; 的粘连字段
                for i in range(len(parts)):
                    rep = repair_symbol_source(parts[i])
                    if rep:
                        symbol, source = rep
                        parts = parts[:i] + [symbol, source] + parts[i+1:]
                        break

            # 常见问题2：ncbi_gene_id + sequence 粘连（多发生在中间文件/某些输出）
            if len(parts) == 17:
                for i in range(len(parts)):
                    rep = repair_ncbi_sequence(parts[i])
                    if rep:
                        ncbi, seq = rep
                        parts = parts[:i] + [ncbi, seq] + parts[i+1:]
                        break

            # 常见问题3：sequence + sequence_len 粘连
            if len(parts) == 17:
                for i in range(len(parts)):
                    rep = repair_sequence_len(parts[i])
                    if rep:
                        seq, seqlen = rep
                        parts = parts[:i] + [seq, seqlen] + parts[i+1:]
                        break

            # 最终：必须 18 列，否则丢到 bad
            if len(parts) != 18:
                bad += 1
                continue

            fixed.append(parts)

    with open(outp, 'w') as g:
        g.write('\t'.join(EXPECTED) + '\n')
        for row in fixed:
            # 统一补齐空值
            row = [x if x is not None else '' for x in row]
            g.write('\t'.join(row) + '\n')

    print(f"[OK] {inp} -> {outp}")
    print(f"[OK] 修复后记录: {len(fixed)}，丢弃异常行: {bad}")

def main():
    mirna_in = Path("data/output/rna_master_mirna_v1.tsv")
    mrna_in  = Path("data/output/rna_master_mrna_v1.tsv")
    mirna_out = Path("data/output/rna_master_mirna_v1.fixed.tsv")
    mrna_out  = Path("data/output/rna_master_mrna_v1.fixed.tsv")

    repair_file(mirna_in, mirna_out)
    repair_file(mrna_in, mrna_out)

if __name__ == "__main__":
    main()