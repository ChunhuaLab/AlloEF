import pandas as pd

# 定义氨基酸三字母代码到单字母代码的映射
three_to_one = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
}

# 将三字母代码转换为单字母代码
def convert_three_to_one(residues):
    return ''.join(three_to_one[residue] for residue in residues if residue in three_to_one)

# 描述性标题映射
index_to_description = {
    'LEVM780105': 'Normalized Accessibility (Na)',
    'CHOP780215': 'Normalized Average Non-bonded Energy per Atom (Nec)',
    'HOPT810101': 'Normalized Hydrophobicity (Nphb)',
    'KYTJ820101': 'Hydrophilicity (Hdrpo)',
    'CHOP780216': 'Normalized Average Polarizability (Hdrpi)',
    'GRAR740102': 'Propensity (Prop)',
    'JANJ790102': 'Isoelectric Point (Isoep)',
    'FUKS010101': 'Relative Molecular Mass (Mass)',
    'ZHOH040101': 'Entropy Encoding (Enc)',
    'ZHOH040102': 'Electron-ion Interaction Potential (Eiip)'
}

# 排除的索引
excluded_indexes = {'AVBF000109', 'YANJ020101', 'GUYH850103', 'ROSM880105'}

# 解析Aaindex文件
def parse_aaindex1(aaindex1_file):
    aaindex_dict = {}
    with open(aaindex1_file, 'r') as file:
        lines = file.readlines()

    current_index = None
    amino_acids = "ARNDCEQGHILKMFPSTWYV"
    current_values = []
    within_index_data = False

    for line in lines:
        if line.startswith("H "):
            if current_index and current_values:
                if len(current_values) == len(amino_acids):
                    aaindex_dict[current_index]['index'] = dict(zip(amino_acids, current_values))
                else:
                    print(f"注意：索引 {current_index} 的值数目与氨基酸数目不匹配，已跳过该索引。")
            current_index = line[2:].strip()
            aaindex_dict[current_index] = {}
            current_values = []
            within_index_data = False
        elif line.startswith("D "):
            aaindex_dict[current_index]['description'] = line[2:].strip()
        elif line.startswith("R "):
            aaindex_dict[current_index]['authors'] = line[2:].strip()
        elif line.startswith("A "):
            aaindex_dict[current_index]['pubmed_id'] = line[2:].strip()
        elif line.startswith("I "):
            within_index_data = True
            continue  # 跳过 "I " 行而直接处理下一行开始的数据
        elif within_index_data:
            # 尝试将当前行转换为数值列表，忽略不可转换的内容
            try:
                current_values.extend(list(map(float, line.strip().split())))
            except ValueError:
                continue

    # 确保最后一个索引的数据也被解析
    if current_index and current_values:
        if len(current_values) == len(amino_acids):
            aaindex_dict[current_index]['index'] = dict(zip(amino_acids, current_values))
        else:
            print(f"注意：索引 {current_index} 的值数目与氨基酸数目不匹配，已跳过该索引。")

    return aaindex_dict

# 读取Excel文件
def read_excel(file_path):
    return pd.read_excel(file_path)

# 将残基序列连接成一个字符串，并确保数据清理
def concatenate_sequence(df):
    residues = df['ResidueInfo'].tolist()
    three_letter_residues = [residue.strip().upper() for residue in residues if residue and isinstance(residue, str)]
    one_letter_sequence = convert_three_to_one(three_letter_residues)

    if len(one_letter_sequence) != len(df):
        print(f"注意：有效残基数 ({len(one_letter_sequence)}) 与 DataFrame 行数 ({len(df)}) 不同。有 {len(df) - len(one_letter_sequence)} 行没有转换，将跳过这些行。")

    print(f"提取的序列: {one_letter_sequence} (长度: {len(one_letter_sequence)})")
    return one_letter_sequence, three_letter_residues

# 对每个氨基酸残基应用aaindex
def apply_aaindex(sequence, aaindex_dict, index_key):
    index_data = aaindex_dict.get(index_key, {}).get('index', {})
    if not index_data:
        raise ValueError(f"在 aaindex 数据中找不到索引 {index_key}。")

    feature_values = [index_data.get(residue, -999) for residue in sequence]
    return feature_values

# 主函数
def main(aaindex1_file, input_excel, output_excel, index_keys):
    aaindex_dict = parse_aaindex1(aaindex1_file)
    df = read_excel(input_excel)

    # 去除没有残基数的行
    df_clean = df.dropna(subset=['ResidueInfo']).copy()
    sequence, valid_residues = concatenate_sequence(df_clean)

    # 初始化结果DataFrame
    result_df = df_clean.copy()

    for index_key in index_keys:
        if index_key in excluded_indexes or index_key not in aaindex_dict:
            print(f"跳过索引：{index_key}")
            continue

        print(f"计算索引：{index_key}")
        try:
            feature_values = apply_aaindex(sequence, aaindex_dict, index_key)
            description = index_to_description.get(index_key, index_key)
            result_df[description] = feature_values
        except ValueError as ve:
            print(ve)

    result_df.to_excel(output_excel, index=False)
    print(f"结果已经成功保存到 {output_excel}")

# 使用示例
aaindex1_file = r'.\aaindex1.txt'
input_excel = r'.\protein.xlsx'
output_excel = r'.\aa.xlsx'

index_keys = ['LEVM780105', 'HOPT810101', 'KYTJ820101', 'MONM990101', 'CHOP780215', 'CHOP780216', 'GRAR740102', 'JANJ790102',
              'ARGP820103', 'DAYM780201', 'DESM900102', 'HUTJ700101', 'KLEP840101', 'KRIW790103', 'NAKH920106', 'NAKH920107',
              'NAKH920108', 'RACS770101', 'RACS770102', 'RACS770103', 'TAKK010101', 'FUKS010101', 'FUKS010102', 'FUKS010103',
              'FUKS010105', 'FUKS010106', 'FUKS010107', 'FUKS010109', 'FUKS010110', 'FUKS010111', 'COSI940101', 'ZHOH040101',
              'ZHOH040102', 'KARS160102', 'KARS160103', 'KARS160108', 'KARS160114', 'KARS160115', 'KARS160116', 'KARS160117',
              'KARS160118', 'KARS160105', 'MAXF760106', 'PONJ960101', 'CHAM820101', 'GRAR740103', 'FASG890101', 'GEIM800108',
              'HOPA770101', 'NAGK730103']

main(aaindex1_file, input_excel, output_excel, index_keys)
