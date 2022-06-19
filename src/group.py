import os
import csv
import sys
import tempfile
import fileinput
import subprocess
from Bio import SeqIO
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from subprocess import DEVNULL, PIPE
from dna_features_viewer import GraphicFeature, GraphicRecord

import pandas as pd

ZH_EXECUTABLE = Path("/home/leonid/zhunt3")
assert ZH_EXECUTABLE.is_file()

# Функция, предсказывающая участки Z-ДНК
def zhunt(query: str, windowsize: int = 6, minsize: int = 3, maxsize: int = 6):
    assert set(query).issubset({"A", "C", "G", "T", "N"})
    fd, temp = tempfile.mkstemp()
    os.close(fd)
    with open(temp, 'w') as stream:
        stream.write(query)

    subprocess.run(
        [ZH_EXECUTABLE,
         str(windowsize), str(minsize), str(maxsize), temp],
        check=True, stdout=PIPE, stderr=DEVNULL,
        input=query, encoding='ascii'
    )
    with open(temp + ".Z-SCORE", 'r') as stream:
        df = pd.read_csv(stream,
                         names=['Start', 'End', 'nu-1', 'nu-2', 'nu-3',
                                'ZH-Score', 'Sequence', 'Conformation'],
                         skiprows=1, sep='\s+')
    os.remove(temp)
    os.remove(temp + ".Z-SCORE")
    return df[['Start', 'End', 'ZH-Score', 'Sequence', 'Conformation']]

# Функция, создающая .bed файл
def make_bed(data, fileName, genome_file):
    to_write = data.copy()
    to_write.insert(0, 'chr', str(genome_file).split('.')[0])
    to_write['name'] = to_write.chr.astype(str) + '_' + to_write.Start.astype(str)
    to_write[['chr', 'Start', 'End', 'name', 'ZH-Score']].to_csv(fileName, sep='\t', header=False, index=False)

genome_list = ["babesia_bigemina", "babesia_bovis_T2Bo",
              "babesia_microti_strain_RI", "babesia_ovata",
              "babesia_sp._Xinjiang", "besnoitia_besnoiti",
              "cryptosporidium_felis", "cryptosporidium_hominis_TU502",
              "cryptosporidium_muris_RN66", "cryptosporidium_parvum_Iowa_II",
              "cryptosporidium_ubiquitum", "eimeria_acervulina",
              "eimeria_brunetti", "eimeria_maxima",
              "eimeria_mitis", "eimeria_necatrix",
              "eimeria_praecox", "eimeria_tenella",
              "gregarina_niphandrodes", "neospora_caninum_Liverpool",
              "plasmodium_berghei_ANKA", "plasmodium_chabaudi",
              "plasmodium_coatneyi", "plasmodium_cynomolgi",
              "plasmodium_falciparum", "plasmodium_fragile",
              "plasmodium_gaboni", "plasmodium_gallinaceum",
              "plasmodium_inui_San_Antonio_1", "plasmodium_knowlesi",
              "plasmodium_malariae", "plasmodium_reichenowi",
              "plasmodium_relictum", "plasmodium_sp._gorilla_clade_G2",
              "plasmodium_vinckei", "plasmodium_vivax",
              "plasmodium_yoelii", "porospora_cf._gigantea_A",
              "porospora_cf._gigantea_B", "theileria_annulata",
              "theileria_equi_strain_WA", "theileria_orientalis_strain_Shintoku",
              "theileria_parva_strain_Muguga", "toxoplasma_gondii_ME49"]

'''
# Этот блок заменяет все символы во .fna файлах на A, C, G, T, N
for i in genome_list:
    path_to_fna = "/home/leonid/PycharmProjects/python/minor_2022/project/supergroup/" + i + "_genomic.fna"
    for line in fileinput.input(path_to_fna, inplace=1):
        if line[0] == '>':
            new_line = line
        else:
            new_line = ''
            for i in line.rstrip():
                if i in ["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"]:
                    new_line += i
                else:
                    new_line += "N"
            new_line = new_line.upper()
            new_line += '\n'
        sys.stdout.write(new_line)
'''
'''
# Этот блок извлекает координаты промоторов из файлов feature table
feature_tables = {}
for i in genome_list:
    path_to_ft = "/home/leonid/PycharmProjects/python/minor_2022/project/supergroup/" + i + "_feature_table.txt"
    feature_tables[i] = pd.read_csv(path_to_ft, sep='\t', dtype=str)
tss_tables = {}
for i in genome_list:
    chromosome_list = []
    start_list = []
    end_list = []
    ID_list = []
    for index, row in feature_tables[i].iterrows():
        if row['# feature'] == 'CDS':
            chromosome_list.append(row["genomic_accession"])
            if type(row["product_accession"]) is not str:
                ID_list.append("Unknown")
            else:
                ID_list.append(row["product_accession"])
            if row['strand'] == '+':
                if (int(row["start"]) - 300) > 0:
                    start_list.append(int(row["start"]) - 300)
                    end_list.append(int(row["start"]))
                else:
                    start_list.append((int(row["start"])))
                    end_list.append(int(row["start"]) + 300)
            else:
                if (int(row["end"]) - 300) > 0:
                    start_list.append(int(row["end"]))
                    end_list.append(int(row["end"]) + 300)
                else:
                    start_list.append(int(row["end"]))
                    end_list.append(int(row["end"]) + 300)
    tss_tables[i] = pd.DataFrame({'chromosome': chromosome_list, 'start': start_list, 'end': end_list, 'ID': ID_list})
    print(tss_tables)
for i in genome_list:
    tss_tables[i].to_csv((i + "_tss" ".bed"), sep='\t', index=False, header=None)
'''
'''
# Этот блок создаёт .bed файлы, содержащие Z-SCORE для 44 геномов
# При этом zhunt запускается для каждого сиквенса во .fna файле отдельно
chr_dict = {}
for i in ["plasmodium_yoelii", "porospora_cf._gigantea_A",
              "porospora_cf._gigantea_B", "theileria_annulata",
              "theileria_equi_strain_WA", "theileria_orientalis_strain_Shintoku",
              "theileria_parva_strain_Muguga", "toxoplasma_gondii_ME49"]:
    path_to_fna = "/home/leonid/PycharmProjects/python/minor_2022/project/supergroup/fna_files/" + i + "_genomic.fna"
    with open(path_to_fna) as f:
        gene_lines = f.readlines()
    start_reading = True
    for j in gene_lines:
        if j[0] == '>':
            if i in chr_dict:
                chr_dict[i].append(j.split(" ")[0][1:])
            else:
                chr_dict[i] = []
                chr_dict[i].append(j.split(" ")[0][1:])
            if start_reading is True:
                start_reading = False
            else:
                fna_file.close()
            fna_file = open((j.split(" ")[0][1:] + ".fna"), "x")
            fna_file.write(j)
        else:
            fna_file.write(j)
    fna_file.close()
for key in chr_dict:
    general_bed_file = open((key + "_z.bed"), "x")
    for curr_chr in chr_dict[key]:
        name_for_fna = curr_chr + ".fna"
        with open(name_for_fna, 'r') as stream:
            sequences = list(SeqIO.parse(stream, format='fasta'))

        # for s in sequences:
        #    print(f"{s.description} \n\t=> {len(s)}")

        seq = sequences[0]
        df = zhunt(str(seq.seq))
        z = df[df['ZH-Score'] >= 500]
        name_for_bed = curr_chr + "_z.bed"
        make_bed(z, name_for_bed, name_for_fna)
        os.remove(name_for_fna)
        with open(name_for_bed, "r") as curr_bed_file:
            list_of_lines = curr_bed_file.readlines()
            for common_line in list_of_lines:
                general_bed_file.write(common_line)
        os.remove(name_for_bed)
    general_bed_file.close()
'''
'''
# Этот блок дописывает потерянные символы для genomic_accession
for i in genome_list:
    name_for_bed = i + "_z.bed"
    curr_df = pd.read_csv(name_for_bed, names=["chromosome", "start", "end", "chr+start", "Z-SCORE"], sep='\t')
    new_chromosome = []
    for index, row in curr_df.iterrows():
        if i == "plasmodium_berghei_ANKA":
            if row["chromosome"] == "CABFNT010000001" or row["chromosome"] == "CABFNT010000002" \
                    or row["chromosome"] == "CABFNT010000003" or row["chromosome"] == "LK023131":
                new_chromosome.append(row["chromosome"] + ".1")
            else:
                new_chromosome.append(row["chromosome"] + ".2")
        elif i == "plasmodium_falciparum":
            if row["chromosome"] == "LN999943" or row["chromosome"] == "LN999944" \
                    or row["chromosome"] == "LN999945" or row["chromosome"] == "LN999947" \
                    or row["chromosome"] == "LN999946" or row["chromosome"] == "LR605957" \
                    or row["chromosome"] == "LR605956":
                new_chromosome.append(row["chromosome"] + ".1")
            elif row["chromosome"] == "AL844506" or row["chromosome"] == "AL844507" \
                    or row["chromosome"] == "AL844509":
                new_chromosome.append(row["chromosome"] + ".3")
            else:
                new_chromosome.append(row["chromosome"] + ".2")
        elif i == "plasmodium_knowlesi":
            if row["chromosome"][0:4] == "CAD" or row["chromosome"] == "LR761316" \
                    or row["chromosome"] == "LR761315":
                new_chromosome.append(row["chromosome"] + ".1")
            elif row["chromosome"] == "AM910988" or row["chromosome"] == "AM910996":
                new_chromosome.append(row["chromosome"] + ".3")
            else:
                new_chromosome.append(row["chromosome"] + ".2")
        elif i == "plasmodium_yoelii":
            new_chromosome.append(row["chromosome"] + ".2")
        elif i == "theileria_annulata":
            if row["chromosome"] == "CR940352":
                new_chromosome.append(row["chromosome"] + ".2")
            else:
                new_chromosome.append(row["chromosome"] + ".1")
        else:
            new_chromosome.append(row["chromosome"] + ".1")
    curr_df = curr_df.drop(["chromosome"], axis=1)
    curr_df.insert(loc=0, column='chromosome', value=new_chromosome)
    curr_df.to_csv(name_for_bed, sep='\t', index=False, header=None)
'''
# На этом месте запускается bash скрипт, объединяющий участки Z-ДНК
'''
for i in genome_list:
    name_for_bed = "final_" + i + "_z.bed"
    curr_df = pd.read_csv(name_for_bed, names=["chromosome", "start", "end", "chr+start", "Z-SCORE"], sep=' ')
    curr_df.to_csv(name_for_bed, sep='\t', index=False, header=None)
'''
# На этом месте запускается bash скрипт, пересекающий tss и final_z файлы
'''
# Дополняем таблицу Z-SCORE
z_data = {}
for i in genome_list:
    path_to_z_data = "/home/leonid/PycharmProjects/python/minor_2022/project/supergroup/" + i + "_intersection.bed"
    z_data[i] = pd.read_csv(path_to_z_data, sep='\t', names=["chromosome", "start", "end", "ID", "chromosome_2",
                                                             "start_2", "end_2", "chr+start", "Z-SCORE"])
    z_data[i] = z_data[i].drop(["chromosome_2", "start_2", "end_2", "chr+start"], axis=1)
score_dict = {}
for i in genome_list:
    score_dict[i] = []
clust = pd.read_csv("full_gca_cluster_stat_with_function.tsv", sep='\t')
# print(z_data)
# print(clust)

for index, row in clust.iterrows():
    current_score = {}
    for i in genome_list:
        current_score[i] = 0
    for i in genome_list:
        label_name = i + "_protein.faa"
        if row[label_name] == '*':
            score_dict[i].append('*')
        else:
            for this_ID in row[label_name].split(','):
                number_of_ID = len(row[label_name].split(','))
                if this_ID in list(z_data[i]["ID"]):
                    z_index = z_data[i].index[z_data[i]["ID"] == this_ID].tolist()[0]
                    current_score[i] += z_data[i]["Z-SCORE"][z_index]
                else:
                    if number_of_ID >= 2:
                        number_of_ID -= 1
            score_dict[i].append(current_score[i] / number_of_ID)
clust = clust.join(pd.DataFrame(score_dict))
clust.to_csv("clusters_with_z_score.tsv", sep='\t', index=False)
'''
'''
# Добавляем колонку с суммарным Z-SCORE
clust = pd.read_csv("clusters_with_z_score.tsv", sep='\t')
z_sum = []
for index, row in clust.iterrows():
    current_sum = 0
    for i in genome_list:
        if row[i] != '*':
            current_sum += float(row[i])
    z_sum.append(current_sum)
clust["z_sum"] = z_sum
clust.sort_values(by=["z_sum"], ascending=False).to_csv("clusters_with_z_score.tsv", sep='\t', index=False)
'''

clust = pd.read_csv("clusters_with_z_score.tsv", sep='\t')