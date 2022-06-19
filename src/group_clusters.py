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

genome_list = ["babesia_bigemina_protein", "babesia_bovis_T2Bo_protein",
              "babesia_microti_strain_RI_protein", "babesia_ovata_protein",
              "babesia_sp._Xinjiang_protein", "besnoitia_besnoiti_protein",
              "cryptosporidium_felis_protein", "cryptosporidium_hominis_TU502_protein",
              "cryptosporidium_muris_RN66_protein", "cryptosporidium_parvum_Iowa_II_protein",
              "cryptosporidium_ubiquitum_protein", "eimeria_acervulina_protein",
              "eimeria_brunetti_protein", "eimeria_maxima_protein",
              "eimeria_mitis_protein", "eimeria_necatrix_protein",
              "eimeria_praecox_protein", "eimeria_tenella_protein",
              "gregarina_niphandrodes_protein", "neospora_caninum_Liverpool_protein",
              "plasmodium_berghei_ANKA_protein", "plasmodium_chabaudi_protein",
              "plasmodium_coatneyi_protein", "plasmodium_cynomolgi_protein",
              "plasmodium_falciparum_protein", "plasmodium_fragile_protein",
              "plasmodium_gaboni_protein", "plasmodium_gallinaceum_protein",
              "plasmodium_inui_San_Antonio_1_protein", "plasmodium_knowlesi_protein",
              "plasmodium_malariae_protein", "plasmodium_reichenowi_protein",
              "plasmodium_relictum_protein", "plasmodium_sp._gorilla_clade_G2_protein",
              "plasmodium_vinckei_protein", "plasmodium_vivax_protein",
              "plasmodium_yoelii_protein", "porospora_cf._gigantea_A_protein",
              "porospora_cf._gigantea_B_protein", "theileria_annulata_protein",
              "theileria_equi_strain_WA_protein", "theileria_orientalis_strain_Shintoku_protein",
              "theileria_parva_strain_Muguga_protein", "toxoplasma_gondii_ME49_protein"]
'''
# Замена опечатки
clust = pd.read_csv("full_gca_cluster_stat_with_function.tsv", sep='\t')
clust.rename(columns = {"plasmosium_berghei_ANKA_protein.faa": "plasmodium_berghei_ANKA_protein.faa",
                        "plasmosium_chabaudi_protein.faa": "plasmodium_chabaudi_protein.faa",
                        "plasmosium_coatneyi_protein.faa": "plasmodium_coatneyi_protein.faa",
                        "plasmosium_cynomolgi_protein.faa": "plasmodium_cynomolgi_protein.faa",
                        "plasmosium_falciparum_protein.faa": "plasmodium_falciparum_protein.faa",
                        "plasmosium_fragile_protein.faa": "plasmodium_fragile_protein.faa",
                        "plasmosium_gaboni_protein.faa": "plasmodium_gaboni_protein.faa",
                        "plasmosium_gallinaceum_protein.faa": "plasmodium_gallinaceum_protein.faa",
                        "plasmosium_inui_San_Antonio_1_protein.faa": "plasmodium_inui_San_Antonio_1_protein.faa",
                        "plasmosium_knowlesi_protein.faa": "plasmodium_knowlesi_protein.faa",
                        "plasmosium_malariae_protein.faa": "plasmodium_malariae_protein.faa",
                        "plasmosium_reichenowi_protein.faa": "plasmodium_reichenowi_protein.faa",
                        "plasmosium_relictum_protein.faa": "plasmodium_relictum_protein.faa",
                        "plasmosium_sp._gorilla_clade_G2_protein.faa": "plasmodium_sp._gorilla_clade_G2_protein.faa",
                        "plasmosium_vinckei_protein.faa": "plasmodium_vinckei_protein.faa",
                        "plasmosium_vivax_protein.faa": "plasmodium_vivax_protein.faa",
                        "plasmosium_yoelii_protein.faa" : "plasmodium_yoelii_protein.faa"}, inplace = True)
clust.to_csv("full_gca_cluster_stat_with_function.tsv", sep='\t', index=False)
'''
'''
# Статистика по кластерам
clust = pd.read_csv("full_gca_44.tsv", sep='\t')
number_of_genus = []
number_of_genoms = []
number_of_proteins = []
for index, row in clust.iterrows():
    curr_gene_num = 0
    curr_prot_num = 0
    curr_genus_num = 0
    current_genus = ""
    for j in ["babesia_bigemina_protein.faa", "babesia_bovis_T2Bo_protein.faa",
              "babesia_microti_strain_RI_protein.faa", "babesia_ovata_protein.faa",
              "babesia_sp._Xinjiang_protein.faa", "besnoitia_besnoiti_protein.faa",
              "cryptosporidium_felis_protein.faa", "cryptosporidium_hominis_TU502_protein.faa",
              "cryptosporidium_muris_RN66_protein.faa", "cryptosporidium_parvum_Iowa_II_protein.faa",
              "cryptosporidium_ubiquitum_protein.faa", "eimeria_acervulina_protein.faa",
              "eimeria_brunetti_protein.faa", "eimeria_maxima_protein.faa",
              "eimeria_mitis_protein.faa", "eimeria_necatrix_protein.faa",
              "eimeria_praecox_protein.faa", "eimeria_tenella_protein.faa",
              "gregarina_niphandrodes_protein.faa", "neospora_caninum_Liverpool_protein.faa",
              "plasmodium_berghei_ANKA_protein.faa", "plasmodium_chabaudi_protein.faa",
              "plasmodium_coatneyi_protein.faa", "plasmodium_cynomolgi_protein.faa",
              "plasmodium_falciparum_protein.faa", "plasmodium_fragile_protein.faa",
              "plasmodium_gaboni_protein.faa", "plasmodium_gallinaceum_protein.faa",
              "plasmodium_inui_San_Antonio_1_protein.faa", "plasmodium_knowlesi_protein.faa",
              "plasmodium_malariae_protein.faa", "plasmodium_reichenowi_protein.faa",
              "plasmodium_relictum_protein.faa", "plasmodium_sp._gorilla_clade_G2_protein.faa",
              "plasmodium_vinckei_protein.faa", "plasmodium_vivax_protein.faa",
              "plasmodium_yoelii_protein.faa", "porospora_cf._gigantea_A_protein.faa",
              "porospora_cf._gigantea_B_protein.faa", "theileria_annulata_protein.faa",
              "theileria_equi_strain_WA_protein.faa", "theileria_orientalis_strain_Shintoku_protein.faa",
              "theileria_parva_strain_Muguga_protein.faa", "toxoplasma_gondii_ME49_protein.faa"]:
        if row[j] != '*':
            curr_gene_num += 1
            curr_prot_num += len(row[j].split(','))
            if j.split('_')[0] != current_genus:
                curr_genus_num += 1
                current_genus = j.split('_')[0]
    number_of_genus.append(curr_genus_num)
    number_of_genoms.append(curr_gene_num)
    number_of_proteins.append(curr_prot_num)
clust["type"] = ["full_gca_44"] * len(clust)
clust["number_of_genus"] = number_of_genus
clust["number_of_genoms"] = number_of_genoms
clust["number_of_proteins"] = number_of_proteins
clust.to_csv("full_gca_cluster_stat.tsv", sep='\t', index=False)
# clust.loc[(clust["number_of_genoms"] >= 22)].to_csv("full_gca_cluster_stat_22_plus.tsv", sep='\t', index=False)
'''
'''
# Гистограммы для кластеров
fig, ax = plt.subplots()
sns.histplot(data=clust, x="number_of_genoms", discrete=True, hue="type")
ax.locator_params(axis='x', integer=True)
plt.tight_layout()
plt.savefig("full_gca_44_clust_by_number_of_genoms.png", dpi=300)
fig, ax = plt.subplots()
sns.histplot(data=clust, x="number_of_proteins", discrete=True, hue="type")
ax.locator_params(axis='x', integer=True)
plt.tight_layout()
plt.savefig("full_gca_44_clust_by_number_of_proteins.png", dpi=300)
'''

'''
clust = pd.read_csv("full_gca_cluster_stat_22_plus.tsv", sep='\t')
fig, ax = plt.subplots()
sns.histplot(data=clust, x="number_of_proteins", discrete=True, hue="type")
ax.locator_params(axis='x', integer=True)
plt.tight_layout()
plt.savefig("full_gca_cluster_p_stat_22_plus.png", dpi=300)
'''
'''
# Определение функций белков
# clust.to_csv("glimpse.tsv", sep='\t', index=False)
all_prot_list = []
for j in ["babesia_bigemina_protein.faa", "babesia_bovis_T2Bo_protein.faa",
              "babesia_microti_strain_RI_protein.faa", "babesia_ovata_protein.faa",
              "babesia_sp._Xinjiang_protein.faa", "besnoitia_besnoiti_protein.faa",
              "cryptosporidium_felis_protein.faa", "cryptosporidium_hominis_TU502_protein.faa",
              "cryptosporidium_muris_RN66_protein.faa", "cryptosporidium_parvum_Iowa_II_protein.faa",
              "cryptosporidium_ubiquitum_protein.faa", "eimeria_acervulina_protein.faa",
              "eimeria_brunetti_protein.faa", "eimeria_maxima_protein.faa",
              "eimeria_mitis_protein.faa", "eimeria_necatrix_protein.faa",
              "eimeria_praecox_protein.faa", "eimeria_tenella_protein.faa",
              "gregarina_niphandrodes_protein.faa", "neospora_caninum_Liverpool_protein.faa",
              "plasmodium_berghei_ANKA_protein.faa", "plasmodium_chabaudi_protein.faa",
              "plasmodium_coatneyi_protein.faa", "plasmodium_cynomolgi_protein.faa",
              "plasmodium_falciparum_protein.faa", "plasmodium_fragile_protein.faa",
              "plasmodium_gaboni_protein.faa", "plasmodium_gallinaceum_protein.faa",
              "plasmodium_inui_San_Antonio_1_protein.faa", "plasmodium_knowlesi_protein.faa",
              "plasmodium_malariae_protein.faa", "plasmodium_reichenowi_protein.faa",
              "plasmodium_relictum_protein.faa", "plasmodium_sp._gorilla_clade_G2_protein.faa",
              "plasmodium_vinckei_protein.faa", "plasmodium_vivax_protein.faa",
              "plasmodium_yoelii_protein.faa", "porospora_cf._gigantea_A_protein.faa",
              "porospora_cf._gigantea_B_protein.faa", "theileria_annulata_protein.faa",
              "theileria_equi_strain_WA_protein.faa", "theileria_orientalis_strain_Shintoku_protein.faa",
              "theileria_parva_strain_Muguga_protein.faa", "toxoplasma_gondii_ME49_protein.faa"]:
    path_to_faa = "/home/leonid/PycharmProjects/python/minor_2022/project/faa_files/" + j
    with open(path_to_faa) as f:
        contents = f.readlines()
        for line in contents:
            if line[0] == '>':
                all_prot_list.append(line.rstrip())

function_list = []
for index, row in clust.iterrows():
    current_function_dict = {}
    for j in ["babesia_bigemina_protein.faa", "babesia_bovis_T2Bo_protein.faa",
              "babesia_microti_strain_RI_protein.faa", "babesia_ovata_protein.faa",
              "babesia_sp._Xinjiang_protein.faa", "besnoitia_besnoiti_protein.faa",
              "cryptosporidium_felis_protein.faa", "cryptosporidium_hominis_TU502_protein.faa",
              "cryptosporidium_muris_RN66_protein.faa", "cryptosporidium_parvum_Iowa_II_protein.faa",
              "cryptosporidium_ubiquitum_protein.faa", "eimeria_acervulina_protein.faa",
              "eimeria_brunetti_protein.faa", "eimeria_maxima_protein.faa",
              "eimeria_mitis_protein.faa", "eimeria_necatrix_protein.faa",
              "eimeria_praecox_protein.faa", "eimeria_tenella_protein.faa",
              "gregarina_niphandrodes_protein.faa", "neospora_caninum_Liverpool_protein.faa",
              "plasmodium_berghei_ANKA_protein.faa", "plasmodium_chabaudi_protein.faa",
              "plasmodium_coatneyi_protein.faa", "plasmodium_cynomolgi_protein.faa",
              "plasmodium_falciparum_protein.faa", "plasmodium_fragile_protein.faa",
              "plasmodium_gaboni_protein.faa", "plasmodium_gallinaceum_protein.faa",
              "plasmodium_inui_San_Antonio_1_protein.faa", "plasmodium_knowlesi_protein.faa",
              "plasmodium_malariae_protein.faa", "plasmodium_reichenowi_protein.faa",
              "plasmodium_relictum_protein.faa", "plasmodium_sp._gorilla_clade_G2_protein.faa",
              "plasmodium_vinckei_protein.faa", "plasmodium_vivax_protein.faa",
              "plasmodium_yoelii_protein.faa", "porospora_cf._gigantea_A_protein.faa",
              "porospora_cf._gigantea_B_protein.faa", "theileria_annulata_protein.faa",
              "theileria_equi_strain_WA_protein.faa", "theileria_orientalis_strain_Shintoku_protein.faa",
              "theileria_parva_strain_Muguga_protein.faa", "toxoplasma_gondii_ME49_protein.faa"]:
        if row[j] != '*':
            for protein in row[j].split(','):
                function = ""
                for protein_string in all_prot_list:
                    if protein_string.split(" ", 1)[0][1:] == protein:
                        function = protein_string.split(" ", 1)[1].rsplit(" [", 1)[0]
                if function in current_function_dict:
                    current_function_dict[function] += 1
                else:
                    if function != "":
                        current_function_dict[function] = 1
    function_list.append(max(current_function_dict, key=current_function_dict.get))
clust["function"] = function_list
clust.to_csv("full_gca_cluster_stat_with_function.tsv", sep='\t', index=False)
print(clust)
'''
'''
clust = pd.read_csv("full_gca_cluster_stat_with_function.tsv", sep='\t')
z_data = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_z_data = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "_intersection.bed"
    z_data[i] = pd.read_csv(path_to_z_data, sep='\t', names=["chromosome", "start", "end", "ID", "chromosome_2",
                                                             "start_2", "end_2", "chr+start", "Z-SCORE"])
    z_data[i] = z_data[i].drop(["chromosome_2", "start_2", "end_2", "chr+start"], axis=1).sort_values(by=["Z-SCORE"],
                                                                                                      ascending=False)

score_dict = {"plasmodium_falciparum_Z-SCORE": [],
              "plasmodium_gaboni_Z-SCORE": [],
              "plasmodium_knowlesi_Z-SCORE": [],
              "plasmodium_vivax_Z-SCORE": [],
              "plasmodium_yoelii_Z-SCORE": []}
for index, row in clust.iterrows():
    current_score = {"falciparum": 0, "gaboni": 0, "knowlesi": 0, "vivax": 0, "yoelii": 0}
    for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        label_name = "plasmodium_" + i + "_protein.faa"
        if row[label_name] == '*':
            score_dict[("plasmodium_" + i + "_Z-SCORE")].append('*')
        else:
            for this_ID in row[label_name].split(','):
                number_of_ID = len(row[label_name].split(','))
                if this_ID in list(z_data[i]["ID"]):
                    z_index = z_data[i].index[z_data[i]["ID"] == this_ID].tolist()[0]
                    current_score[i] += z_data[i]["Z-SCORE"][z_index]
                else:
                    if number_of_ID >= 2:
                        number_of_ID -= 1
            score_dict[("plasmodium_" + i + "_Z-SCORE")].append(current_score[i] / number_of_ID)
clust = clust.join(pd.DataFrame(score_dict))
clust.to_csv("clusters_with_z_score.tsv", sep='\t', index=False)
print(clust.loc[(clust["plasmodium_yoelii_Z-SCORE"] != 0) & (clust["plasmodium_yoelii_Z-SCORE"] != '*')])
'''
'''
clust = pd.read_csv("clusters_with_z_score.tsv", sep='\t')
clust.rename(columns = {"plasmodium_falciparum_Z-SCORE": "plasmodium_falciparum",
                        "plasmodium_vivax_Z-SCORE": "plasmodium_vivax",
                        "plasmodium_gaboni_Z-SCORE": "plasmodium_gaboni",
                        "plasmodium_knowlesi_Z-SCORE": "plasmodium_knowlesi",
                        "plasmodium_yoelii_Z-SCORE": "plasmodium_yoelii"}, inplace = True)
clust.to_csv("clusters_with_z_score.tsv", sep='\t', index=False)
# print(z_data)
print(clust)
'''