import pandas as pd
import ast
import numpy as np
import csv
import matplotlib.pyplot as plt
import datetime
import numpy as np

class Hic_proj_GWAS():

    def __init__(self):
        self.promoter_overlaps_bin = []
        self.counter = 0
        self.promoters_bins = []
        self.non_promoters_bins = []
        self.both_bins = []
        self.bin_dict = []
        self.GWAS_seprated_by_chr = []
        self.element_data_seprated_by_chr = []
        self.GWAS_element_dict = {}
        self.element_GWAS_dict = {}
        self.results = []
        self.all = []
        self.temp = set()

    def seprate_schizophrni_GWAS(self):
        GWAS_file_name = "GWAS/GWAS_NDD_MINE.xlsx"
        data = pd.read_excel(GWAS_file_name, header=None)
        data.columns = list(data.loc[0])
        # print(list(data.loc[0]))
        # print(data)
        Schizophrenia_data = data.query("HPO_TERM == 'Schizophrenia'")
        not_Schizophrenia_data = data.query("HPO_TERM != 'Schizophrenia'")
        Schizophrenia_data = Schizophrenia_data.reset_index(drop=True)
        not_Schizophrenia_data = not_Schizophrenia_data.reset_index(drop=True)

        # print(len(data))
        # print(len(Schizophrenia_data) , len(not_Schizophrenia_data))
        start = data.loc[:]["start_hg19"]
        end = data.loc[:]["end_hg19"]
        chr = data.loc[:]["chr"]
        # new_df = list(Schizophrenia_data.loc[:]["chr"]) + list(Schizophrenia_data.loc[:]["start_hg19"])
        # print(new_df)
        # new_df.to_csv("GWAS/Schizophrenia_GWAS.csv",index=False)
        Schizophrenia_data = data[data['HPO_TERM'].str.contains("Schizophrenia")]
        Schizophrenia_data = Schizophrenia_data.reset_index(drop=True)
        not_Schizophrenia_data = data[~data['HPO_TERM'].str.contains("Schizophrenia")]
        not_Schizophrenia_data = not_Schizophrenia_data.reset_index(drop=True)

        chr = ["chr"] + list(Schizophrenia_data.loc[:]["chr"])
        s = ["pos"] + list(Schizophrenia_data.loc[:]["start_hg19"])
        record_list = [list(item) for item in
                        list(zip(chr, s))]
        with open("GWAS/Schizophrenia_GWAS2.csv", "w",newline='') as fp:
            writer = csv.writer(fp)
            writer.writerows(record_list)

        chr = ["chr"] + list(not_Schizophrenia_data.loc[1:]["chr"])
        s = ["pos"] + list(not_Schizophrenia_data.loc[1:]["start_hg19"])
        record_list = [list(item) for item in
                       list(zip(chr, s))]

        with open("GWAS/not_Schizophrenia_GWAS2.csv", "w", newline='') as fp:
            writer = csv.writer(fp)
            writer.writerows(record_list)

    def seprate_autism_GWAS(self):
        GWAS_data = pd.read_excel("GWAS/GWAS_NDD_MINE.xlsx", header=None)
        GWAS_data.columns = list(GWAS_data.loc[0])
        # print(list(data.loc[0]))
        # print(data)
        autism_data = GWAS_data[GWAS_data['HPO_TERM'].str.contains("Autism|autism|Autistic|autistic")]
        chr = ["chr"] + list(autism_data.loc[1:]["chr"])
        s = ["pos"] + list(autism_data.loc[1:]["start_hg19"])
        record_list = [list(item) for item in
                       list(zip(chr, s))]

        with open("GWAS/autism_data.csv", "w", newline='') as fp:
            writer = csv.writer(fp)
            writer.writerows(record_list)

    def make_bin_dictionary_interactions(self,file_name):

        bins_dicts = {}
        unique = []
        temp = []
        interaction_file_name = file_name
        interaction_data = pd.read_csv(interaction_file_name, header=None)
        interaction_data.columns = ['INT_ID', 'F1_ID', 'F1_chr', 'F1_start', 'F1_end', 'F2_ID', 'F2_chr', 'F2_start',
                                    'F2_end', 'distance', 'interacting_read', 'logpvalue']
        interaction_data = interaction_data.astype("str")

        F1_bin = list(interaction_data.loc[:]["F1_ID"])
        F1_start = list(interaction_data.loc[:]["F1_start"])
        F1_end = list(interaction_data.loc[:]["F1_end"])
        F1_chr = list(interaction_data.loc[:]["F1_chr"])

        F2_bin = list(interaction_data.loc[:]["F2_ID"])
        F2_start = list(interaction_data.loc[:]["F2_start"])
        F2_end = list(interaction_data.loc[:]["F2_end"])
        F2_chr = list(interaction_data.loc[:]["F2_chr"])

        for n, bin in enumerate(F1_bin):
            temp.append(bin)
            if n == 0:
                continue
            if bin not in bins_dicts:
                unique.append(
                    {"id": int(F1_bin[n]),
                     "start": int(F1_start[n]),
                     "end": int(F1_end[n]), "chr": int(F1_chr[n])})
                if int(F1_chr[n]) < 24 :
                    bins_dicts[bin] = [int(F2_bin[n])]
            else:
                if int(F1_chr[n]) < 24:
                    bins_dicts[bin].append(int(F2_bin[n]))

        for n, bin in enumerate(F2_bin):
            temp.append(bin)
            if n == 0:
                continue
            if bin not in bins_dicts:
                unique.append(
                    {"id": int(F2_bin[n]),
                     "start": int(F2_start[n]),
                     "end": int(F2_end[n]), "chr": int(F2_chr[n])})
                if int(F1_chr[n]) < 24:
                    bins_dicts[bin] = [int(F1_bin[n])]
            else:
                if int(F1_chr[n]) < 24:
                    bins_dicts[bin].append(int(F1_bin[n]))

        # print(len(bins_dicts))
        return bins_dicts

    def call_promoters_related_bins(self,file_name):

        promoter_overlap_file_name = file_name
        promoter_overlap_data = pd.read_csv(promoter_overlap_file_name, header=None, low_memory=False)
        promoter_overlap_data.columns = list(promoter_overlap_data.loc[0])
        print(list(promoter_overlap_data.loc[0]))
        overlaped_bins = promoter_overlap_data.loc[1:]["overlaped_element_F_ID"]
        v = np.vectorize(self.promoters_related_bins)
        v(overlaped_bins)
        self.promoter_overlaps_bin = set(self.promoter_overlaps_bin)
        # print("here : ")
        # print(len(self.promoter_overlaps_bin))
        # print("=========================")

    def promoters_related_bins(self, bin):
        self.counter += 1
        x = ast.literal_eval(bin)
        for i in x:
            self.promoter_overlaps_bin.append(int(i))

    def call_seprate_related_promoters_bins(self,file_name,promoter_file_name):

        self.bin_dict = self.make_bin_dictionary_interactions(file_name)
        self.call_promoters_related_bins(promoter_file_name)
        self.counter = 0
        v = np.vectorize(self.seprate_related_promoters_bins)
        v(list(self.bin_dict))


    def seprate_related_promoters_bins(self, bin):

        self.counter += 1
        promoter_related = 0

        for i in self.bin_dict[bin]:
            if i in self.promoter_overlaps_bin:
                promoter_related += 1
        if promoter_related == len(self.bin_dict[bin]):
            self.promoters_bins.append(int(bin))
        elif promoter_related < len(self.bin_dict[bin]) and promoter_related != 0:
            self.both_bins.append(int(bin))
        elif promoter_related == 0:
            self.non_promoters_bins.append(int(bin))


    # mode = 0 for GWAS else for others
    def seprate_chr(self, data, mode=0):
        seprated_by_chr = []
        if mode == 0:
            for chr in range(1, 24):
                if chr < 23:
                    chro = "chr" + str(chr)
                else:
                    chro = "chrX"
                temp = data.query("chr == @chro")
                temp = temp.sort_values(by=['pos'])
                temp = temp.reset_index(drop=True)
                seprated_by_chr.append(temp)


        if mode == 1:
            for chr in range(1, 24):
                # chro = "chr" + str(chr)
                chr_str = str(chr)
                temp = data.query("chr == @chr or chr == @chr_str")
                temp = temp.sort_values(by=['start'])
                temp = temp.reset_index(drop=True)
                seprated_by_chr.append(temp)
            s = 0
        return seprated_by_chr

    def GWAS_annotation(self,file_name, data ,element_file_name, permutation = False):
        self.element_GWAS_dict = {}
        if not permutation :
            path = "GWAS/"
            if file_name == "breast_cancer_GWAS" :
                GWAS_file_name = path + file_name + ".tsv"
                GWAS_data = pd.read_csv(GWAS_file_name, sep="\t", header=None)

            else:
                GWAS_file_name = path + file_name + ".csv"
                GWAS_data = pd.read_csv(GWAS_file_name, header=None)

        if permutation:
            GWAS_data = data
        GWAS_data.columns = ["chr","pos"]
        GWAS_data = GWAS_data.loc[1:]
        GWAS_data = GWAS_data.astype({'pos': 'int64'})
        element_data = pd.read_csv(element_file_name, header=None,low_memory=False)
        element_data.columns = list(element_data.loc[0])
        # print(list(element_data.loc[0]))

        element_data = element_data.loc[1:]
        element_data = element_data.astype({'start': 'int64'})
        element_data = element_data.astype({'end': 'int64'})
        GWAS_chr = GWAS_data["chr"].tolist()
        GWAS_pos = GWAS_data["pos"].tolist()

        self.GWAS_seprated_by_chr = self.seprate_chr(GWAS_data)
        self.element_data_seprated_by_chr = self.seprate_chr(element_data, mode=1)
        t1 = datetime.datetime.now()
        w = np.vectorize(self.check_GWAS_element)
        w(GWAS_chr, GWAS_pos)
        t2 = datetime.datetime.now()
        # print(len(self.element_GWAS_dict))
        # print(t2 - t1)
        self.all.append(self.element_GWAS_dict)
        self.map_GWAS_to_bins(file_name, element_file_name)

    def check_element_GWAS(self, element_id, element_start, element_end, element_chr):

        t1 = datetime.datetime.now()
        if element_chr == "X" or element_chr == "Y":
            element_chr = 23
        element_chr = int(element_chr)
        mutation_pos = list(self.GWAS_seprated_by_chr[element_chr - 1].loc[:]["pos"])
        start_index = np.searchsorted(mutation_pos,element_start)
        end_index = np.searchsorted(mutation_pos,element_end)

        for index in range(int(start_index), end_index):

            # if mutation_pos[index] >= element_end :
            #     break

            if mutation_pos[index] > element_start and mutation_pos[index] < element_end:
                self.temp.add(int(element_id))
                # if element_id not in set(self.element_GWAS_dict):
                #     self.element_GWAS_dict[int(element_id)] = [mutation_pos[index]]
                # else:
                #     self.element_GWAS_dict[int(element_id)].append(mutation_pos[index])
        # print(m)


    def check_GWAS_element(self, GWAS_chr, GWAS_pos):

        if GWAS_chr == "chrX" or GWAS_chr == "chrY":
            GWAS_chr = "chr23"

        GWAS_chr = int(GWAS_chr.replace("chr", ""))
        element_data = self.element_data_seprated_by_chr[GWAS_chr - 1]

        element_id = element_data['id'].tolist()
        element_start = element_data['start'].tolist()
        element_end = element_data['end'].tolist()

        start_index = np.searchsorted(element_end, GWAS_pos)
        end_index = np.searchsorted(element_start, GWAS_pos)

        for index in range(start_index, end_index):

            if GWAS_pos > element_start[index] and GWAS_pos < element_end[index]:
                # if GWAS_pos not in self.GWAS_element_dict:
                #     self.GWAS_element_dict[GWAS_pos] = [element_id[index]]
                # else:
                #     self.GWAS_element_dict[GWAS_pos].append(element_id[index])

                if element_id[index] not in self.element_GWAS_dict:
                    self.element_GWAS_dict[int(element_id[index])] = [GWAS_pos]
                else:
                    self.element_GWAS_dict[int(element_id[index])].append(GWAS_pos)

    def make_dict_of_bins_and_element_overlap(self, file_name):
        dict = {}
        element_data = pd.read_csv(
            "enhancer/bin_id_interaction_overlap_" + file_name + "_file" + ".csv", header=None)
        element_data.columns = ['bin_id', 'element_id']

        element_bins = element_data["bin_id"]
        element_ids = element_data["element_id"]

        for i, j in zip(element_bins, element_ids):
            dict[int(i)] = ast.literal_eval(j)

        return dict

    def map_GWAS_to_bins(self,f,file_name):
        element_data = pd.read_csv(file_name)
        element_data.columns = list(element_data)
        element_promoter_bins = []
        element_non_promoter_bins = []
        element_both_bins = []
        element_data.astype({'id': 'int64'})
        element_bins = list(set(element_data.loc[:]["id"]))
        n = 0

        s_promoters_bin = set(self.promoters_bins)

        s_non_promoters_bin = set(self.non_promoters_bins)

        s_both_bin = set(self.both_bins)

        promoter_GWAS = 0
        non_promoter_GWAS = 0
        both_GWAS = 0
        for bin in element_bins:

            if bin in s_promoters_bin:
                element_promoter_bins.append(bin)
                x = self.element_GWAS_dict
                if bin in self.element_GWAS_dict:
                    promoter_GWAS += 1

            elif bin in s_non_promoters_bin:
                element_non_promoter_bins.append(bin)
                if bin in self.element_GWAS_dict:
                    non_promoter_GWAS += 1

            elif bin in s_both_bin:
                element_both_bins.append(bin)
                if bin in self.element_GWAS_dict:
                    both_GWAS += 1

        n = non_promoter_GWAS / len(element_non_promoter_bins)
        p = (promoter_GWAS + both_GWAS) / (len(element_promoter_bins) + len(element_both_bins))
        all = (promoter_GWAS + non_promoter_GWAS + both_GWAS) / (
              len(element_promoter_bins) + len(element_non_promoter_bins) + len(element_both_bins))
        # print(n)
        # print(p)
        # print(all)
        # self.draw_bar_chart(p, n, all, "CNON " + f, f)

        # print(len(element_promoter_bins))

        # print(len(element_promoter_bins) - promoter_GWAS)
        # print(len(element_non_promoter_bins) - non_promoter_GWAS)
        # print(len(element_both_bins) - both_GWAS)
        #
        # print("only promoter interactions : ", promoter_GWAS / len(element_promoter_bins))
        # print(non_promoter_GWAS)
        # print(promoter_GWAS)
        # print(both_GWAS)
        # print(len(element_non_promoter_bins))
        # print(len(element_promoter_bins))
        # print(len(element_both_bins))
        #
        print("only non-promoter interactions : ", non_promoter_GWAS / len(element_non_promoter_bins))

        print("only promoter interactions and both : ",
              (promoter_GWAS + both_GWAS) / (len(element_promoter_bins) + len(element_both_bins)))

        print("all : ", (promoter_GWAS + non_promoter_GWAS + both_GWAS),
              len(element_promoter_bins) + len(element_non_promoter_bins) + len(element_both_bins),
              (promoter_GWAS + non_promoter_GWAS + both_GWAS) / (
              len(element_promoter_bins) + len(element_non_promoter_bins) + len(element_both_bins)))

        self.results.append(p)
        return  p

    def draw_bar_chart(self, m, n, k, name,GWAS_type):

        objects = ("Bin that \n  have interactions \n with promoters",
                   "Bin that \n  have interactions \n with non-promoters",
                   "All bins ")
        y_pos = np.arange(3)
        performance = [m * 100, n * 100, k * 100]
        bar_list = plt.bar(y_pos, performance, align='center', alpha=0.9, width=0.2)
        bar_list[0].set_color('lightgreen')
        bar_list[1].set_color('darkseagreen')
        bar_list[2].set_color('lawngreen')
        plt.xticks(y_pos, objects)
        plt.ylabel('Percentage of bins which overlap with '+GWAS_type)
        plt.title(name)
        plt.savefig("GWAS/chart" + "/" + name + ".png")
        plt.savefig(name + ".png")


p = Hic_proj_GWAS()
p.call_seprate_related_promoters_bins("Neun-/Neun-_Neun_HiC.csv","promoter/Neun-_overlap_promoter.csv")
for i in ["breast_cancer_GWAS" , "Schizophrenia_GWAS" , "not_Schizophrenia_GWAS", "autism_GWAS" ] :
    p.GWAS_annotation(i ,"","bin_data/Neun-_new_unique_bin.csv",False)
    print("finish", i)
#
# p.draw_bar_chart(0.0016, 0.0014, 0.0015, "CNON " + "autism GWAS", "autism GWASS")
