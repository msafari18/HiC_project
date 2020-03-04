import pandas as pd
import ast
import numpy as np
import csv
import matplotlib.pyplot as plt


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

    def seprate_schizophrni_GWAS(self):
        GWAS_file_name = "GWAS/GWAS_NDD_MINE.xlsx"
        data = pd.read_excel(GWAS_file_name, header=None)
        data.columns = list(data.loc[0])
        # print(list(data.loc[0]))
        print(data)
        Schizophrenia_data = data.query("HPO_TERM == 'Schizophrenia'")
        not_Schizophrenia_data = data.query("HPO_TERM != 'Schizophrenia'")
        Schizophrenia_data = Schizophrenia_data.reset_index(drop=True)
        not_Schizophrenia_data = not_Schizophrenia_data.reset_index(drop=True)

        print(len(data))
        print(len(Schizophrenia_data) , len(not_Schizophrenia_data))
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


    def make_bin_dictionary_interactions(self):

        bins_dicts = {}
        unique = []
        temp = []
        interaction_file_name = "HMEC_SRR1658680_5k_0.01_5reads_2.05.19.txt"
        interaction_data = pd.read_csv(interaction_file_name, sep="\t", header=None)
        interaction_data.columns = ['INT_ID', 'F1_ID', 'F1_chr', 'F1_start', 'F1_end', 'F2_ID', 'F2_chr', 'F2_start',
                                    'F2_end', 'distance', 'interacting_read', 'logpvalue']
        # print(list(interaction_data.loc[0]))
        # interaction_data = interaction_data.loc[1:]
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
                bins_dicts[bin] = [int(F2_bin[n])]
            else:
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
                bins_dicts[bin] = [int(F1_bin[n])]
            else:
                bins_dicts[bin].append(int(F1_bin[n]))
        # unique = sorted(unique, key=lambda i: i['id'])
        # csv_columns = ["id","start","end","chr"]
        # try:
        #     with open("bin_data/new_unique_bin.csv", 'w', newline='') as csvfile:
        #         writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        #         writer.writeheader()
        #         for data in unique:
        #             writer.writerow(data)
        # except IOError:
        #     print("I/O error / can not write on csv file")
        return bins_dicts

    def call_promoters_related_bins(self):

        promoter_overlap_file_name = "promoter/new_overlap_promoter_version3.csv"
        promoter_overlap_data = pd.read_csv(promoter_overlap_file_name, header=None, low_memory=False)
        promoter_overlap_data.columns = list(promoter_overlap_data.loc[0])
        # print(list(promoter_overlap_data.loc[0]))
        overlaped_bins = promoter_overlap_data.loc[1:]["overlaped_bin_F_ID"]

        v = np.vectorize(self.promoters_related_bins)
        v(overlaped_bins)
        self.promoter_overlaps_bin = set(self.promoter_overlaps_bin)
        # self.promoter_overlaps_bin = sorted(self.promoter_overlaps_bin)
        # print(len(self.promoter_overlaps_bin))

    def promoters_related_bins(self, bin):
        self.counter += 1
        x = ast.literal_eval(bin)
        for i in x:
            self.promoter_overlaps_bin.append(int(i))

    def call_seprate_related_promoters_bins(self):

        self.bin_dict = self.make_bin_dictionary_interactions()
        self.call_promoters_related_bins()
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

            # find elements which interact by promoter

    # return id of region that promoter interacted with
    def interacted_id_promoter(self, csv_file_name="new_overlap_promoter_version3.csv"):

        promoter_data = pd.read_csv(csv_file_name, header=None)
        promoter_data.columns = list(promoter_data)
        promoter_dict = {}

        for n in range(1, len(promoter_data)):
            interacted_id = promoter_data.at[n, "interacted_id"]
            interacted_ids = ast.literal_eval(interacted_id)

            promoter_id = promoter_data.at[n, "promoter_id"]
            for id in interacted_ids:
                id = int(id)
                if id in promoter_dict:
                    if promoter_id not in promoter_dict[id]:
                        promoter_dict[id].append(promoter_id)
                else:
                    promoter_dict[id] = [promoter_id]

        print(len(set(promoter_dict)))
        return promoter_dict

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
                temp = temp.reset_index(drop=True)
                temp = temp.sort_values(by=['pos'])
                temp = temp.reset_index(drop=True)
                seprated_by_chr.append(temp)

        if mode == 1:
            for chr in range(1, 24):
                # chro = "chr" + str(chr)
                chr_str = str(chr)
                temp = data.query("chr == @chr or chr == @chr_str")
                temp = temp.reset_index(drop=True)
                temp = temp.sort_values(by=['start'])
                temp = temp.reset_index(drop=True)
                seprated_by_chr.append(temp)
            s = 0

            # m = [len(i) for i in seprated_by_chr]
            # print(sum(m))

        return seprated_by_chr

    def GWAS_annotation(self,file_name, data ,element_file_name, permutation = False):
        if not permutation :
            # GWAS_file_name = "GWAS/breast_cancer_GWAS.tsv"
            GWAS_file_name = file_name
            # GWAS_data = pd.read_csv(GWAS_file_name, sep="\t", header=None)
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

        element_id = element_data.loc[:]["id"]
        element_start = element_data.loc[:]["start"]
        element_end = element_data.loc[:]["end"]
        element_chr = element_data.loc[:]["chr"]

        GWAS_chr = GWAS_data.loc[:]["chr"]
        GWAS_pos = GWAS_data.loc[:]["pos"]

        self.GWAS_seprated_by_chr = self.seprate_chr(GWAS_data)
        self.element_data_seprated_by_chr = self.seprate_chr(element_data, mode=1)

        v = np.vectorize(self.check_element_GWAS)
        v(element_id, element_start, element_end, element_chr)
        # print(self.element_GWAS_dict)
        # print("number of overlap : ",len(self.element_GWAS_dict))
        self.all.append(self.element_GWAS_dict)

        # w = np.vectorize(self.check_GWAS_element)
        # w(GWAS_chr, GWAS_pos)
        # print(self.GWAS_element_dict)
        # print(len(self.GWAS_element_dict))

    def check_element_GWAS(self, element_id, element_start, element_end, element_chr):
        if element_chr == "X" or element_chr == "Y":
            element_chr = 23
        try:
            element_chr = int(element_chr)
        except:
            element_chr = 24
        if element_chr < 24:
            mutation_pos = list(self.GWAS_seprated_by_chr[element_chr - 1].loc[:]["pos"])
            start_index = np.searchsorted(mutation_pos,element_start)
            for index in range(start_index, len(mutation_pos)):
                if mutation_pos[index] > element_end :
                    break

                if mutation_pos[index] >= element_start and mutation_pos[index] <= element_end:
                    if element_id not in self.element_GWAS_dict:
                        self.element_GWAS_dict[int(element_id)] = [mutation_pos[index]]
                    else:
                        self.element_GWAS_dict[int(element_id)].append(mutation_pos[index])

    def check_GWAS_element(self, GWAS_chr, GWAS_pos):

        if GWAS_chr == "chrX" or GWAS_chr == "chrY":
            GWAS_chr = "chr23"

        GWAS_chr = int(GWAS_chr.replace("chr", ""))
        element_data = self.element_data_seprated_by_chr[GWAS_chr - 1]

        element_id = list(element_data.loc[:]["id"])
        element_start = list(element_data.loc[:]["start"])
        element_end = list(element_data.loc[:]["end"])

        start_index = np.searchsorted(element_end, GWAS_pos)
        # end_index = np.searchsorted(bin_start, end)

        for index in range(start_index, len(element_data)):
            if GWAS_pos > element_end[index]:
                break
            if GWAS_pos >= element_start[index] and GWAS_pos <= element_end[index]:
                if GWAS_pos not in self.GWAS_element_dict:
                    self.GWAS_element_dict[GWAS_pos] = [element_id[index]]
                else:
                    self.GWAS_element_dict[GWAS_pos].append(element_id[index])

    def make_dict_of_bins_and_element_overlap(self, file_name):
        dict = {}
        element_data = pd.read_csv(
            "enhancer/bin_id_interaction_overlap_" + file_name + "_file" + ".csv", header=None)
        element_data.columns = ['bin_id', 'element_id']

        element_bins = element_data["bin_id"]
        element_ids = element_data["element_id"]

        for i, j in zip(element_bins, element_ids):
            dict[int(i)] = ast.literal_eval(j)

        len(dict)
        return dict

    def map_GWAS_to_bins(self):
        # element_data = pd.read_csv(
        #     "enhancer/bin_id_interaction_overlap_" + file_name + "_file" + ".csv", header=None)
        # element_data.columns = ['bin_id', 'element_id']
        # x = []
        # for i in self.GWAS_element_dict :
        #     for k in  self.GWAS_element_dict[i] :
        #         if k not in x :
        #             x.append(k)
        # print(len(x))
        element_data = pd.read_csv("bin_data/new_unique_bin.csv")
        element_data.columns = list(element_data)
        element_promoter_bins = []
        element_non_promoter_bins = []
        element_both_bins = []
        element_data.astype({'id': 'int64'})
        # element_bins = element_data.loc[:]["bin_id"]
        element_bins = list(set(element_data.loc[:]["id"]))
        n = 0
        # print(type(self.promoters_bins[0]))
        # print(type(element_bins[0]))

        s_promoters_bin = set(self.promoters_bins)
        # print(len(self.promoters_bins), len(s_promoters_bin))

        s_non_promoters_bin = set(self.non_promoters_bins)
        # print(len(self.non_promoters_bins), len(s_non_promoters_bin))

        s_both_bin = set(self.both_bins)
        # print(len(self.both_bins), len(s_both_bin))

        promoter_GWAS = 0
        non_promoter_GWAS = 0
        both_GWAS = 0
        pG = []
        npG = []
        bG = []
        z = []
        # bin_element_dict = self.make_dict_of_bins_and_element_overlap(file_name)
        m = 0
        # print(" ============== ")
        # print(len(set(self.GWAS_element_dict)))
        # # print(len(element_bins))
        # x = 0
        # for gwas in self.GWAS_element_dict :
        #     for j in self.GWAS_element_dict[gwas] :
        #         x+=1
        #         if int(j) in s_promoters_bin :
        #             promoter_GWAS+=1
        #         if int(j) in s_non_promoters_bin :
        #             non_promoter_GWAS+=1
        #         if int(j) in s_both_bin:
        #             both_GWAS+=1



        # for bin in element_bins :
        #     if bin in self.element_GWAS_dict :
        #         m+= len(self.element_GWAS_dict[bin])
        #
        #     if n % 10000 == 0:
        #         print(n)
        #         n+=1
        #
        #     if bin in s_promoters_bin :
        #         element_promoter_bins.append(bin)
        #         if bin in self.element_GWAS_dict :
        #             print(bin)
        #             for j in self.element_GWAS_dict[bin] :
        #                 if j not in z :
        #                     z.append(j)
        #                 if j not in pG :
        #                     pG.append(j)
        #
        #     elif bin in s_non_promoters_bin:
        #         element_non_promoter_bins.append(bin)
        #         if bin in self.element_GWAS_dict:
        #             for j in self.element_GWAS_dict[bin] :
        #                 if j not in z :
        #                     z.append(j)
        #                 if j not in npG :
        #                     npG.append(j)
        #
        #
        #     elif bin in s_both_bin:
        #         element_both_bins.append(bin)
        #         if bin in self.element_GWAS_dict:
        #             for j in self.element_GWAS_dict[bin] :
        #                 if j not in z :
        #                     z.append(j)
        #                 if j not in bG :
        #                     bG.append(j)
        #

        print(len(element_bins))
        for bin in element_bins:
            # if n % 10000 == 0:
            #     print(n)
            #     n += 1

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

        print(len(element_promoter_bins))
                    # if n % 10000 == 0:
                    #     print(n)
                    # n+=1
                    # if bin in s_promoters_bin :
                    #     element_promoter_bins.append(bin)
                    #     elements_ids = [int(i) for i in bin_element_dict[bin]]
                    #     for id in elements_ids :
                    #         if id in self.element_GWAS_dict :
                    #             promoter_GWAS += 1
                    #             # promoter_GWAS += len(self.element_GWAS_dict[id])
                    #
                    # elif bin in s_non_promoters_bin:
                    #     element_non_promoter_bins.append(bin)
                    #     elements_ids = [int(i) for i in bin_element_dict[bin]]
                    #     for id in elements_ids:
                    #         if id in self.element_GWAS_dict:
                    #             non_promoter_GWAS += 1
                    #         #     non_promoter_GWAS += len(self.element_GWAS_dict[id])
                    #
                    # elif bin in s_both_bin:
                    #     element_both_bins.append(bin)
                    #     elements_ids = [int(i) for i in bin_element_dict[bin]]
                    #     for id in elements_ids:
                    #         if id in self.element_GWAS_dict:
                    #             promoter_GWAS += 1
                    # both_GWAS += 1
                    #     both_GWAS += len(self.element_GWAS_dict[id])
        #
        # s_element_promoter_bins = set(element_promoter_bins)
        # s_element_non_promoter_bins = set(element_non_promoter_bins)
        # s_element_both_bins = set(element_both_bins)
        # print("results : ")
        # print(promoter_GWAS)
        # print(non_promoter_GWAS)
        # print(both_GWAS)

        # print(promoter_GWAS + non_promoter_GWAS + both_GWAS)
        # print(m)
        # print(len(pG))
        # print(len(npG))
        # print(len(bG))
        # print(len(pG) + len(npG) + len(bG))
        # print(len(z))

        # print(len(element_promoter_bins) - promoter_GWAS)
        # print(len(element_non_promoter_bins) - non_promoter_GWAS)
        # print(len(element_both_bins) - both_GWAS)
        #
        # print("only promoter interactions : ", promoter_GWAS / len(element_promoter_bins))
        print("only non-promoter interactions : ", non_promoter_GWAS / len(element_non_promoter_bins))

        print("only promoter interactions and both : ",
              (promoter_GWAS + both_GWAS) / (len(element_promoter_bins) + len(element_both_bins)))

        # print((promoter_GWAS + both_GWAS) / (len(element_promoter_bins) + len(element_both_bins)))
        # self.results.append((promoter_GWAS + both_GWAS) / (len(element_promoter_bins) + len(element_both_bins)))
        # return (promoter_GWAS + both_GWAS) / (len(element_promoter_bins) + len(element_both_bins))
        print("all : ", (promoter_GWAS + non_promoter_GWAS + both_GWAS),
              len(element_promoter_bins) + len(element_non_promoter_bins) + len(element_both_bins),
              (promoter_GWAS + non_promoter_GWAS + both_GWAS) / (
              len(element_promoter_bins) + len(element_non_promoter_bins) + len(element_both_bins)))



        # print(len(element_promoter_bins) + len(element_non_promoter_bins) + len(element_both_bins))
        # print(len(element_promoter_bins) + len(element_both_bins))

    def draw_bar_chart(self, m, n, k, name,GWAS_type):

        # objects = ("Bin that \n only have interactions \n with promoters",
        #            "Bin that \n only have interactions \n with non-promoters",
        #            "All bins ")

        objects = ("Bin that \n  have interactions \n with promoters",
                   "Bin that \n  have interactions \n with non-promoters",
                   "All bins ")


        y_pos = np.arange(3)
        performance = [m * 100, n * 100, k * 100]

        bar_list = plt.bar(y_pos, performance, align='center', alpha=0.5, width=0.2)
        bar_list[0].set_color('c')
        bar_list[1].set_color('turquoise')

        plt.xticks(y_pos, objects)
        plt.ylabel('Percentage of bins which overlap with '+GWAS_type+' GWAS')
        # plt.xlim(-0.8, 1 + .5, 2 )
        plt.title(name + ' and '+GWAS_type +' GWAS ')
        plt.savefig("all_new_charts_3" + "/" + name + ".png")
        # plt.show()
        plt.savefig(name + ".png")

# gwas finding
# bins finding

p = Hic_proj_GWAS()
# p.seprate_schizophrni_GWAS()
# p.call_seprate_related_promoters_bins()
# print(" =============================== ")
# p.GWAS_annotation("GWAS/not_Schizophrenia_GWAS.csv","","bin_data/new_unique_bin.csv",False)
# print("================================")
# x = 0
# for i in p.element_GWAS_dict :
#     x += len(p.element_GWAS_dict[i])
# print(x)
p.map_GWAS_to_bins()
# print("finish")

# p.draw_bar_chart(0.12,0.39,0.49,"all bins")
# p.draw_bar_chart(0.00574,0.00397,0.00495,"Chart_v2 : bins ","breast cancer")
# p.draw_bar_chart(0.00574,0.00397,0.00495,"Chart_v2 : bins ","breast cancer")
p.draw_bar_chart(0.0132,0.0128,0.0130,"v3 bins ","not Schizophrenia")
