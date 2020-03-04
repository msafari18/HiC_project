import pandas as pd
import scipy
from gtfparse import read_gtf
import ast
import numpy as np
import csv
import random
import time
import matplotlib.pyplot as plt
from scipy.stats import norm

default_dict = {1: 249230780, 2: 243160372, 3: 197955065, 4: 191013594, 5: 180898964, 6: 171054606, 7: 159025831,
                8: 146279544,
                9: 141150045, 10: 135515571, 11: 134945505, 12: 133815041, 13: 115099310, 14: 107287770, 15: 102518943,
                16: 90288797, 17: 81187695, 18: 78005232, 19: 59110467, 20: 62944081, 21: 48110795, 22: 51237083,
                'X': 155257495, 'Y': 59001391, 23:155257495}

class EQTL():
    def __init__(self):
        pass

        self.chr_length_dict = default_dict
        self.result = []
        self.new = {}

    def read_EQTL_data(self, file_name):

        EQTL_data = pd.read_csv(file_name, sep='\t', header=None)
        EQTL_data.columns = ['variant_id', 'gene_id', 'chr_one', 'chr_two', 'start', 'end', 'variant_pos',
                             'pval_nominal']

        m = EQTL_data.loc[:]["variant_id"].unique()
        # print("here : unique variant id : ", len(m))
        new = {}
        n = 0
        for id in m :
            if id == "variant_id":
                continue
            # if n % 100000 == 0:
            #     print(n)
            n+=1
            ids = id.split("_")
            try:
                chr = int(ids[0])
            except:
                chr = 'X'

            new_pos = random.randrange(0, self.chr_length_dict[chr] + 1)
            new[id] = new_pos
            # new[id] = ids[1]

        # print(len(new))
        self.new = new
        return EQTL_data

    def make_bin_dict_start_end_dict(self, file_name):

        bin_dict_start_end = {}
        bin_data = pd.read_csv(file_name, header=None)
        bin_data.columns = list(bin_data.loc[0][:])
        for i in range(1, len(bin_data)):
            bin_dict_start_end[str(bin_data.at[i, "id"])] = (int(bin_data.at[i, "start"]), int(bin_data.at[i, "end"]),int(bin_data.at[i, "chr"]))

        return bin_dict_start_end

    def make_promoter_id_variant_dict(self, promoter_file_name, EQTL_file_name):

        promoter_id_variant_dict = {}
        promoter_data = pd.read_csv(promoter_file_name, header=None)
        promoter_data.columns = list(promoter_data.loc[0])
        promoter_data = promoter_data
        EQTL_data = self.read_EQTL_data(EQTL_file_name)

        EQTL_dict = {}
        for i in range(0, len(EQTL_data)):
            gene_id = EQTL_data.at[i, "gene_id"]
            if gene_id not in EQTL_dict:
                EQTL_dict[gene_id] = [EQTL_data.at[i, "variant_id"]]
            else:
                EQTL_dict[gene_id].append(EQTL_data.at[i, "variant_id"])

        # m = [len(EQTL_dict[i]) for i in EQTL_dict]
        # print("here :  number of variant ", sum(m))
        # print("========================= >")
        # print(len(EQTL_data.loc[:]["gene_id"].unique()))
        # print(len(EQTL_dict))
        # m = [len(EQTL_dict[i]) for i in EQTL_dict]
        # print(sum(m))/

        for i in range(1, len(promoter_data)):
            # if i % 100000 == 0:
            #     print(i)
            gene_id = promoter_data.at[i, "gene_id"]
            promoter_id = promoter_data.at[i, "promoter_id"]
            if gene_id in EQTL_dict:
                variant_id = EQTL_dict[gene_id]
                promoter_id_variant_dict[promoter_id] = variant_id

        # print(len(promoter_id_variant_dict))

        return promoter_id_variant_dict

    def make_promoter_interacted_bin(self, file_name):

        promoter_interacted_bin = {}
        promoter_overlap_data = pd.read_csv(file_name, header=None)
        promoter_overlap_data.columns = list(promoter_overlap_data.loc[0][:])

        for i in range(1, len(promoter_overlap_data)):
            id = promoter_overlap_data.at[i, "promoter_ID"]
            interacted = promoter_overlap_data.at[i, "related_interactions_id"]
            promoter_interacted_bin[id] = ast.literal_eval(interacted)

        return promoter_interacted_bin

    def interaction_EQTL_annotation(self, promoter_file_name, EQTL_file_name, bin_data_file_name,
                                    promoter_overlap_file_name):

        bin_data_dict = self.make_bin_dict_start_end_dict(bin_data_file_name)
        print("first_finish")
        promoter_interacted_bins_dict = self.make_promoter_interacted_bin(promoter_overlap_file_name)
        print("second_finish")
        promoter_id_variant_dict = self.make_promoter_id_variant_dict(promoter_file_name, EQTL_file_name)
        print("third_finish")
        for counter in range(100) :
            all_visited_bins = []
            EQTL_in_bin_promoter_interacting_region = set()
            bins = []
            bins_variant_dict = {}
            n = 0
            for id in promoter_id_variant_dict:
                n += 1
                # variants = [int(i.split("_")[1]) for i in promoter_id_variant_dict[id]]
                variants = [int(self.new[i]) for i in promoter_id_variant_dict[id]]
                chr_var = [i.split("_")[0] for i in promoter_id_variant_dict[id] ]
                new_chr_var = []
                for c in chr_var :
                    try :
                        c = int(c)
                        new_chr_var.append(c)
                    except:
                        new_chr_var.append(23)
                variants_id = promoter_id_variant_dict[id]
                variants = set(variants)
                if id in promoter_interacted_bins_dict:
                    interacted_bins = promoter_interacted_bins_dict[id]
                    interacted_bins = set(interacted_bins)

                    for i in interacted_bins :
                        all_visited_bins.append(i)

                    for v_id, v, chr in zip(variants_id, variants, new_chr_var):
                        for bin in interacted_bins:
                            start, end, chr_bin = bin_data_dict[str(bin)]
                            start = int(start)
                            end = int(end)
                            chr_bin = int(chr_bin)
                            # print(type(chr),type(chr_bin))
                            if v >= start and v <= end and chr_bin == int(chr):
                                # if v not in EQTL_in_bin_promoter_interacting_region :
                                EQTL_in_bin_promoter_interacting_region.add(v_id)
                                bins.append(bin)
                                if bin not in bins_variant_dict :
                                    bins_variant_dict[bin] = [v_id]

                                else :
                                    if v_id not in bins_variant_dict[bin] :
                                        bins_variant_dict[bin].append(v_id)

            # print(len(EQTL_in_bin_promoter_interacting_region))
            # print("here :  all visited bins :", len(set(all_visited_bins)))
            # print("here : bins that overlap with variants : ", len(set(bins)))
            EQTL_data = self.read_EQTL_data(EQTL_file_name)
            self.result.append(len(set(bins)))
            print("finish : ",counter,len(set(bins)))
        print(self.result)


    def normal_dist_chart(self, name, data):

        real_val = 11000
        # data = [i * 100 for i in len_of_interactions]
        # print("all results : ", data)
        length_of_data = len(data)
        # Fit a normal distribution to the data:
        mean, std = norm.fit(data)
        # Plot the histogram.
        plt.hist(data, bins=length_of_data, density=True, alpha=1, color='c')
        # Plot the PDF.
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, length_of_data)
        p = norm.pdf(x, mean, std)
        plt.plot(x, p, 'k', linewidth=2)
        title = name + "\n Fit results: mu = %.4f,  std = %.4f" % (mean, std)
        plt.title(title)
        # plt.axvline(x=real_val, ymax=250, color='r', alpha=0.5)
        plt.savefig(name + ".png")
        plt.axvline(x=real_val, ymax=250, color='r', alpha=0.5)
        plt.savefig(name + "2.png")
        onesample_results = scipy.stats.ttest_1samp(data, real_val)
        print(name + "p-value : ", onesample_results[1])
        better = 0
        for i in data:
            if i > real_val:
                better += 1


p = EQTL()
d = p.make_promoter_interacted_bin("promoter/final_promoter_overlap.csv")
bins = []
print("finish")
for i in d :
    for j in d[i] :
        # if j not in bins :
        bins.append(j)
# print(d)
bins = set(bins)
print(len(bins))
# p.interaction_EQTL_annotation("promoter/promoters_with_chr2.csv", "Breast_Mammary_Tissue.signifpairs_GTEX.txt",
#                               "bin_data/new_unique_bin.csv", "promoter/final_promoter_overlap.csv")
# p.normal_dist_chart("EQTL",p.result)