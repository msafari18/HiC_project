from gtfparse import read_gtf
from random import randrange
import numpy as np
import csv
from pandas import read_csv
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from scipy.stats import norm
import random
import time
import ast
import datetime
from multiprocessing import Pool

RNA_NAME = "rRNA"
PROMOTERS_INTERACTIONS_REGIONS = 144859
default_dict = {1: 249230780, 2: 243160372, 3: 197955065, 4: 191013594, 5: 180898964, 6: 171054606, 7: 159025831,
                8: 146279544,
                9: 141150045, 10: 135515571, 11: 134945505, 12: 133815041, 13: 115099310, 14: 107287770, 15: 102518943,
                16: 90288797, 17: 81187695, 18: 78005232, 19: 59110467, 20: 62944081, 21: 48110795, 22: 51237083,
                'X': 155257495, 'Y': 59001391}


class Permutation():
    def __init__(self, gene_gtf_file_name, ):
        self.interacted_len = []
        self.gtf_data = None
        self.RNA_data = None
        self.promoter_data = None
        self.bin_data = None
        self.promoter_seprated = []
        self.chr_length_dict = default_dict
        self.gene_gtf_file_name = gene_gtf_file_name
        self.interaction = []
        self.bin_seprated = []
        self.promoter_dict = {}

    # reading data
    def read_data(self, RNA_file_name, is_RNA=False, is_promoter=False, is_bins=False, is_enhancer=False):
        if not is_RNA and not is_enhancer:
            self.gtf_data = read_gtf(self.gene_gtf_file_name,
                                     usecols=['chr', 'start'])
        if is_RNA:
            path = "elements/"
            self.RNA_data = read_csv(path + RNA_file_name)
            self.RNA_data.columns = ["index", "start", "end", "gene_name", "chr", "gene_biotype"]
            self.RNA_data = self.RNA_data[1:]

        if is_enhancer:
            path = "elements/"
            self.RNA_data = read_csv(path + RNA_file_name)
            self.RNA_data.columns = ["id", "chr", "start", "end", "element"]
            self.RNA_data = self.RNA_data[1:]

        if is_promoter:
            self.promoter_data = read_csv("promoter/promoters_with_chr.csv", header=None)
            self.promoter_data.columns = ["promoter_id", "start", "end", "gene_name", "chr"]
            self.promoter_data = self.promoter_data[1:]
            self.seprate_chr(0)
        # print(self.chr_length_dict)

        if is_bins:
            self.bin_data = pd.read_csv("bin_data/new_unique_bin.csv", header=None)
            self.bin_data.columns = ["bin_id", "bin_start", "bin_end", "chr"]
            self.bin_data = self.bin_data.loc[1:].astype('int64')
            self.bin_data = self.bin_data.sort_values(by=['bin_start'])
            self.seprate_chr(1)

    # seprate data base on chromosome nuber
    def seprate_chr(self, mode=0):
        for chr in range(1, 24):
            temp = self.promoter_data.query("chr == @chr")
            temp = temp.reset_index(drop=True)
            temp = temp.sort_values(by=['start'])
            self.promoter_seprated.append(temp)
            if mode == 1:
                temp = self.bin_data.query("chr == @chr")
                temp = temp.reset_index(drop=True)
                temp = temp.sort_values(by=['bin_start'])
                self.bin_seprated.append(temp)

    # find length of every chr
    def find_chr_length(self):
        for i in range(1, 24):
            if i == 23:

                temp = self.gtf_data.query("chr == 'X'")
                if not temp.empty: self.chr_length_dict['X'] = max(temp.loc[:]["start"])
                temp = self.gtf_data.query("chr == 'Y'")
                if not temp.empty: self.chr_length_dict['Y'] = max(temp.loc[:]["start"])
                continue

            i_str = str(i)
            temp = self.gtf_data.query("chr == @i_str or chr == @i")
            if not temp.empty: self.chr_length_dict[i] = max(temp.loc[:]["start"])

    # make random data
    def make_random_RNA(self):

        new_RNA_list = []
        for RNA_index in range(1, len(self.RNA_data)):
            next = False
            while (not next):
                start = self.RNA_data.at[RNA_index, "start"]
                end = self.RNA_data.at[RNA_index, "end"]
                chr = self.RNA_data.at[RNA_index, "chr"]
                if not (chr == 'X' or chr == 'Y'):
                    try:
                        chr = int(chr)
                    except:
                        chr = chr
                if chr == 23:
                    chr = 'X'
                RNA_length = end - start
                if type(chr) != int:
                    if not (chr == 'X' or chr == 'Y'):
                        next = True
                        continue
                seed = time.time()
                random.seed(seed)
                new_start = randrange(0, self.chr_length_dict[chr] + 1)
                print(new_start)
                new_end = new_start + RNA_length
                if (not self.overlap_with_promoter(new_start, new_end, chr)):
                    next = True
                    new_RNA_list.append({"index": RNA_index, "start": new_start, "end": new_end, "chr": chr})
        # print(new_RNA_list)
        new_RNA_data = pd.DataFrame(new_RNA_list)
        new_RNA_data = new_RNA_data[["index", "start", "end", "chr"]]
        # print(new_RNA_data)
        # print("finish")
        return new_RNA_data

    # new data have overlap with promoter or not ?
    def overlap_with_promoter(self, RNA_start, RNA_end, chr):

        if chr == "X" or chr == "Y":
            chr = 23
        try:
            chr = int(chr)
        except:
            chr = 24

        if chr < 24:
            promoter_data = self.promoter_seprated[chr - 1]
            promoter_start = promoter_data.loc[:]["start"]
            promoter_end = promoter_data.loc[:]["end"]
            length = RNA_end - RNA_start
            start_index = np.searchsorted(promoter_end, RNA_start - length)
            end_index = np.searchsorted(promoter_start, RNA_end)
            i = start_index
            if start_index >= 2:
                i = i - 2
            while i <= end_index:
                if end_index >= len(promoter_end):
                    break
                if RNA_start > promoter_start[i] - length and RNA_end < promoter_end[i] + length:
                    return True
            return False

    # find RNA overlap
    def find_RNAs_overlp(self, data, is_RNA=True, is_exon=False, is_repeat=False):
        if is_RNA:
            self.interaction = []
            RNA_data = data
            # RNA_data.columns = ["index", "start", "end","ge","chr","b"]
            RNA_data.columns = ["index", "start", "end", "chr"]
            RNA_data = RNA_data.loc[1:]
            RNA_data = RNA_data.astype({'start': 'int64'})
            RNA_data = RNA_data.astype({'end': 'int64'})
            RNA_data = RNA_data.sort_values(by=['start'])
            v = np.vectorize(self.check_RNAs_overlap)

            v(RNA_data.loc[:]["start"], RNA_data.loc[:]["end"], RNA_data.loc[:]["index"],
              RNA_data.loc[:]["chr"], True, False, False)
            self.interaction = list(set(self.interaction))
            # print(len(self.interaction))

    # check RNA overlap with interactions
    def check_RNAs_overlap(self, start, end, id, chr, is_RNA, is_exon, is_repeat):

        if is_repeat:
            chr = chr.replace("chr", "")
        if chr == "X" or chr == "Y":
            chr = 23
        try:
            chr = int(chr)
        except:
            chr = 24
        if chr < 24:
            bin_data = self.bin_seprated[int(chr) - 1]
            bin_start = bin_data.loc[:]["bin_start"]
            bin_end = bin_data.loc[:]["bin_end"]
            bin_id = bin_data.loc[:]["bin_id"]
            length = end - start - 10
            start_index = np.searchsorted(bin_end, start - length)
            end_index = np.searchsorted(bin_start, end)
            i = start_index
            if start_index >= 2:
                i = i - 2

            while i <= end_index:
                if end_index >= len(bin_end):
                    break

                if start > bin_start[i] - length and end < bin_end[i] + length:
                    id_BIN = bin_id[i]
                    self.interaction.append(id_BIN)

                i += 1

    # annotate prompter
    def promoter_element_interaction_annotation(self, file_name, is_RNA=True,
                                                is_exon=False, is_repeat=False):
        x = []
        for bin_id in self.interaction:
            if bin_id in self.promoter_dict:
                x.append(bin_id)
        promoter_related = len(set(x))
        m = promoter_related / PROMOTERS_INTERACTIONS_REGIONS

        return m

    # return id of region that promoter interacted with
    def interacted_id_promoter(self, csv_file_name="promoter/final_promoter_overlap.csv"):

        promoter_data = pd.read_csv(csv_file_name, header=None)
        promoter_data.columns = list(promoter_data.loc[0])
        promoter_dict = {}

        for n in range(1, len(promoter_data)):
            interacted_id = promoter_data.at[n, "related_interactions_id"]
            interacted_ids = ast.literal_eval(interacted_id)

            promoter_id = promoter_data.at[n, "promoter_ID"]
            for id in interacted_ids:
                id = int(id)
                if id in promoter_dict:
                    if promoter_id not in promoter_dict[id]:
                        promoter_dict[id].append(promoter_id)
                else:
                    promoter_dict[id] = [promoter_id]

        self.promoter_dict = promoter_dict
        return promoter_dict

    def run_permutation_test(self, data):
        self.find_RNAs_overlp(data)
        res = self.promoter_element_interaction_annotation(RNA_NAME)
        return res

    # fir normal to histogram
    def normal_dist_chart(self, name, len_of_interactions):

        real_val = {"rRNA": 0.1, "miRNA": 0.66, "snRNA": 0.45, "lincRNA": 3, "pseudogene": 3, "enhancer": 44,
                    "Strong_Enhancer": 21, "Insulator": 10, "new_heterochromatin": 30, "new_repressed": 3.7,
                    "new_repetitive": 0.45}
        data = [i * 100 for i in len_of_interactions]
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
        plt.axvline(x=real_val[name], ymax=250, color='r', alpha=0.5)
        plt.savefig(name + "permutation" + ".png")
        onesample_results = scipy.stats.ttest_1samp(data, real_val[name])
        print(name + "p-value : ", onesample_results[1])
        better = 0
        for i in data:
            if i > real_val[name]:
                better += 1

        print("better results : ", better)
        plt.show()


if __name__ == '__main__':

    proj = Permutation("Homo_sapiens.GRCh37.87.gtf")
    # proj.read_data(RNA_NAME + "_file.csv", True, True, True)
    proj.interacted_id_promoter()
    # proj.find_RNAs_overlp("")
    # res = proj.promoter_element_interaction_annotation(RNA_NAME)
    # print(res)
    RNA_NAME = "rRNA"
    # RNA_NAME = "miRNA"
    # RNA_NAME = "snRNA"
    # RNA_NAME = "lincRNA"
    # RNA_NAME = "pseudogene"
    # RNA_NAME = "enhancer"
    # RNA_NAME = "Insulator"
    # RNA_NAME = "new_heterochromatin"
    # RNA_NAME = "new_repressed"
    # RNA_NAME = "new_repetitive"

    # for RNA_NAME in ["rRNA","miRNA","snRNA","lincRNA","pseudogene"] :
    proj.read_data(RNA_NAME + "_file.csv", True, True, True, False)
    iteration = 50
    t1 = datetime.datetime.now()
    datasets = []
    all_res = []
    counter = 0

    for i in range(20):
        print("start")
        counter += 1
        datasets = []
        for j in range(iteration):
            datasets.append(proj.make_random_RNA())

        with Pool(iteration) as p:
            results = p.map(proj.run_permutation_test, datasets)

        for j in results:
            all_res.append(j)
        print("finish round ", counter)

    print(len(all_res))
    # print(all_res)
    t2 = datetime.datetime.now()
    proj.normal_dist_chart(RNA_NAME, all_res)
    print(RNA_NAME, "finish , time : ", t2 - t1)
