import pandas as pd
import random
import time
import csv
from codes.GWAS_analysis import Hic_proj_GWAS
from scipy.stats import norm
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as npRandom
import seaborn as sns
import statistics
import datetime
from multiprocessing import Pool

default_dict = {1: 249230780, 2: 243160372, 3: 197955065, 4: 191013594, 5: 180898964, 6: 171054606, 7: 159025831,
                8: 146279544,
                9: 141150045, 10: 135515571, 11: 134945505, 12: 133815041, 13: 115099310, 14: 107287770, 15: 102518943,
                16: 90288797, 17: 81187695, 18: 78005232, 19: 59110467, 20: 62944081, 21: 48110795, 22: 51237083,
                'X': 155257495, 'Y': 59001391}


class Permutation_for_GWAS():
    def __init__(self):
        self.p2 = None
        pass

    def seprate_chr(self, data):
        seprated_by_chr = []
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
        length = [len(item) for item in seprated_by_chr]
        return length

    def make_random_GWAS(self, file_name, is_breast=True):
        path = "GWAS/"
        if is_breast:
            GWAS_data = pd.read_csv(path + file_name + ".tsv", sep="\t", header=None,low_memory=False)
        else:
            GWAS_data = pd.read_csv(path + file_name + ".csv", header=None,low_memory=False)

        GWAS_data.columns = ["chr", "pos"]
        GWAS_data = GWAS_data.loc[1:]
        GWAS_data = GWAS_data.astype({'pos': 'int64'})
        length_chr = self.seprate_chr(GWAS_data)
        new_GWAS = []
        s = 0
        for chr in default_dict:
            if chr == 'Y': break
            if chr == 'X' or chr == 'Y':
                chr_index = 23
            else:
                chr_index = chr

            high = default_dict[chr]
            s += length_chr[chr_index - 1]
            new_positions = npRandom.randint(0, high, size=(length_chr[chr_index - 1],))
            for pos in new_positions:
                if chr == 23: chr = 'X'
                new_GWAS.append({"chr": "chr" + str(chr), "pos": pos})

        new_GWAS_data = pd.DataFrame(new_GWAS)
        return new_GWAS_data

    def normal_dist_chart(self, name, len_of_interactions):
        real_val = 0.76
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
        plt.axvline(x=real_val, ymax=250, color='r', alpha=0.5)
        plt.savefig(name + ".png")
        onesample_results = scipy.stats.ttest_1samp(data, real_val)
        print(name + "p-value : ", onesample_results[1])
        better = 0
        for i in data:
            if i > real_val:
                better += 1

        print("better results : ", better)
        plt.show()

    def run_permutation_test(self,data) :
        # p2 = Hic_proj_GWAS()
        # p2.call_seprate_related_promoters_bins()
        print("start")
        self.p2.element_GWAS_dict = {}
        self.p2.GWAS_annotation("", data,"bin_data/new_unique_bin.csv", True)
        res = self.p2.map_GWAS_to_bins()
        print("finish")
        return res


# file_name = "Schizophrenia_GWAS"


if __name__ == '__main__':
    file_name = "Schizophrenia_GWAS"
    p1 = Permutation_for_GWAS()
    p1.p2 = Hic_proj_GWAS()
    iteration = 2
    t1 = datetime.datetime.now()
    p1.p2.call_seprate_related_promoters_bins()
    datasets = []
    all_res = []
    counter = 0
    for i in range(1) :
        counter +=1
        m = [False for i in range(iteration)]
        n = [file_name for i in range(iteration)]

        f = np.vectorize(p1.make_random_GWAS)
        datasets = f(n,m)

        counter = [i for i in range(iteration)]
        with Pool(iteration) as p:
            results = p.map(p1.run_permutation_test ,datasets)

        for j in results :
            all_res.append(j)
        # print("finish : ", counter)

    print(len(all_res))
    print(all_res)
    t2 = datetime.datetime.now()
    p1.normal_dist_chart(file_name, all_res)
    print("finish , time : ", t2 - t1)
