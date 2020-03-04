import pandas as pd
import csv
import ast
import numpy as np
import matplotlib.pyplot as plt

class non_promoter_analysis():
    def __init__(self):
        self.m = 0
        self.bin_seprated = []
        self.overlaped = []

    def read_data(self):

        data = pd.read_csv("FANTOM_CAT.lv3_robust.info_table.gene.tsv", header=None, sep="\t")
        data.columns = data.loc[0]
        lnc_data = data.query(
            'geneClass == "lncRNA_divergent" or geneClass == "lncRNA_antisense" or geneClass == "lncRNA_intergenic"or geneClass == "lncRNA_sense_intronic"')
        lnc_data = lnc_data.reset_index(drop=True)
        new_lncRNA = []
        counter = 0
        for i in range(len(lnc_data)):
            pos = lnc_data.at[i, "loc"]
            chr = pos.split(":")[0]
            chr = chr.replace("chr", "")
            pos = pos.split(":")[1]
            start = pos.split("-")[0]
            end = pos.split("-")[1]
            end = end.replace(",", "")
            end = end.replace("+", "")
            type = lnc_data.at[i, "geneClass"]
            id = "lncRNA_" + str(counter)
            new_lncRNA.append({"id": id, "start": start, "end": end, "chr": chr, "type": type})
            counter += 1
        csv_file = "elements/lncRNA.csv"
        csv_columns = ['id', 'start', 'end', 'chr', 'type']

        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in new_lncRNA:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

        return new_lncRNA

    def find_lncRNA_promoter(self):

        lnc_data = self.read_data()
        probable_promoters = []
        for i in range(len(lnc_data)):
            start = lnc_data[i]["start"]
            chr = lnc_data[i]["chr"]
            type = lnc_data[i]["type"]
            id = lnc_data[i]["id"] + "_promoter"

            probable_promoters.append(
                {"id": id, "start": int(start) - 1000, "end": int(start) + 1000,
                 "chr": chr, "type": type})

        csv_file = "promoter/lncRNA_promoter.csv"
        csv_columns = ['id', 'start', 'end', 'chr', 'type']

        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in probable_promoters:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

    def non_promoter_regions(self, promoter_file_name, bin_file_name):
        promoter_bins = []
        non_promoter_bins = []
        promoter_interacted_bin = {}

        promoter_overlap_data = pd.read_csv(promoter_file_name, header=None)
        promoter_overlap_data.columns = list(promoter_overlap_data.loc[0][:])

        bin_data = pd.read_csv(bin_file_name, header=None)
        bin_data.columns = list(bin_data.loc[0][:])
        bin_data_id = bin_data.loc[:]["id"]

        for i in range(1, len(promoter_overlap_data)):
            id = promoter_overlap_data.at[i, "promoter_ID"]
            interacted = promoter_overlap_data.at[i, "related_interactions_id"]
            promoter_interacted_bin[id] = ast.literal_eval(interacted)

        for i in promoter_interacted_bin:
            for j in promoter_interacted_bin[i]:
                promoter_bins.append(j)
        promoter_bins = set(promoter_bins)
        all = set(bin_data_id)
        all = set([str(i) for i in all])
        promoter_bins = set([str(i) for i in promoter_bins])

        for i in all:
            if i not in promoter_bins:
                non_promoter_bins.append(i)
        non_promoter_bins = set(non_promoter_bins)
        print("here  : ")
        non_promoter_bins_dict = []
        for i in range(len(bin_data)):
            if i % 1000 == 0:
                print(i)
            id = str(bin_data.at[i, "id"])
            start = str(bin_data.at[i, "start"])
            end = str(bin_data.at[i, "end"])
            chr = str(bin_data.at[i, "chr"])
            if id in non_promoter_bins:
                non_promoter_bins_dict.append({"id": id, "start": start, "end": end, "chr": chr})

        print(len(non_promoter_bins_dict))

        csv_file = "bin_data/non_promoter_bins.csv"
        csv_columns = ['id', 'start', 'end', 'chr']

        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in non_promoter_bins_dict:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

        return non_promoter_bins_dict

    def protein_coding_gene(self,file_name):
        data = pd.read_csv(file_name,header=None)
        data.columns = list(data.loc[0])
        gene_data = data.query("gene_biotype == 'protein_coding' ")
        gene_data = gene_data.reset_index(drop=True)
        gene_data.to_csv("elements/protein_coding.csv",index=False)

    def seprate_by_chr(self, data):

        for chr in range(1, 24):
            s = str(chr)
            temp = data.query("chr == @chr or chr == @s")
            temp = temp.reset_index(drop=True)
            temp = temp.sort_values(by=['start'])
            self.bin_seprated.append(temp)
        m = [len(i) for i in self.bin_seprated]
        print(sum(m))

    def overlap_with_lncRNA(self, lnc_file_name, bin_file_name):

        non_promoter_bins = pd.read_csv(bin_file_name, header=None)
        non_promoter_bins.columns = list(non_promoter_bins.loc[0][:])
        non_promoter_bins = non_promoter_bins[1:]
        non_promoter_bins = non_promoter_bins.astype({'start': 'int64'})
        non_promoter_bins = non_promoter_bins.astype({'end': 'int64'})

        lnc_RNA_promoter = pd.read_csv(lnc_file_name, header=None)
        lnc_RNA_promoter.columns = list(lnc_RNA_promoter.loc[0][:])
        lnc_RNA_promoter = lnc_RNA_promoter[1:]
        lnc_RNA_promoter = lnc_RNA_promoter.astype({'start': 'int64'})
        lnc_RNA_promoter = lnc_RNA_promoter.astype({'end': 'int64'})

        self.seprate_by_chr(non_promoter_bins)

        v = np.vectorize(self.check_lnc_overlap)
        v(lnc_RNA_promoter.loc[:]["start"], lnc_RNA_promoter.loc[:]["end"],lnc_RNA_promoter.loc[:]["id"],
          lnc_RNA_promoter.loc[:]["chr"])

        print(len(set(self.overlaped)))
        print(len(set(self.overlaped)) / 125726)
        print((125726 - len(set(self.overlaped))) / 125726)

    def check_lnc_overlap(self, start, end, id, chr):
        self.m += 1
        if self.m % 10000 == 0:
            print(self.m)
        if chr == "X" or chr == "Y":
            chr = 23
        try:
            chr = int(chr)
        except:
            chr = 24
        if chr < 24:
            bin_data = self.bin_seprated[int(chr) - 1]
            bin_start = bin_data.loc[:]["start"]
            bin_end = bin_data.loc[:]["end"]
            bin_id = bin_data.loc[:]["id"]
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
                    self.overlaped.append(id_BIN)
                i += 1
    def draw_bar_chart(self, m, n, name, type):
        if type == 1:
            objects = ("non-promoter bins that \n overlap with " + name,
                       "non-promoter bins that \n don't have overlap with " + name)
        y_pos = np.arange(len(objects))
        performance = [m * 100, n * 100]
        plt.clf()
        bar_list = plt.bar(y_pos, performance, align='center', alpha=0.5, width=0.2)
        bar_list[0].set_color('c')
        plt.xticks(y_pos, objects)
        plt.ylabel('percentage')
        plt.xlim(-0.6, 1 + .6)
        plt.title(name + ' and non_promoters')
        plt.savefig("non_promoter_chart" + str(type) + "/" + name + ".png")
        # plt.show()

p = non_promoter_analysis()
# p.read_data()
# p.find_lncRNA_promoter()
# p.non_promoter_regions("promoter/final_promoter_overlap.csv", "bin_data/new_unique_bin.csv")
# p.overlap_with_lncRNA("promoter/lncRNA_promoter.csv", "bin_data/non_promoter_bins.csv")
p.draw_bar_chart(0.06,0.93,"protein coding gene",1)
# p.protein_coding_gene("elements/gene_file.csv")