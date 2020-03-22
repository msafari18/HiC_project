from gtfparse import read_gtf
import csv
import pandas as pd
import numpy as np
import ast


# from datetime import datetime
class Hic_proj():
    def __init__(self, gene_gtf_file_name, interaction_txt_file_name, prmoters_csv_file_name, overlaps_csv_file_name):
        self.probable_promoters = []
        self.gene_gtf_file_name = gene_gtf_file_name
        self.interaction_txt_file_name = interaction_txt_file_name
        self.prmoters_csv_file_name = prmoters_csv_file_name
        self.exons_data = None
        self.interaction_data = None
        self.promoters_data = None
        self.promoters_overlap_data = None
        self.overlaps = []
        self.interaction_data_seprated_chr_F1 = []
        self.interaction_data_seprated_chr_F2 = []
        self.overlaps_csv_file_name = overlaps_csv_file_name
        self.bin_data = None
        self.bin_seprated = []
        self.interaction = []
        self.is_exist = set()
        self.bin_dict = {}

    def find_intron(self):

        data_frame = read_gtf("t.gtf",
                              usecols=['start', 'end', 'gene_id', 'seqname', 'feature'])
        exons_pos = {}
        gene_exons = {}
        gene_pos = {}
        gene_id = data_frame.gene_id.unique()
        exon_data = data_frame.query('feature == "exon"')
        exon_data = exon_data.reset_index(drop=True)
        print(len(exon_data))
        gene_data = data_frame.query('feature == "gene"')
        gene_data = gene_data.reset_index(drop=True)
        print(len(gene_data))
        for i in range(len(gene_data)):
            id = gene_data.at[i, "gene_id"]
            if id not in gene_pos:
                gene_pos[id] = [gene_data.at[i, "start"], gene_data.at[i, "end"]]

        print(len(gene_pos))
        counter = 0
        print("here")
        for n, id in enumerate(gene_id):
            if n % 1000 == 0:
                print(n)
            gene_ex = exon_data.query('gene_id == @id')
            gene_ex = gene_ex.reset_index(drop=True)
            pos = []
            for e in range(len(gene_ex)):
                position = [gene_ex.at[e, "start"], gene_ex.at[e, "end"]]
                if position not in pos:
                    pos.append(position)
                    ex_id = "exon_" + str(counter) + "_" + str(gene_ex.at[e, "seqname"])
                    exons_pos[ex_id] = position
                    if id not in gene_exons:
                        gene_exons[id] = [ex_id]
                    else:
                        gene_exons[id].append(ex_id)
                    counter += 1

        print(len(gene_exons))
        print(len(exons_pos))

        with open('exons/gene_pos.csv', 'w', newline="") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in gene_pos.items():
                writer.writerow([key, value])

        with open('exons/exons_pos.csv', 'w', newline="") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in exons_pos.items():
                writer.writerow([key, value])

        with open('exons/gene_exons.csv', 'w', newline="") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in gene_exons.items():
                writer.writerow([key, value])

        with open('exons/gene_exons.csv') as csv_file:
            reader = csv.reader(csv_file)
            gene_exons = dict(reader)
        with open('exons/exons_pos.csv') as csv_file:
            reader = csv.reader(csv_file)
            exons_pos = dict(reader)
        with open('exons/gene_pos.csv') as csv_file:
            reader = csv.reader(csv_file)
            gene_pos = dict(reader)

        # print(gene_exons)
        # print(exons_pos)
        # print(gene_pos)
        introns_pos = {}
        counter = 0
        nn = 0

        for gene in gene_exons:
            if nn % 1000 == 0:
                print(nn)
            nn += 1
            gene_position = ast.literal_eval(gene_pos[gene])
            gene_ex = ast.literal_eval(gene_exons[gene])
            positions = []
            for exon in gene_ex:
                pos = ast.literal_eval(exons_pos[exon])
                positions.append((int(pos[0]), int(pos[1])))

            positions = sorted(positions, key=lambda x: x[0])

            for i in range(len(positions) - 2):
                if positions[i][1] > positions[i + 1][0]:
                    continue
                p = [positions[i][1], positions[i + 1][0]]
                intron_id = "intron_" + str(counter) + "_" + exon.split("_")[2]
                introns_pos[intron_id] = p
                counter += 1

            if gene_position[0] != positions[0][0]:
                intron_id = "intron_" + str(counter) + "_" + exon.split("_")[2]
                introns_pos[intron_id] = [gene_position[0], positions[0][0]]
                counter += 1
            intron_id = "intron_" + str(counter) + "_" + exon.split("_")[2]
            introns_pos[intron_id] = [positions[-1][1], gene_position[1]]
            counter += 1

        with open('exons/gene_introns.csv', 'w', newline="") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in introns_pos.items():
                writer.writerow([key, value])

    def pre_process_intron_exons(self):

        with open('exons/exons_pos.csv') as csv_file:
            reader = csv.reader(csv_file)
            exons = dict(reader)

        with open('exons/introns_pos.csv') as csv_file:
            reader = csv.reader(csv_file)
            introns = dict(reader)

        new_introns = []
        new_exons = []
        for exon in exons:
            pos = ast.literal_eval(exons[exon])
            start = pos[0]
            end = pos[1]
            chr = exon.split("_")[2]
            new_exons.append({"id": exon, "start": start, "end": end, "chr": chr})

        for intron in introns:
            pos = ast.literal_eval(introns[intron])
            start = pos[0]
            end = pos[1]
            chr = intron.split("_")[2]
            new_introns.append({"id": intron, "start": start, "end": end, "chr": chr})

        csv_file = "elements/introns.csv"
        csv_columns = ['id', 'start', 'end', 'chr']

        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in new_introns:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

        csv_file = "elements/exons.csv"
        csv_columns = ['id', 'start', 'end', 'chr']

        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in new_exons:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

    def read_gene_gtf_file(self):
        # ['start','end','gene_name','gene_id','seqname','exon_number','feature']
        data_frame = read_gtf(self.gene_gtf_file_name,
                              usecols=['start', 'end', 'seqname', 'feature', 'gene_biotype'])

        data = data_frame.query("feature == 'gene'")
        data = data.reset_index(drop=True)
        data.to_csv("gene_file.csv")
        # data_frame.to_csv(s   elf.gene_gtf_file_name)
        # exon_data1 = exon_data.query('exon_number == "1"')
        # x = data_frame.gene_name.unique()
        # print(len(x))
        # x = data_frame.feature.unique()
        # print(x)

        # pgene = data_frame.query('gene_biotype == "pseudogene"')
        # pgene = pgene.reset_index(drop=True)
        # pgene.drop_duplicates(keep="last", inplace=True)
        # pgene = pgene.reset_index(drop=True)
        # pgene.to_csv("pseudogene_file.csv", header=True, index=True)

        # rRNA = data_frame.query('gene_biotype == "rRNA"')
        # rRNA = rRNA.reset_index(drop=True)
        # rRNA.drop_duplicates(keep="last", inplace=True)
        # rRNA = rRNA.reset_index(drop=True)
        # rRNA.to_csv("rRNA_file.csv", header=True, index=True)
        #
        # lincRNA = data_frame.query('gene_biotype == "lincRNA"')
        # lincRNA = lincRNA.reset_index(drop=True)
        # lincRNA.drop_duplicates(keep="last", inplace=True)
        # lincRNA = lincRNA.reset_index(drop=True)
        # lincRNA.to_csv("lincRNA_file.csv", header=True, index=True)
        # #
        # snRNA = data_frame.query('gene_biotype == "snRNA"')
        # snRNA = snRNA.reset_index(drop=True)
        # snRNA.drop_duplicates(keep="last", inplace=True)
        # snRNA = snRNA.reset_index(drop=True)
        # snRNA.to_csv("snRNA_file.csv", header=True, index=True)
        # #
        # miRNA = data_frame.query('gene_biotype == "miRNA"')
        # miRNA = miRNA.reset_index(drop=True)
        # miRNA.drop_duplicates(keep="last", inplace=True)
        # miRNA = miRNA.reset_index(drop=True)
        # miRNA.to_csv("miRNA_file.csv", header=True, index=True)
        #
        # misc_RNA = data_frame.query('gene_biotype == "misc_RNA"')
        # misc_RNA = misc_RNA.reset_index(drop=False)
        # misc_RNA.drop_duplicates(keep="last", inplace=True)
        # misc_RNA.to_csv("miscRNA_file.csv", header=True, index=False)
        #
        # snoRNA = data_frame.query('gene_biotype == "snoRNA"')
        # snoRNA = snoRNA.reset_index(drop=True)
        # snoRNA.drop_duplicates(keep="last", inplace=True)
        # snoRNA = snoRNA.reset_index(drop=True)
        # snoRNA.to_csv("snoRNA_file.csv", header=True, index=True)

        # exons = data_frame.query('feature == "exon"')
        # exons = exons.reset_index(drop=True)
        # exons.drop_duplicates(keep="last", inplace=True)
        # exons = exons.reset_index(drop=True)
        # exons.to_csv("exons_file.csv", header=True, index=True)
        #
        # gene = data_frame.query('feature == "gene"')
        # gene = gene.reset_index(drop=True)
        # gene.to_csv("gene_file.csv", header=True, index=False)

        # print("12")

        # self.write_parts_to_csv("rRNA_file",rRNA)
        # self.write_parts_to_csv("lincRNA_file",lincRNA)
        # self.write_parts_to_csv("protein_coding_file",protein_coding_part)
        # self.write_parts_to_csv("snRNA_file",snRNA)
        # print("1")
        # self.write_parts_to_csv("miRNA_file",miRNA)
        # self.write_parts_to_csv("miscRNA_file",misc_RNA)
        # self.write_parts_to_csv("snoRNA_file", snoRNA)
        # self.write_parts_to_csv("exon_file",exons)
        # self.write_parts_to_csv("gene_file", gene)
        # print("1")

        # self.exons_data = exon_data1.reset_index(drop=True)

    def read_overlaps_csv_file(self):

        self.overlaps_data = pd.read_csv(self.overlaps_csv_file_name, header=None)
        self.overlaps_data.columns = ["promoter_ID", "promoter_gene_name", "overlaped_element_F_ID",
                                      "overlaped_element_start", "overlaped_element_end", "overlaped_element_INT_ID",
                                      "number_of_overlap_bp"]
        m = 0
        print(self.overlaps_data.shape)
        for i in range(len(self.overlaps_data)):
            m += len(self.overlaps_data.loc[i]["overlaped_element_F_ID"])
            print(self.overlaps_data.loc[i]["overlaped_element_F_ID"])

        print("number of overlap : ", m)

    def read_promter_csv_file(self):
        print(self.prmoters_csv_file_name)
        self.promoters_data = pd.read_csv(self.prmoters_csv_file_name, header=None)
        self.promoters_data.columns = ["promoter_id", "start", "end", "gene_name", "gene_id", "chr"]

        # self.promoters_data = self.promoters_data.loc[:100][:]

    def read_interaction_txt_file(self, file_name="", mode=0):

        if mode == 1:
            interaction_data = pd.read_csv(self.interaction_txt_file_name, sep="\t", header=None)
            interaction_data.columns = ['INT_ID', 'F1_ID', 'F1_chr', 'F1_start', 'F1_end', 'F5promoter', 'F5enhancer',
                                        'F5_cell_enhancer', 'F5_tissue_enhancer', 'h3k27ac', 'h3k9me3', 'CTCF',
                                        'hmm_fraction', 'hmm', 'promoter_region_of_gene', 'gene_symbol', 'gene_type',
                                        'F2_ID', 'F2_chr', 'F2_start', 'F2_end', 'F5promoter', 'F5enhancer', 'F5_cell_enhancer',
                                        'F5_tissue_enhancer', 'h3k27ac', 'h3k9me3', 'CTCF', 'hmm_fraction', 'hmm',
                                        'promoter_region_of_gene', 'gene_symbol', 'gene_type', 'distance',
                                        'interacting_read', 'logpvalue']
            print(list(interaction_data.loc[0]))
            interaction_dict = []
            interaction_data = interaction_data.loc[1:]
            new_dataframe = interaction_data[["INT_ID", "F1_ID", "F1_chr", "F1_start", "F1_end", "F2_ID", "F2_chr",
                                              "F2_start",
                                              "F2_end", "distance",
                                              "interacting_read", "logpvalue"]]

            new_dataframe.to_csv("Neun-/Neun-_Neun_HiC.csv",index=False)


        if mode == 0:
            bin_data = pd.read_csv(file_name, header=None)
            bin_data.columns = ["bin_id", "start", "end", "chr"]
            bin_data = bin_data.loc[1:].astype('int64')
            bin_data = bin_data.sort_values(by=['start'])
            self.seprate_by_chr(bin_data)
            self.interaction_data = bin_data

            # self.seprate_by_chr(self.interaction_data)

    def seprate_by_chr(self, data):

        print(len(data))
        for chr in range(1, 25):
            s = str(chr)
            temp = data.query("chr == @chr or chr == @s")
            temp = temp.reset_index(drop=True)
            temp = temp.sort_values(by=['start'])
            self.interaction_data_seprated_chr_F1.append(temp)

            # temp = data.query("F2_chr == @chr or F2_chr == @s")
            # temp = temp.reset_index(drop=True)
            # temp = temp.sort_values(by=['F2_start'])
            # self.interaction_data_seprated_chr_F2.append(temp)

        m = [len(i) for i in self.interaction_data_seprated_chr_F1]
        n = [len(j) for j in self.interaction_data_seprated_chr_F2]

        # print(sum(m), sum(n))

    def find_probable_promoters(self):
        self.read_gene_gtf_file()
        for i in range(len(self.exons_data)):
            start = self.exons_data.at[i, "start"]
            gene_name = self.exons_data.at[i, "gene_name"]
            seq_name = self.exons_data.at[i, "seqname"]
            gene_id = self.exons_data.at[i, "gene_id"]
            self.probable_promoters.append(
                {"promoter_id": "promoter_" + str(i), "start": start - 5000, "end": start + 1000,
                 "gene_name": gene_name, "gene_id": gene_id, "chr": seq_name})

        self.write_promoter_places_to_csv()

    def write_promoter_places_to_csv(self):

        csv_file = "promoter/promoters_with_chr_v3.csv"
        csv_columns = ['promoter_id', 'start', 'end', 'gene_name', 'gene_id', 'chr']
        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in self.probable_promoters:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

    def make_bin_dictionary_interactions(self):
        p.read_promter_csv_file()

        bins_dicts = {}
        unique = []
        temp = []
        interaction_file_name = self.interaction_txt_file_name

        interaction_data = pd.read_csv(interaction_file_name, header=None)
        interaction_data.columns = ['INT_ID', 'F1_ID', 'F1_chr', 'F1_start', 'F1_end', 'F2_ID', 'F2_chr',
                                    'F2_start',
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
            if int(bin) not in bins_dicts:
                bins_dicts[int(bin)] = [int(F2_bin[n])]
            else:
                bins_dicts[int(bin)].append(int(F2_bin[n]))
                bins_dicts[int(bin)] = list(set(bins_dicts[int(bin)]))
        for n, bin in enumerate(F2_bin):
            temp.append(bin)
            if n == 0:
                continue
            if int(bin) not in bins_dicts:
                bins_dicts[int(bin)] = [int(F1_bin[n])]
            else:
                bins_dicts[int(bin)].append(int(F1_bin[n]))
                bins_dicts[int(bin)] = list(set(bins_dicts[int(bin)]))
        self.bin_dict = bins_dicts
        print(len(bins_dicts))
        return bins_dicts

    # def make_bin_dictionary_interactions(self, HiC_file_name = ""):
    #
    #     bins_dicts = {}
    #     unique = []
    #     temp = []
    #     interaction_file_name = self.interaction_txt_file_name
    #     interaction_data = pd.read_csv(interaction_file_name , header=None)
    #     interaction_data.columns = ['INT_ID', 'F1_ID', 'F1_chr', 'F1_start', 'F1_end', 'F2_ID', 'F2_chr', 'F2_start',
    #                                 'F2_end', 'distance', 'interacting_read', 'logpvalue']
    #     interaction_data = interaction_data.astype("str")
    #
    #     F1_bin = list(interaction_data.loc[:]["F1_ID"])
    #     F1_start = list(interaction_data.loc[:]["F1_start"])
    #     F1_end = list(interaction_data.loc[:]["F1_end"])
    #     F1_chr = list(interaction_data.loc[:]["F1_chr"])
    #
    #     F2_bin = list(interaction_data.loc[:]["F2_ID"])
    #     F2_start = list(interaction_data.loc[:]["F2_start"])
    #     F2_end = list(interaction_data.loc[:]["F2_end"])
    #     F2_chr = list(interaction_data.loc[:]["F2_chr"])
    #
    #     for n, bin in enumerate(F1_bin):
    #         temp.append(bin)
    #         if n == 0:
    #             continue
    #         if int(bin) not in bins_dicts:
    #             unique.append(
    #                 {"id": int(F1_bin[n]),
    #                  "start": int(F1_start[n]),
    #                  "end": int(F1_end[n]), "chr": int(F1_chr[n])})
    #             bins_dicts[int(bin)] = [int(F2_bin[n])]
    #         else:
    #             bins_dicts[int(bin)].append(int(F2_bin[n]))
    #     for n, bin in enumerate(F2_bin):
    #         temp.append(bin)
    #         if n == 0:
    #             continue
    #         if int(bin) not in bins_dicts:
    #             unique.append(
    #                 {"id": int(F2_bin[n]),
    #                  "start": int(F2_start[n]),
    #                  "end": int(F2_end[n]), "chr": int(F2_chr[n])})
    #             bins_dicts[int(bin)] = [int(F1_bin[n])]
    #         else:
    #             bins_dicts[int(bin)].append(int(F1_bin[n]))
    #     unique = sorted(unique, key=lambda i: i['id'])
    #     print(len(unique))
    #     # csv_columns = ["id", "start", "end", "chr"]
    #     # try:
    #     #     with open("bin_data/Neun-_new_unique_bin.csv", 'w', newline='') as csvfile:
    #     #         writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    #     #         writer.writeheader()
    #     #         for data in unique:
    #     #             writer.writerow(data)
    #     # except IOError:
    #     #     print("I/O error / can not write on csv file")
    #
    #     return bins_dicts

    def find_overlap_with_promoter_vectorize(self, promoter_id, promoter_start, promoter_end, promoter_chr,
                                             promoter_gene_name, promoter_gene_id):

        if not promoter_start == "start":

            chr = promoter_chr

            p = promoter_id.split("_")
            if int(p[1]) % 10000 == 0:
                print(p[1])

            if chr == "X" or chr == "Y":
                chr = 23
            try:
                chr = int(chr)
            except:
                chr = 24

            if chr >= 1 and chr <= 23:

                temp_interaction = self.interaction_data_seprated_chr_F1[chr - 1]
                promoter_start = int(promoter_start)
                promoter_end = int(promoter_end)
                start = temp_interaction.loc[:]['start']
                end = temp_interaction.loc[:]['end']
                id = temp_interaction.loc[:]['bin_id']
                length = 2000  # or 2000
                start_index = np.searchsorted(end, promoter_start - length)
                end_index = np.searchsorted(start, promoter_end)
                i = start_index

                # if len(temp_interaction) != 0:
                # for i in range(0, len(temp_interaction)):
                while i <= end_index:
                    if end_index >= len(end):
                        break
                    if promoter_start > start[i] - 1990 and promoter_end < end[i] + 1990:
                        if promoter_start < start[i]:
                            bp_overlap = promoter_end - start[i]
                        elif promoter_start >= start[i] and promoter_end <= end[i]:
                            bp_overlap = promoter_end - promoter_start
                        elif promoter_end > end[i]:
                            bp_overlap = end[i] - promoter_start
                        if promoter_id not in self.is_exist:
                            # print(int(id[i]))
                            # print(self.bin_dict[143])
                            # print("=======> ",self.bin_dict[int(id[i])])
                            self.overlaps.append({"promoter_ID": promoter_id,
                                                  "promoter_gene_name": promoter_gene_name,
                                                  "promoter_gene_id": promoter_gene_id,
                                                  "overlaped_element_F_ID": [id[i]],
                                                  "overlaped_element_start": [int(start[i])],
                                                  "overlaped_element_end": [int(end[i])],
                                                  # "overlaped_element_INT_ID": [INT_ID],
                                                  "number_of_overlap_bp": [bp_overlap],
                                                  "related_interactions_id": [k for k in self.bin_dict[int(id[i])]]
                                                  })

                            self.is_exist.add(promoter_id)
                        else:
                            exist_promoter = [sub_dict for sub_dict in self.overlaps if
                                              sub_dict["promoter_ID"] == promoter_id]

                            if id[i] not in exist_promoter[0]["overlaped_element_F_ID"]:
                                exist_promoter[0]["overlaped_element_F_ID"].append(int(id[i]))
                                exist_promoter[0]["overlaped_element_start"].append(start[i])
                                exist_promoter[0]["overlaped_element_end"].append(end[i])
                                exist_promoter[0]["number_of_overlap_bp"].append(bp_overlap)
                                for j in self.bin_dict[id[i]]:
                                    exist_promoter[0]["related_interactions_id"].append(j)
                                exist_promoter[0]["related_interactions_id"] = list(
                                    set(exist_promoter[0]["related_interactions_id"]))
                                # exist_promoter[0]["overlaped_element_INT_ID"].append(INT_ID)

                    i += 1

                    # F2_start = temp_interaction.loc[:]['F2_start']
                    # F2_end = temp_interaction.loc[:]['F2_end']
                    # temp_interaction = self.interaction_data_seprated_chr_F2[chr - 1]
                    # start_index = np.searchsorted(F2_end, promoter_start - length)
                    # end_index = np.searchsorted(F2_start, promoter_end)
                    # i = start_index
                    #
                    # # indexes = temp_interaction.loc[:].index
                    # while i <= end_index:
                    #     # index = indexes[i]
                    #     F2_start = int(temp_interaction.at[i, 'F2_start'])
                    #     F2_end = int(temp_interaction.at[i, 'F2_end'])
                    #     if promoter_end < F2_start :
                    #         break
                    #     elif promoter_start > F2_start - 1990 and promoter_end < F2_end + 1990:
                    #         if F1_interaction[i] == 1:
                    #             F1_interaction[i] = 0
                    #         elif F1_interaction[i] == 0:
                    #             F1_interaction[i] = 2
                    #     i+=1
                    #
                    # F1_interactions = [i for i, n in enumerate(F1_interaction) if n == 1]
                    # F2_interactions = [i for i, n in enumerate(F1_interaction) if n == 2]
                    # temp_interaction = self.interaction_data_seprated_chr_F1[chr - 1]
                    # for i in F1_interactions:
                    #     start = int(temp_interaction.at[i, "F1_start"])
                    #     end = int(temp_interaction.at[i, "F1_end"])
                    #     INT_ID = temp_interaction.at[i, "INT_ID"]
                    #     id = temp_interaction.at[i, "F1_ID"]
                    #
                    #     if promoter_start < start:
                    #         bp_overlap = promoter_end - start
                    #     elif promoter_start >= start and promoter_end <= end:
                    #         bp_overlap = promoter_end - promoter_start
                    #     elif promoter_end > end:
                    #         bp_overlap = end - promoter_start
                    #
                    #     if promoter_id not in self.is_exist :
                    #         self.overlaps.append({"promoter_ID": promoter_id,
                    #                               "promoter_gene_name": promoter_gene_name,
                    #                               "promoter_gene_id": promoter_gene_id,
                    #                               "overlaped_element_F_ID": [id],
                    #                               "overlaped_element_start": [int(start)],
                    #                               "overlaped_element_end": [int(end)],
                    #                               "overlaped_element_INT_ID": [INT_ID],
                    #                               "number_of_overlap_bp": [bp_overlap]
                    #                               })
                    #         self.is_exist.add(promoter_id)
                    #
                    #     else:
                    #         exist_promoter = [sub_dict for sub_dict in self.overlaps if
                    #                           sub_dict["promoter_ID"] == promoter_id]
                    #
                    #         if id not in exist_promoter[0]["overlaped_element_F_ID"]:
                    #             exist_promoter[0]["overlaped_element_F_ID"].append(id)
                    #             exist_promoter[0]["overlaped_element_start"].append(start)
                    #             exist_promoter[0]["overlaped_element_end"].append(end)
                    #             exist_promoter[0]["number_of_overlap_bp"].append(bp_overlap)
                    #         exist_promoter[0]["overlaped_element_INT_ID"].append(INT_ID)
                    #
                    # for i in F2_interactions:
                    #     start = int(temp_interaction.at[i, "F2_start"])
                    #     end = int(temp_interaction.at[i, "F2_end"])
                    #     id = temp_interaction.at[i, "F2_ID"]
                    #     INT_ID = temp_interaction.at[i, "INT_ID"]
                    #
                    #     if promoter_start < start:
                    #         bp_overlap = promoter_end - start
                    #     elif promoter_start >= start and promoter_end <= end:
                    #         bp_overlap = promoter_end - promoter_start
                    #     elif promoter_end > end:
                    #         bp_overlap = end - promoter_start

                    # if promoter_id not in self.is_exist:
                    #     self.overlaps.append({"promoter_ID": promoter_id,
                    #                           "promoter_gene_name": promoter_gene_name,
                    #                           "promoter_gene_id": promoter_gene_id,
                    #                           "overlaped_element_F_ID": [id],
                    #                           "overlaped_element_start": [int(start)],
                    #                           "overlaped_element_end": [int(end)],
                    #                           "overlaped_element_INT_ID": [INT_ID],
                    #                           "number_of_overlap_bp": [bp_overlap]
                    #                           })
                    #     self.is_exist.add(promoter_id)
                    # else:
                    #     exist_promoter = [sub_dict for sub_dict in self.overlaps if
                    #                       sub_dict["promoter_ID"] == promoter_id]
                    #
                    #     if id not in exist_promoter[0]["overlaped_element_F_ID"]:
                    #         exist_promoter[0]["overlaped_element_F_ID"].append(id)
                    #         exist_promoter[0]["overlaped_element_start"].append(start)
                    #         exist_promoter[0]["overlaped_element_end"].append(end)
                    #         exist_promoter[0]["number_of_overlap_bp"].append(bp_overlap)
                    #     exist_promoter[0]["overlaped_element_INT_ID"].append(INT_ID)

    def find_interaction_starts(self):
        # preprocessing ovelap data
        self.bin_data = pd.read_csv("unique_bin.csv", header=None)
        self.bin_data.columns = ["bin_id", "bin_start", "bin_end", "chr"]
        self.bin_data = self.bin_data.loc[1:].astype('int64')
        self.bin_data = self.bin_data.sort_values(by=['bin_start'])
        self.bin_seprated = []

        for chr in range(1, 24):
            temp = self.bin_data.query("chr == @chr")
            temp = temp.reset_index(drop=True)
            temp = temp.sort_values(by=['bin_start'])
            self.bin_seprated.append(temp)

        RNA_file_name = "miRNA_file.csv"
        RNA_data = pd.read_csv(RNA_file_name, header=None)

        RNA_data.columns = ["index", "start", "end", "gene_name", "seqname", "gene_biotype"]
        RNA_data = RNA_data.loc[1:]
        RNA_data = RNA_data.astype({'start': 'int64'})
        RNA_data = RNA_data.astype({'end': 'int64'})
        RNA_data = RNA_data.sort_values(by=['start'])

        v = np.vectorize(self.find_RNA_overlap)
        v(RNA_data.loc[:]["start"], RNA_data.loc[:]["end"], RNA_data.loc[:]["index"], RNA_data.loc[:]["seqname"])

        csv_file_name = "interaction_overlap_miRNA_2.csv"
        csv_col = ['bin_id', 'RNA_id', 'bin_start', 'bin_end', 'RNA_start', 'RNA_end', 'chr']

        try:
            with open(csv_file_name, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_col)
                writer.writeheader()
                for data in self.interaction:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

    def write_promoter_overlaps_to_csv(self):

        csv_file = "promoter/SEP45_overlap_promoter.csv"
        csv_columns = ['promoter_ID', 'promoter_gene_name', 'promoter_gene_id', 'overlaped_element_F_ID',
                       'overlaped_element_start',
                       'overlaped_element_end',
                       'overlaped_element_INT_ID',
                       'number_of_overlap_bp', 'related_interactions_id']
        try:
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in self.overlaps:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")

    def check(self):
        p1 = pd.read_csv("promoter/overlaps_vectorize2.csv", header=None)
        p2 = pd.read_csv("final_promoter_overlap.csv", header=None)
        # print(p1.loc[0][:])
        p1.columns = list(p1.loc[0][:])
        p2.columns = list(p2.loc[0][:])
        id1 = p1.loc[:]["promoter_ID"]
        id2 = p2.loc[:]["promoter_ID"]
        m = set(id1) - set(id2)
        n = set(id2) - set(id1)
        print(m)


gtf_file_name = "Homo_sapiens.GRCh37.87.gtf"
# txt_file_name = "smal_interaction_data.txt"
txt_file_name = "SEP45/SEP45_CNON_HiC.csv"
csv_file_name = "promoter/promoters_with_chr.csv"
overlaps_csv_file_name = "promoter/final_promoter_overlaps_small_v1.csv"

p = Hic_proj(gtf_file_name, txt_file_name, csv_file_name, overlaps_csv_file_name)
# p.read_interaction_txt_file(mode=1)
# p.make_bin_dictionary_interactions(txt_file_name)
# p.find_probable_promoters()
#
# p.read_gene_gtf_file()
# print("finish")
# p.pre_process_intron_exons()
# p.find_probable_promoters()
# p.write_promoter_places_to_csv()

# now = datetime.now()
# current_time = datetime.now()
# print("Current Time =", current_time)



# print(current_time2 - current_time)
###########################################

p.read_promter_csv_file()
print("step1 finished")
p.read_interaction_txt_file("bin_data/SEP45_new_unique_bin.csv")
print("step2 finished")
# p.find_interaction_starts()
s = p.promoters_data.loc[:]["start"]
e = p.promoters_data.loc[:]["end"]
gn = p.promoters_data.loc[:]["gene_name"]
id = p.promoters_data.loc[:]["promoter_id"]
ch = p.promoters_data.loc[:]["chr"]
gn_id = p.promoters_data.loc[:]["gene_id"]
p.make_bin_dictionary_interactions()
vfunc = np.vectorize(p.find_overlap_with_promoter_vectorize)
vfunc(id, s, e, ch, gn, gn_id)

print("step3 finished")
p.write_promoter_overlaps_to_csv()
print("all finished")

# p.check()
# p.make_bin_dictionary_interactions()
# p.read_gene_gtf_file()
# p.read_overlaps_csv_file()
