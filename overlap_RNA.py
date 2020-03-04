import pandas as pd
import csv
import datetime
import matplotlib.pyplot as plt;
import numpy as np
import matplotlib.pyplot as plt
import ast
from scipy.stats import chi2_contingency,stats

UNIQUE_BINS_NUMBER = 270352
# PROMOTERS_INTERACTIONS_REGIONS = 144639
PROMOTERS_INTERACTIONS_REGIONS = 150482

NON_PROMOTER_INTERACTIONS_REGIONS = UNIQUE_BINS_NUMBER - PROMOTERS_INTERACTIONS_REGIONS
ALL_INTERACTIONS = 563870
ALL_PROMOTER_INTERACTIONS = 286560
NON_PROMOTER_INTERACTIONS = ALL_INTERACTIONS - ALL_PROMOTER_INTERACTIONS  # 277,310


class Hic_proj_overlap():
    def __init__(self):
        self.data_hetrochromatin_overlap = []
        self.enhancers_data_seprated = []
        self.bin_seprated = []
        self.interaction = []
        self.interaction_vise_versa = []
        self.bin_id_interacted_with_elements = {}
        self.elements_interacted_with_bin_id = {}
        self.m = 0
        self.duplicate = []
        self.temp_bin = []
        self.repeat_bin = set()

    # it use to splite huge data in seprated data
    # record_per_file by default is 500000
    def split_big_data(self, csv_file_name="interaction_overlap_repeat_file.csv", record_per_file=500000):

        data = open(csv_file_name, 'r').readlines()
        header = data[0]
        data.pop(0)
        file = 1
        record_per_file = 500000

        for j in range(len(data)):
            if j % record_per_file == 0:
                write_file = data[j:j + record_per_file]
                write_file.insert(0, header)
                open(str(csv_file_name) + str(file) + '.csv', 'w+').writelines(write_file)
                file += 1

    # read repeat element
    def read_repeat_element(self, csv_file_name="repeat_element_hg19.out",
                            repeat_file_name="repeat/repeatMaskers.csv"):

        csv_col = ["id", "chr", "start_repeat", "end_repeat", "type_of_repeat"]
        i = 0
        with open(repeat_file_name, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_col)
            writer.writeheader()

            with open(csv_file_name, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    line = line.strip().split()
                    dict = {}
                    dict["id"] = "R_" + str(i)
                    dict["chr"] = str(line[4])
                    dict["start_repeat"] = str(line[5])
                    dict["end_repeat"] = str(line[6])
                    dict["type_of_repeat"] = str(line[10])
                    writer.writerow(dict)
                    i += 1

    # seprate repeat element by type
    def seprate_repeat_data(self, repeat_file_name="repeat/repeatMaskers.csv"):
        # csv_col = ["id", "p", "type"]
        # csv_col = ['bin_id', 'repeat_id', 'bin_start', 'bin_end', 'repeat_start', 'repeat_end',
        #                                 'chr',
        #                                 'type_of_repeat', 'overlap_bp']
        repeat_data = pd.read_csv(repeat_file_name, header=None)
        repeat_data.columns = ["id", "chr", "start_repeat", "end_repeat", "type_of_repeat"]

        # print(i,end=", ")
        # repeat_data.columns = csv_col
        # print(len(repeat_data))
        # result = repeat_data.query("type_of_repeat == 'LINE/L1' ")
        result = repeat_data[repeat_data['type_of_repeat'].str.contains("LTR")]
        # print(len(result))
        # x = result.loc[:]["id"]
        # x = set(list(x))
        # print(len(x))
        # # result = result.reset_index(drop=True)
        # # result.to_csv("repeat/LINE_"+repeat_file_name,index=False)
        # result = repeat_data[repeat_data['type'].str.contains("SINE")]
        # print(len(result))
        # # result = result.reset_index(drop=True)
        # # result.to_csv("repeat/SINE_"+repeat_file_name,index=False)
        # result = repeat_data[repeat_data['type'].str.contains("Simple_repeat")]
        # print(len(result))
        result = result.reset_index(drop=True)
        result.to_csv("repeat/LTR_repeatMaskers.csv", index=False)

    def seprate_by_chr(self, data):

        for chr in range(1, 24):
            s = str(chr)
            temp = data.query("chr == @chr or chr == @s")
            temp = temp.reset_index(drop=True)
            temp = temp.sort_values(by=['bin_start'])
            self.bin_seprated.append(temp)

    # find all RNAs overlaps
    def find_RNAs_overlp(self, RNA_file_name, csv_file_name="bin_data/new_unique_bin.csv", is_RNA=True, is_exon=False,
                         is_repeat=False, is_enhancer=False):
        # preprocessing ovelap data
        bin_data = pd.read_csv(csv_file_name, header=None)
        bin_data.columns = ["bin_id","bin_start","bin_end","chr"]
        bin_data = bin_data.loc[1:].astype('int64')
        bin_data = bin_data.sort_values(by=['bin_start'])
        self.seprate_by_chr(bin_data)

        if is_enhancer :
            path = "enhancer/"
        if is_RNA or is_exon :
            path = "elements/"
        data = pd.read_csv(path + RNA_file_name + "_file.csv", header=None)
        data.columns = list(data.loc[0])

        data = data.loc[1:]
        data = data.astype({'start': 'int64'})
        data = data.astype({'end': 'int64'})
        data = data.sort_values(by=['start'])
        print("name of resded data : ", RNA_file_name)
        print("lenth of reading data : ", len(data))
        print("finish reading ! ")
        print("start processing !")
        self.interaction = []
        self.duplicate = []
        self.bin_id_interacted_with_elements = {}
        self.temp_bin = []
        t1 = datetime.datetime.now()
        v = np.vectorize(self.check_RNAs_overlap)
        v(data.loc[:]["start"], data.loc[:]["end"], data.loc[:]["id"],
          data.loc[:]["chr"], "", is_RNA, is_exon, is_repeat)

        #---------------------------------------------------------------------------------
        interaction_dict = self.make_bin_dictionary_interactions()
        print("here to see what happen to remove wrong interactions : ")
        # for i in interaction_dict :
        #     print(type(i))
        # print(len(interaction_dict))
        # print(len(set(interaction_dict)))
        # print(len(self.bin_id_interacted_with_elements))
        # removed_bin = []
        # for j in self.bin_id_interacted_with_elements :
        #     for i in interaction_dict[j] :
        #         if int(i) in self.bin_id_interacted_with_elements :
        #             removed_bin.append(j)
        #             removed_bin.append(int(i))
        # removed_bin = list(set(removed_bin))
        # print("should removed",len(removed_bin))
        # print(removed_bin)
        # print(len(self.bin_id_interacted_with_elements))
        # for i in removed_bin :
        #     self.bin_id_interacted_with_elements.pop(i)
        # print(len(self.bin_id_interacted_with_elements))
        # x = set(self.bin_id_interacted_with_elements)- set(removed_bin)
        # print(len(x))
        # self.bin_id_interacted_with_elements = list(x)


        print("finish finding overlaps ! ")
        csv_file_name = "interaction_overlap/interaction_overlap_" + RNA_file_name + "_file.csv"
        csv_col = ['element_id', 'bin_id', 'bin_start', 'bin_end', 'element_start', 'element_end', 'chr',
                   'overlap_bp']
        csv_file_name2 = "interaction_overlap/bin_id_interaction_overlap_" + RNA_file_name + "_file.csv"

        print("number of bin_id : ", len(self.bin_id_interacted_with_elements))

        if is_repeat:
            # for i in range(1, 12):
            self.interaction = []
            path = "repeat/"
            repeat_data = pd.read_csv(path + RNA_file_name + ".csv", header=None)
            # repeat_file_name = "repeatMasker.csv" + str(i) + ".csv"
            # repeat_data = pd.read_csv(repeat_file_name, header=None, error_bad_lines=False, engine="c")
            repeat_data.columns = ["id", "chr", "start_repeat", "end_repeat", "type_of_repeat"]
            repeat_data = repeat_data.loc[1:]
            repeat_data = repeat_data.astype({'start_repeat': 'int64'})
            repeat_data = repeat_data.astype({'end_repeat': 'int64'})
            repeat_data = repeat_data.sort_values(by=['start_repeat'])
            print("name of resded data : ", RNA_file_name)
            print("lenth of reading data : ", len(repeat_data))
            print("finish reading ! ")
            print("start processing !")
            t1 = datetime.datetime.now()
            v = np.vectorize(self.check_RNAs_overlap)

            v(repeat_data.loc[:]["start_repeat"], repeat_data.loc[:]["end_repeat"], repeat_data.loc[:]["id"],
              repeat_data.loc[:]["chr"], repeat_data.loc[:]["type_of_repeat"], is_RNA, is_exon, is_repeat)

            csv_file_name = "repeat/interaction_overlap_" + RNA_file_name + ".csv"
            csv_col = ['bin_id', 'element_id', 'bin_start', 'bin_end', 'element_start', 'element_end', 'chr',
                       'type_of_repeat', 'overlap_bp']
            csv_file_name2 = "repeat/bin_id_interaction_overlap_" + RNA_file_name + "_file.csv"

        try:
            with open(csv_file_name, 'w', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=csv_col)
                writer.writeheader()
                for data in self.interaction:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")
        try:
            with open(csv_file_name2, 'w', newline='') as f:
                writer = csv.writer(f)
                for row in self.bin_id_interacted_with_elements.items():
                    writer.writerow(row)
        except IOError:
            print("I/O error / can not write on csv file")
        t2 = datetime.datetime.now()
        print("time: ", t2 - t1)

    # check RNA overlap with interactions
    def check_RNAs_overlap(self, start, end, id, chr, type, is_RNA, is_exon, is_repeat):
        self.m += 1
        if self.m % 1000 == 0:
            print(self.m)
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
                    if start < bin_start[i] and end < bin_end[i]:
                        bp_overlap = end - bin_start[i]
                    elif start >= bin_start[i] and end <= bin_end[i]:
                        bp_overlap = end - start
                    elif start >= bin_start[i] and end > bin_end[i]:
                        bp_overlap = bin_end[i] - start
                    elif start < bin_start[i] and end > bin_end[i]:
                        bp_overlap = 5000  # bin_length
                    else:
                        bp_overlap = 0

                    if is_repeat:
                        if id_BIN not in self.repeat_bin :
                            self.repeat_bin.add(id_BIN)
                        # exist_id = {}
                        # for sub_dict in self.interaction:
                        #     if sub_dict["element_id"] == id:
                        #         exist_id = sub_dict
                        #         break
                        # if len(exist_id) == 0:
                        # self.interaction.append(
                        #     {"element_id": id, "bin_id": [id_BIN], "bin_start": [bin_start[i]],
                        #      "bin_end": [bin_end[i]],
                        #      "element_start": start, "element_end": end, "chr": chr, "overlap_bp": [bp_overlap]})
                        # else:
                        #     exist_id = exist_id
                        #     if id_BIN not in exist_id["bin_id"]:
                        #         exist_id["bin_id"].append(id_BIN)
                        #         exist_id["bin_start"].append(bin_start[i])
                        #         exist_id["bin_end"].append(bin_end[i])
                        #         exist_id["overlap_bp"].append(bp_overlap)

                        # if id_BIN not in self.bin_id_interacted_with_elements:
                        #     self.bin_id_interacted_with_elements[id_BIN] = [id]
                        # else:
                        #     if id not in self.bin_id_interacted_with_elements[id_BIN]:
                        #         self.bin_id_interacted_with_elements[id_BIN].append(id)

                                # self.interaction.append(
                                #     {"bin_id": id_BIN, "repeat_id": id, "bin_start": bin_start[i], "bin_end": bin_end[i],
                                #      "repeat_start": start, "repeat_end": end, "chr": chr, "type_of_repeat": type,
                                #      "overlap_bp": bp_overlap})
                    if is_RNA or is_exon:

                        # exist_id = {}
                        # for sub_dict in self.interaction:
                        #     if sub_dict["element_id"] == id:
                        #         exist_id = sub_dict
                        #         break
                        # if len(exist_id) == 0:
                        #     self.interaction.append(
                        #         {"element_id": id, "bin_id": [id_BIN], "bin_start": [bin_start[i]],
                        #          "bin_end": [bin_end[i]],
                        #          "element_start": start, "element_end": end, "chr": chr, "overlap_bp": [bp_overlap]})
                        # else:
                        #     exist_id = exist_id
                        #     if id_BIN not in exist_id["bin_id"]:
                        #         exist_id["bin_id"].append(id_BIN)
                        #         exist_id["bin_start"].append(bin_start[i])
                        #         exist_id["bin_end"].append(bin_end[i])
                        #         exist_id["overlap_bp"].append(bp_overlap)
                        t = set(self.temp_bin)
                        # print("here : ",t)
                        if id_BIN not in t :
                            self.temp_bin.append(id_BIN)
                            self.bin_id_interacted_with_elements[id_BIN] = [id]
                        else:
                            if id not in self.bin_id_interacted_with_elements[id_BIN]:
                                self.bin_id_interacted_with_elements[id_BIN].append(id)
                i += 1

    # find elements which interact by promoter
    def promoter_element_interaction_annotation(self, file_name,p_file_name, is_RNA=True,
                                                is_exon=False, is_repeat=False):
        total_len = []
        x = []
        promoter_interaction_dict = self.interacted_id_promoter(p_file_name)
        if is_repeat:
            element_data = list(self.repeat_bin)
            for bin_id in element_data:
                if bin_id in promoter_interaction_dict:
                    x.append(bin_id)

            promoter_related = len(set(x))
            # print(promoter_related)
            non_promoter_related = len(element_data) - len(set(x))
            m = promoter_related / PROMOTERS_INTERACTIONS_REGIONS
            n = non_promoter_related / NON_PROMOTER_INTERACTIONS_REGIONS
            self.draw_bar_chart(m, n, file_name, 2)
            print(m, n, "type2")
            # element_data = pd.read_csv("repeat/bin_id_interaction_overlap_" + file_name + "_file.csv",
            #                            header=None)
            # element_data.columns = ['bin_id', 'element_id']
            # element_data_2 = pd.read_csv(
            #     "repeat/interaction_overlap_" + file_name + ".csv", header=None)
        # if is_RNA:
        #     path = ""
        element_data = pd.read_csv(
            "interaction_overlap/bin_id_interaction_overlap_" + file_name + "_file" + ".csv", header=None)
        # element_data = pd.read_csv(
        #     "enhancer/bin_id_interaction_overlap_" + file_name + "_file" + ".csv", header=None)
        element_data.columns = ['bin_id', 'element_id']
        # element_data_2 = pd.read_csv(
            # "interaction_overlap/interaction_overlap_" + file_name + "_file" + ".csv", header=None)
        # element_data_2 = pd.read_csv(
        #     "enhancer/interaction_overlap_" + file_name + "_file" + ".csv", header=None)

        element_promoter_interactions = {}
        # print(len(element_data_2))
        types = []
        x = []
        for n in range(0, len(element_data)):
            bin_id = int(element_data.at[n, "bin_id"])
            ids = element_data.at[n, "element_id"]
            ids = ast.literal_eval(ids)

            if bin_id in promoter_interaction_dict:
                x.append(bin_id)
                for id in ids:
                    if id in element_promoter_interactions:
                        for j in promoter_interaction_dict[bin_id]:
                            element_promoter_interactions[id].append(j)
                    else:

                        m = [j for j in promoter_interaction_dict[bin_id]]
                        element_promoter_interactions[id] = m

        if is_repeat:
            final = []
            for repeat in element_promoter_interactions:
                final.append(
                    {"id": repeat, "interacted_promoters": element_promoter_interactions[repeat]})
            csv_file_name = "promoters_overlap/interaction_" + file_name + "_promoter" + ".csv"
            csv_col = ["id", "interacted_promoters", "type"]

        if is_RNA:
            final = []
            for RNA in element_promoter_interactions:
                final.append(
                    {"id": RNA, "interacted_promoters": element_promoter_interactions[RNA]})
            csv_file_name = "promoters_overlap/interaction_" + file_name + "_promoter" + ".csv"
            csv_col = ["id", "interacted_promoters"]

        try:
            with open(csv_file_name, 'w', newline="") as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=csv_col)
                writer.writeheader()
                for data in final:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")


        # print("===============================")
        # print(len(final))
        # print(len(element_data_2))
        # print(len(element_data_2) - len(final))
        # m = len(final) / len(element_data_2)
        # n = (len(element_data_2) - len(final)) / len(element_data_2)
        # print(m, n, "type 1")
        # self.draw_bar_chart(m, n, file_name, 1)

        # print(total_len)
        # print(len(element_data))
        # print(len(set(x)))
        # print(len(element_data) - len(set(x)))

        promoter_related = len(set(x))

        # print(promoter_related)
        non_promoter_related = len(element_data) - len(set(x))

        # print(promoter_related)
        # print(PROMOTERS_INTERACTIONS_REGIONS - promoter_related)
        # print(non_promoter_related)
        # print(NON_PROMOTER_INTERACTIONS_REGIONS - non_promoter_related)

        m = promoter_related / PROMOTERS_INTERACTIONS_REGIONS
        n = non_promoter_related / NON_PROMOTER_INTERACTIONS_REGIONS
        self.draw_bar_chart(m, n, file_name, 2)
        print(m, n, "type2")

        obs = np.array([[promoter_related, non_promoter_related], [PROMOTERS_INTERACTIONS_REGIONS - promoter_related, NON_PROMOTER_INTERACTIONS_REGIONS - non_promoter_related]])
        print(chi2_contingency(obs))

        print("finish : " + file_name)
        return m

    # return id of region that promoter interacted with
    def interacted_id_promoter(self, csv_file_name):

        promoter_data = pd.read_csv(csv_file_name, header=None)
        promoter_data.columns = list(promoter_data.loc[0])
        print(list(promoter_data.loc[0]))
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

        print(len(set(promoter_dict)))
        return promoter_dict

    def draw_bar_chart(self, m, n, name, type):
        if type == 2:
            objects = ("Promoters related bins \n that have interactions with \n" + name,
                       "non-promoters related bins \n that have interactions with \n" + name)
        if type == 1:
            objects = (name + "s that overlap with \n promoters related \n bins ",
                       name + "s that overlap with \n non-promoter related \n bins ")

        y_pos = np.arange(len(objects))
        performance = [m * 100, n * 100]
        plt.clf()
        bar_list = plt.bar(y_pos, performance, align='center', alpha=0.5, width=0.2)
        bar_list[0].set_color('c')
        plt.xticks(y_pos, objects)
        plt.ylabel('interactions percentage')
        plt.xlim(-0.6, 1 + .6)
        plt.title(name + ' and promoters\' interactions')
        plt.savefig("chart_v3_promoter" + str(type) + "/" + name + ".png")
        # plt.show()

    def read_enhancer(self):

        csv_col = ["id", "chr", "start", "end", "element"]
        i = 0
        with open("enhancer/enhancers.csv", 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_col)
            writer.writeheader()

            with open("enhancer/wgEncodeBroadHmmHmecHMM_hg19.txt", 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if i == 0:
                        i += 1
                        continue
                    line = line.strip().split()
                    dict = {}
                    dict["id"] = str(i)
                    dict["chr"] = str(line[1])
                    dict["start"] = str(line[2])
                    dict["end"] = str(line[3])
                    dict["element"] = str(line[4])
                    # dict["type_of_repeat"] = str(line[10])
                    writer.writerow(dict)
                    i += 1

    def split_enhancer(self):
        enhancer_type_list = ["Strong_Enhancer", "Insulator", "Repetitive/CNV", "Heterochrom/lo",
                              "Repressed"]

        # enhancer = ["Strong_Enhancer", "Weak_Enhancer"]
        enhancer_data = pd.read_csv("enhancer/enhancers.csv", header=None)
        enhancer_data.columns = ["id", "chr", "start", "end", "element"]

        temp = enhancer_data.query('element == "Strong_Enhancer" or element == "Weak_Enhancer"')
        temp = temp.reset_index(drop=True)
        temp.drop_duplicates(keep="last", inplace=True)
        temp = temp.reset_index(drop=True)
        temp.to_csv("enhancer/enhancer_file.csv", header=True, index=False)

        for i in enhancer_type_list:
            temp = enhancer_data.query('element == @i')
            temp = temp.reset_index(drop=True)
            temp.drop_duplicates(keep="last", inplace=True)
            temp = temp.reset_index(drop=True)
            if i == "Repetitive/CNV": i = "Repetitive"
            if i == "Heterochrom/lo": i = "Heterochrom"
            temp.to_csv("enhancer/" + i + "_file.csv", header=True, index=False)

    def read(self):
        element_data = pd.read_csv("repeat/interaction_overlap_LTR_repeatMaskers.csv",
                                   header=None)
        element_data.columns = ['1', 'repeat_id', '3', '4', '5', '6', '7', '8', '9']
        print(len(element_data.repeat_id.unique()))
        print(len(element_data))

    def find_specific_location(self):

        # element_data = pd.read_csv("enhancer/Heterochrom_file.csv")
        # element_data = pd.read_csv("enhancer/Repetitive_file.csv")
        element_data = pd.read_csv("enhancer/Repressed_file.csv")
        element_data.columns = list(element_data)
        element_data = element_data.astype({'start': 'int64'})
        element_data = element_data.astype({'end': 'int64'})
        element_data = element_data.sort_values(by=['start'])
        name_list = ["promoter/promoters_with_chr_v3.csv"]
        for name in name_list:
            enhancers_data = pd.read_csv(name)
            enhancers_data.columns = list(enhancers_data)
            enhancers_data = enhancers_data.astype({'start': 'int64'})
            enhancers_data = enhancers_data.astype({'end': 'int64'})
            self.enhancers_data_seprated = []

            for chr in range(1, 24):
                s = str(chr)
                temp = enhancers_data.query("chr == @chr or chr == @s")
                temp = temp.reset_index(drop=True)
                temp = temp.sort_values(by=['start'])
                temp = temp.reset_index(drop=True)


                self.enhancers_data_seprated.append(temp)

            print(len(enhancers_data))
            print(len(self.enhancers_data_seprated))
            print(sum([len(i) for i in self.enhancers_data_seprated]))
            self.m = 0
            v = np.vectorize(self.check_specific_location)
            v(element_data.loc[:]["id"], element_data.loc[:]["start"], element_data.loc[:]["end"],
              element_data.loc[:]["chr"])

            print("res:", len(self.data_hetrochromatin_overlap))

        ids = set(element_data.loc[:]["id"])
        res_ids = set(self.data_hetrochromatin_overlap)
        unique = ids - res_ids
        print(len(unique))
        new_df = element_data.loc[element_data['id'].isin(unique)]
        new_df.to_csv("enhancer/new_Repressed_v3.csv",index=False)

    def make_bin_dictionary_interactions(self):

        bins_dicts = {}
        unique = []
        temp = []
        interaction_file_name = "HMEC_SRR1658680_5k_0.01_5reads_2.05.19.txt"

        interaction_data = pd.read_csv(interaction_file_name, sep="\t", header=None)
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
            if int(bin) not in bins_dicts:
                unique.append(
                    {"id": int(F1_bin[n]),
                     "start": int(F1_start[n]),
                     "end": int(F1_end[n]), "chr": int(F1_chr[n])})
                bins_dicts[int(bin)] = [int(F2_bin[n])]
            else:
                bins_dicts[int(bin)].append(int(F2_bin[n]))
        for n, bin in enumerate(F2_bin):
            temp.append(bin)
            if n == 0:
                continue
            if int(bin) not in bins_dicts:
                unique.append(
                    {"id": int(F2_bin[n]),
                     "start": int(F2_start[n]),
                     "end": int(F2_end[n]), "chr": int(F2_chr[n])})
                bins_dicts[int(bin)] = [int(F1_bin[n])]
            else:
                bins_dicts[int(bin)].append(int(F1_bin[n]))
        unique = sorted(unique, key=lambda i: i['id'])
        print(len(unique))
        csv_columns = ["id","start","end","chr"]
        try:
            with open("bin_data/new_unique_bin.csv", 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in unique:
                    writer.writerow(data)
        except IOError:
            print("I/O error / can not write on csv file")
        # print(bins_dicts)

        # return bins_dicts


    def check_specific_location(self, id, start, end, chr):

        self.m += 1
        # print(self.m)
        if self.m % 10000 == 0:
            print(self.m)
        if chr == "X" or chr == "Y":
            chr = 23
        try:
            chr = int(chr)
        except:
            chr = 24
        if chr < 24:

            data = self.enhancers_data_seprated[int(chr) - 1]
            bin_start = data.loc[:]["start"]
            bin_end = data.loc[:]["end"]
            # bin_id = data.loc[:]["id"]
            length = end - start - 10
            start_index = np.searchsorted(bin_end, start - length)
            end_index = np.searchsorted(bin_start, end)
            i = start_index
            if start_index > 0:
                if bin_end[start_index - 1] > start:
                    print("ridiiiiiiiiiiiiiiiiiiiiiiii")
            if start_index >= 10:
                i = i - 10
            while i <= end_index:
                if end_index >= len(bin_end):
                    break
                if start > bin_start[i] - length and end < bin_end[i] + length:
                    if id not in self.data_hetrochromatin_overlap:
                        self.data_hetrochromatin_overlap.append(id)
                i += 1

                # len(self.data_hetrochromatin_overlap)


p = Hic_proj_overlap()
# p.interacted_id_promoter()
# p.count_promoters("miRNA")
#["Strong_Enhancer","new_Heterochrom_v3","Insulator","new_Repetitive_v3","new_Repressed_v3"]
# ["miRNA","lincRNA","snRNA","rRNA","pseudogene"]
# ["enhancer","Strong_Enhancer","Heterochrom","Insulator","Repressed","Repetitive"]

for i in ["lncRNA"]:
    # p.find_RNAs_overlp(i, is_RNA = True,is_repeat= False, is_enhancer=False, is_exon=False)
    p.promoter_element_interaction_annotation(i,"promoter/new_overlap_promoter_version2.csv", is_RNA=True, is_repeat=False)
# p.read()
# p.read_enhancer()
# p.split_enhancer()
# p.read_repeat_element()
# p.find_interaction_starts()
# p.split()

# print("all finished")
# # p.interacted_id_promoter()
# file_names = ["lincRNA"]
# for i in file_names:
#     p.promoter_element_interaction_annotation(i, is_repeat=False, is_RNA=True)

# for name in ["repeat"]:
#     m, n = p.count_promoters(name)
#     p.draw_bar_chart(n, m, name)
# p.draw_bar_chart(0.46,0.54,"Simple_repeat")
# p.read_repeat_element()
# p.seprate_repeat_data()
# print("finish")
# p.seprate_repeat_data()
# p.draw_bar_chart(0.54,0.46,"LINE repeat",1)
# p.find_interactions_for_bins()
# p.find_specific_location()
# p.check_specific_location(1,10,10)
# p.make_bin_dictionary_interactions()
# p.draw_bar_chart(0.93,0.87,"SINE_repeatMaskers",2)