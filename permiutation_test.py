from gtfparse import read_gtf

class Permutation():

    def __init__(self,gene_gtf_file_name):
        self.gtf_data = None
        self.chr_length_dict = {}
        self.gene_gtf_file_name = gene_gtf_file_name

    def read_data(self):
        self.gtf_data = read_gtf(self.gene_gtf_file_name,
                              usecols=['seqname', 'start'])

        for i in range(1,24) :
            if i == 23 :
                temp = self.gtf_data.query("chr = X")
                self.chr_length_dict['X'] = max(temp.loc[:]["start"])
                temp = self.gtf_data.query("chr = Y")
                self.chr_length_dict['Y'] = max(temp.loc[:]["start"])

            temp = self.gtf_data.query("chr = @i")
            self.chr_length_dict[i] = max(temp.loc[:]["start"])

        print(self.chr_length_dict)


proj = Permutation("Homo_sapiens.GRCh37.87.gtf")
proj.read_data()