import pandas as pd


class eqtm(object):
    def __init__(self, filepath, sep):
        self.filepath = filepath
        self.sep = sep

        self.data = None
        self.bedtoolsFormat_filepath = ""

    def read_eqtmFile(self):
        self.data = pd.read_csv(self.filepath, sep=self.sep)
        return self.data

    def read_eqtm_chr_startEndPos_name(self):
        return [['chr'+item.split(self.sep)[2],
                 int(item.split(self.sep)[3])-25,
                 int(item.split(self.sep)[3])+25,
                 item.split(self.sep)[1]]
                for item in open(self.filepath, "r").readlines()[1:]]

    def save_in_bedtoolsFormat(self, save_filepath):
        self.bedtoolsFormat_filepath = save_filepath
        eqtm_content = self.read_eqtm_chr_startEndPos_name()
        bedtoolsFormat = open(self.bedtoolsFormat_filepath, "w")
        for line in eqtm_content:
            bedtoolsFormat.write("\t".join([str(item) for item in line]))
            bedtoolsFormat.write("\n")
        print("Saved eqtm in bedtools format in path %s" % self.bedtoolsFormat_filepath)

