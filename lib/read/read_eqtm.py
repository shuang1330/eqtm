class eqtm(object):
    def __init__(self, filepath, sep="\t"):
        self.filepath = filepath
        self.sep = sep

        self.headers = []
        self.data = []
        self.bedtoolsFormat_filepath = ""

        self.cond_colName = ""
        self.cond_operator = ""
        self.cond_threshold = 0
        self.cond_colInd = 0

    def read_eqtmFile(self, cond_colName,
                      cond_operator, cond_threshold,
                      full=False):
        f = open(self.filepath, "r")
        self.headers = f.readline().strip().split(self.sep)
        print("Headers:", self.headers)
        self.cond_colName = cond_colName
        self.cond_colInd = self.headers.index(self.cond_colName)
        self.cond_threshold = cond_threshold
        self.cond_operator = cond_operator
        if full:
            self.data = [line.strip().split(self.sep) for line in f.readlines()[1:]]
        else:
            while True:
                line = f.readline()
                if line:
                    content = line.strip().split(self.sep)
                    compare_value = float(content[self.cond_colInd])
                    if self.cond_operator == "smaller_than":
                        if compare_value <= self.cond_threshold:
                            self.data.append(content)
                    elif self.cond_operator == "larger_than":
                        if compare_value >= self.cond_threshold:
                            self.data.append(content)
                    elif self.cond_operator == "equal_to":
                        if compare_value == self.cond_threshold:
                            self.data.append(content)
                    else:
                        raise IOError("Please insert the right operator name: "
                                      "smaller_than, larger_than, or equal_to")
                else:
                    break
        f.close()
        return self.data

    def read_eqtm_chr_startEndPos_name(self):
        if self.data is None:
            raise IOError("Please call read_eqtmFile function first.")
        return [['chr'+item[2], int(item[3])-25,
                 int(item[3])+25, item[1]]
                for item in self.data]

    def save_in_bedtoolsFormat(self, save_filepath):
        self.bedtoolsFormat_filepath = save_filepath
        eqtm_content = self.read_eqtm_chr_startEndPos_name()
        bedtoolsFormat = open(self.bedtoolsFormat_filepath, "w")
        for line in eqtm_content:
            bedtoolsFormat.write("\t".join([str(item) for item in line]))
            bedtoolsFormat.write("\n")
        print("Saved eqtm in bedtools format in path %s" % self.bedtoolsFormat_filepath)

