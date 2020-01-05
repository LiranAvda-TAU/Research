import xlrd
from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType
import unicodedata
import pandas as pd

chromosomes = [str(i) for i in range(1, 23)] + ['x', 'y']
fileTypesDelimeter = {FileType.TSV: "\t", FileType.CSV: ",", FileType.UNCLEAR: "\" \""}


class FileReader:

    research_path = r"C:\Users\Liran\PycharmProjects\Research"

    def __init__(self, path, fileName, fileType: FileType = FileType.TSV):
        self.path = path
        self.name = fileName
        self.type = fileType
        self.size = -1  # not defined yet

    # extracts all genes with confidence in orthology to C.elegans's genes
    # returns a dictionary with key: human gene id, value: c.elegans gene id
    def readGenesWithOrthologyConfidence(self):
        orthologousGenes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        print("sanity check: " + f.readline())  # sanity check
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])

            if line[10] == '1':  # confidence in orthology
                orthologousGenes[line[0]] = line[4]  # genes
        f.close()
        return orthologousGenes

    # receiving an open file for reading, extracts all genes with alleles and puts them in the dictionary
    # the dictionary contains key: gene id, value:chromosome, start, end, variantName, original bp, alleles, place
    @staticmethod
    def __readGenesWithVariants(file, genesAndVariants, variantType):
        file.readline()  # skip the headlines
        for row in file:
            line = row.rstrip('\n').split('\t')

            if line[1] != "":  # variant exists
                pass
                variantAlleles = line[5]
                if variantAlleles == "COSMIC_MUTATION":
                    alleles = ["cosmic", "cosmic"]
                else:
                    alleles = variantAlleles.split("/")
                    # gene info = chromosome, start, end, original bp, alleles, place, somatic/germ line
                    geneInfo = [line[2], line[3], line[4], alleles[0], alleles[1:], line[7], variantType]
                    genesAndVariants[(line[0], line[1])] = geneInfo

    # receives the path and goes through all number of chromosomes
    def readAllGenesWithVariants(self):
        genesAndVariants = {}
        variantTypes = ["somatic", "germline"]
        for variantType in variantTypes:
            for chromosomeNumber in chromosomes:
                filePath = self.path + "\\" + variantType + self.name + str(chromosomeNumber) + ".txt"
                try:
                    file = open(filePath, FileMode.READ.value)
                    self.__readGenesWithVariants(file, genesAndVariants, variantType)
                except IOError:
                    print("Oops, File " + filePath + " does not exist :O")
        return genesAndVariants

    # half conditions - conditions that are to be filtered only when they are alone
    def readFileFilterConditions(self, condition_column, conditions):
        genes = {}
        if self.type == FileType.XLS:
            wb = xlrd.open_workbook(self.path + self.name)
            sheet = wb.sheet_by_index(0)
        else:
            print("not supported yet")
            exit()
        for row in range(1, sheet.nrows):
            gene_conditions = sheet.cell_value(row, condition_column)
            invalid_conditions = 0
            listed_gene_conditions = FileReader.\
                listCleaner(gene_conditions.replace("?", ",").replace("|", ",").split(","))
            for gene_condition in listed_gene_conditions:
                for condition in conditions:
                    if condition in gene_condition or condition.lower() in gene_condition:
                        invalid_conditions += 1
                        if "modifier" in condition or "susceptibility" in condition or "risk" in condition:
                            invalid_conditions += 1
            if invalid_conditions < len(listed_gene_conditions):
                genes[sheet.cell_value(row, 0)] = sheet.row_values(row, 1)  # column 0 is the gene name
        return genes

    # receives a list and returns the list without items that are empty or made of spaces
    @staticmethod
    def listCleaner(l: list):
        new_list = []
        for item in l:
            if item.strip():
                new_list.append(item)
        return new_list

    # read the ortholist excel page, and returns two dicts: one converting ensembleId of a human gene to its
    # HGNC symbol, and one converting from human gene HGNC symbol to its c-elegans ortholog and number of programs
    # supporting that claim

    def readOrtholist(self):
        humanGenes = {}
        homologousGenes = {}
        wb = xlrd.open_workbook(self.path + self.name)
        sheet = wb.sheet_by_index(0)
        for row in range(1, sheet.nrows):
            HGNCSymbol = sheet.cell_value(row, 5)
            EnsemblId = sheet.cell_value(row, 4)
            if EnsemblId not in humanGenes:
                humanGenes[EnsemblId] = HGNCSymbol  # Ensembl ID
                homologousGenes[HGNCSymbol] = []
            # 0 - WormBase ID, 6 - number of programs out of 6
            homologousGenes[HGNCSymbol].append(tuple([sheet.cell_value(row, 0), int(sheet.cell_value(row, 6))]))
        return humanGenes, homologousGenes

    # receives a table of human genes and returns a dict that converts the human gene name to its id
    def genesIdAndNames(self):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        f.readline()  # headlines
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            genes[line[1]] = line[0]  # geneName = geneId
        f.close()
        return genes

    # returns a dict with gene names as keys and their length in bp as values
    def getGenesLength(self, gene_name_index, start_index, end_index):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        f.readline()  # headlines
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            # line[3] - gene end, line[2] = gene start
            genes[line[gene_name_index]] = int(line[end_index])-int(line[start_index])
        f.close()
        return genes

    # returns a list of all c.elegans genes
    def get_genes_list(self, column: int = 0):
        genes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            genes.append(row.strip("\n").split(fileTypesDelimeter[self.type])[column])
        f.close()
        return genes

    def get_c_elegans_genes_from_output(self, column):
        genes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            genes.extend(FileReader.fromStringToGenesList(line[column], 0))
        f.close()
        return genes

    # receives (1) a key index, (2) value index, (3) boolean value determining whether to disregard the first sentence,
    # and returns a dictionary of all genes and the longest coding sequence length
    def get_genes_cd_length(self, key_index, value_index, delete_first_line: bool = False):
        lengths = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            try:
                gene = line[key_index]
                length = line[value_index]
                if gene in lengths:
                    if length > lengths[gene]:
                        lengths[gene] = length
                else:
                    lengths[gene] = length
            except:
                print(line)
                f.close()
                exit()
        f.close()
        return lengths

    def from_file_to_dict(self, key_index: int = 0, value_index: int = 1, delete_first_line: bool = False):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            try:
                genes[line[key_index]] = line[value_index]
            except:
                print("An error has occured!")
                print(line)
                f.close()
                exit()
        f.close()
        return genes

    def fromFileToDictWithTupleKey(self, first_key_index: int = 0, second_key_index: int = 1, value_index: int = 2,
                                   delete_first_line: bool = False):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            try:
                genes[(line[first_key_index], line[second_key_index])] = line[value_index]
            except:
                print("An error has occured!")
                print(line)
                f.close()
                exit()
        f.close()
        return genes

    def fromFileToList(self, delete_first_line: bool = False):
        keys = []
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            keys.append(row.rstrip('\n'))
        return keys

    def fromFileToDictWithPluralValues(self, key_index, value_index, delete_first: bool = False):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            if line[key_index] not in genes:
                genes[line[key_index]] = {line[value_index]}
            else:
                genes[line[key_index]].add(line[value_index])
        f.close()
        return genes

    def from_MMP_file_to_dict_with_listed_values(self, key_index, list_of_value_indexes: list, delete_first: bool = False):
        dic = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            if line[7] == "intergenic":
                # not in a gene
                continue
            values = []
            for value_index in list_of_value_indexes:
                values.append(line[value_index])
            dic[line[key_index].strip(" ")] = values
        f.close()
        return dic

    def fromFileToDictWithListsAsValues(self, key_index, value_index, delete_first: bool = False):
        dic = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            key = line[key_index]
            values = line[value_index].strip("\n").split(", ")
            dic[key] = values
        f.close()
        return dic

    def readResultsFile(self):
        try:
            f = open(self.path + self.name, FileMode.READ.value)
        except:
            print("file cannot be opened, maybe it doesn't exist")
            exit()
        return f

    def getCelegansIdToHumanNameDict(self, value_index: int = 1, key_index: int = 2):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            cElegansGenesTuples = FileReader.fromStringToTuplesList(line[key_index])
            for t in cElegansGenesTuples:
                genes[t[0]] = line[value_index]
        f.close()
        return genes

    def sizeOfFile(self):
        count = 0
        f = open(self.name, FileMode.READ.value)
        for _ in f:
            count += 1
        f.close()
        self.size = count
        return count

    @staticmethod
    def fromStringToTuplesList(s: str):
        tuples_list = []
        str_list = s.replace("(", "").replace("'", "").replace(" ", "").split("),")
        for item in str_list:
            item = item.replace(")", "")
            tuples_list.append(tuple(item.split(",")))
        return tuples_list

    @staticmethod
    def fromStringToGenesList(s: str, c_elegans_column: int):
        genes = []
        tuples_list = FileReader.fromStringToTuplesList(s)
        for t in tuples_list:
            genes.append(t[c_elegans_column])
        return genes

    def readRedundantFile(self, key_index, value_index, delete_first_line: bool = False):
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            if line[key_index] not in genes:
                genes[line[key_index]] = line[value_index]
        f.close()
        return genes

    def fromPluralFileToOneFile(self, num_of_files, target_file_path):
        target_file = open(target_file_path, FileMode.WRITE.value)
        for i in range(num_of_files):
            f = open(self.path+self.name+str(i))
            for row in f:
                target_file.write(row)
            print("File number " + str(i) + " is now completely downloaded")
            f.close()
        target_file.close()

    def extractCelegansFromData(self, orthologs_column, gene_column):
        genes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            c_elegans_part = line[orthologs_column]
            genes.extend(FileReader.fromStringToGenesList(c_elegans_part, gene_column))
        f.close()
        return genes

    def makeDictFromSummary(self, key_index, value_index, delete_first_line: bool = False):
        genesDescription = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            try:
                genesDescription[line[key_index].strip("\"")] = line[value_index]
            except:
                print(line)
                f.close()
                exit()
        f.close()
        return genesDescription

    def addSummaryToGenesAndSeparateTuples(self, summary_path, new_data_plus_summary, key_index, summary_index,
                                           summary_delimeter, delete_first_line: bool = False):
        summary = open(summary_path)
        summary_dic = FileReader(summary_path, "", summary_delimeter).makeDictFromSummary(key_index, summary_index,
                                                                                          delete_first_line)
        summary.close()

        data = open(self.path + self.name, FileMode.READ.value)
        new_data = open(new_data_plus_summary, FileMode.WRITE.value)

        for row in data:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            human_gene_id = line[0]
            human_gene_name = line[1]
            orthologs = FileReader.fromStringToTuplesList(line[2])
            phenotype = line[3]

            for t in orthologs:
                ortholog_id = t[0]
                if ortholog_id in summary_dic and len(summary_dic[ortholog_id]) > 1:
                    summary_value = summary_dic[ortholog_id]
                else:
                    summary_value = "Not mentioned"
                new_data.write(human_gene_id + "\t" + human_gene_name + "\t" + ", ".join(t) + "\t" +
                                phenotype + "\t" + summary_value + "\n")
        data.close()
        new_data.close()

    @staticmethod
    def lineFixer(line):
        new_line = ''
        line = unicodedata.normalize("NFKD", line)
        for ch in line:
            if ('A' <= ch <= 'z') or ch == " " or '1' <= ch <= '9':
                new_line += ch
        return new_line

    def fromHGMDtoDict(self):
        genes_and_variants = {}
        f = open(self.path + self.name, FileMode.READ.value)
        line = f.readline().rstrip("\n")
        while line:
            if line.endswith(":"):
                gene = line.strip("\n")[:line.find(":")]
                genes_and_variants[gene] = []
                line = f.readline()
                while not line.startswith("\n") and len(line) > 3:
                    variant = line.strip("\n").split("\t")[1]
                    genes_and_variants[gene].append(variant)
                    line = f.readline()
                line = f.readline().rstrip("\n")
        f.close()
        return genes_and_variants

    def readData(self, key_index, delete_first_line: bool = False):
        genes_data = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(fileTypesDelimeter[self.type])
            key = line[key_index]
            genes_data[key] = row.strip("\n")
        f.close()
        return genes_data

    def get_list_from_excel_using_pandas(self, column_name='WormBase Gene ID', sheet_name='kinase'):
        xls = pd.ExcelFile(self.path + self.name)
        df = pd.read_excel(xls, sheet_name)
        return [gene_id for gene_id in df[column_name]]

    def get_dictionary_from_excel_using_pandas(self, key_column_name, value_column_name, sheet_name='kinase'):
        xls = pd.ExcelFile(self.path + self.name)
        df = pd.read_excel(xls, sheet_name)
        two_list = [(key, value) for key, value in zip(df[key_column_name], df[value_column_name])]
        return dict(two_list)






