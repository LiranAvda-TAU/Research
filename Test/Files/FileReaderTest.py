from Code.Files.FileReader import FileReader
from Code.Enum.FileType import FileType
from Test.TestFunctions import TestFunctions


class FileReaderTest:

    @staticmethod
    def readGenesWithOrthologyConfidenceTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = "\human-c.elegans.txt"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        orthologousGenes = fd.readGenesWithOrthologyConfidence()
        return orthologousGenes

    @staticmethod
    def readAllGenesWithVariantsTest():
        filePath = FileReader.research_path + r"\Data\GenesVariants"
        fileName = r"-human-c.elegans-variant-"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genesAndVariants = fd.readAllGenesWithVariants()
        return genesAndVariants

    @staticmethod
    def readFileFilterConditionsTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = r"\maybe DM HGMD.xlsx"
        fileType = FileType.XLS
        unwanted_conditions = ["Autism", "Schizophrenia", "Bipolar", "Obsessive",
                               "Alzheimer", "Parkinson", "Schizoaffective", "Hyperactivity", "ADHD",
                               "Amyotrophic lateral sclerosis", "Tourette", "X-inactivation", "Skewed",
                               "Speech delay", "Premature ovarian failure", "Primary ovarian insufficiency",
                               "Prion disease", "Miscarriage", "Cancer", "Late onset", "modifier",
                               "susceptibility", "risk", "Intellectual disability"]

        fd = FileReader(filePath, fileName, fileType)
        genesWithoutConditions = fd.readFileFilterConditions(1, unwanted_conditions)
        return genesWithoutConditions

    @staticmethod
    def readOrtholist():
        filePath = FileReader.research_path + r"\Data"
        fileName = r"\ortholist_master.xlsx"
        fileType = FileType.XLS

        fd = FileReader(filePath, fileName, fileType)
        humanGenes, homologousGenes = fd.readOrtholist()
        return humanGenes, homologousGenes

    @staticmethod
    def genesIdAndNamesTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = r"\gene_id-gene_name.txt"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        fromGeneNameToIdDict = fd.genesIdAndNames()
        return fromGeneNameToIdDict

    @staticmethod
    def getCElegansGenesTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = r"\c_elegans_genes-start-end.txt"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genes = fd.get_genes_list()
        return genes

    @staticmethod
    def getCElegansGenesFromOutputTest():
        filePath = FileReader.research_path + r"\Test\Files"
        fileName = r"\genes-homolougs-phenotypes-filteredBySize"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genes = fd.get_c_elegans_genes_from_output(column=1)
        return genes

    @staticmethod
    def getGenesLengthTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = r"\gene_id-gene_name-start-end.txt"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genesLength = fd.getGenesLength(1,2,3)
        return genesLength

    @staticmethod
    def fromFileToDictTest():
        filePath = FileReader.research_path + r"\Test\Files"
        fileName = r"\relevant-c-elegans-genes-and-their-id"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genes = fd.from_file_to_dict(key_index=0, value_index=1)
        return genes

    @staticmethod
    def fromFileToDictWithPluralValuesTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = r"\gene_id_accession_number.txt"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genes = fd.fromFileToDictWithPluralValues(key_index=0, value_index=1)
        return genes

    @staticmethod
    def fromPluralFileToOneFileTest():
        filePath = FileReader.research_path + r"\Test\Files\accessions-and-sequences"
        fileName = r"\results-part-"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        fd.fromPluralFileToOneFile(19, "blast_results")

    @staticmethod
    def extractCelegansFromDataTest():
        filePath = FileReader.research_path + r"\Executors"
        fileName = r"\filterized-data-by-true-homologs"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genes = fd.extractCelegansFromData(2, 0)
        return genes

    @staticmethod
    def addSummaryToGenesAndSeparateTuplesTest():
        filePath = FileReader.research_path + r"\Executors"
        fileName = r"\filterized-data-190619"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        fd.addSummaryToGenesAndSeparateTuples(FileReader.research_path + r"\Data\summary.csv",
                                              "filterized_data_and_gene_summary-190619",
                                              0,
                                              8,
                                              FileType.UNCLEAR,
                                              True)


# TESTING

frt = FileReaderTest()


# d: dict = frt.readGenesWithOrthologyConfidenceTest()
# testFunc = TestFunctions("readGenesWithOrthologyConfidence", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()
#
# d: dict = frt.readAllGenesWithVariantsTest()
# testFunc = TestFunctions("readGenesWithVariants", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()

# d: dict = frt.readFileFilterConditionsTest()
# testFunc = TestFunctions("readFileFilterConditions", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()
#
# humanGenes, homologousGenes = frt.readOrtholist()
# testFunc = TestFunctions("readOrtholist - humanGenes", humanGenes)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()
#
# testFunc = TestFunctions("readOrtholist - humanGenes", homologousGenes)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()
#
# d: dict = frt.genesIdAndNamesTest()
# testFunc = TestFunctions("genesIdAndNames", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()
#
# d: dict = frt.getGenesLengthTest()
# testFunc = TestFunctions("getGenesLength", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()

# genesList = frt.getCElegansGenesTest()
# testFunc = TestFunctions("getCElegansGenes", lst=genesList)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.printFirstLinesInLst(5)
# testFunc.printRandomLinesInLst(5)


# genesList = frt.getCElegansGenesFromOutputTest()
# testFunc = TestFunctions("get_c_elegans_genes_from_output", lst=genesList)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.printFirstLinesInLst(5)
# testFunc.printRandomLinesInLst(5)

# genesDict = frt.fromFileToDictTest()
# testFunc = TestFunctions("from_file_to_dict", dictionary=genesDict)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)

# genesDict = frt.fromFileToDictWithPluralValuesTest()
# testFunc = TestFunctions("fromFileToDictWithPluralValues", dictionary=genesDict)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)

# frt.fromPluralFileToOneFileTest()

# genesList = frt.extractCelegansFromDataTest()
# testFunc = TestFunctions("extractCelegansFromDataTest", lst=genesList)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.printFirstLinesInLst(5)
# testFunc.printRandomLinesInLst(5)
#
# print(", ".join(genesList))

# frt.addSummaryToGenesAndSeparateTuplesTest()

