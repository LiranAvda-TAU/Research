from Code.Enum.FileType import FileType
from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileMaker import FileMaker
from Code.Files.FileReader import FileReader
from Code.Utils.BioPython import BioPython
from Test.Files.FileReaderTest import FileReaderTest
from Test.TestFunctions import TestFunctions


class DataExtracterTest:
    @staticmethod
    def findGenesVariantsAndHomologousTest():
        filePath = FileReader.research_path + r"\Data"
        fileName = "\human-c.elegans.txt"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        orthologousGenes = fd.read_genes_with_orthology_confidence()

        filePath = FileReader.research_path + r"\Data\GenesVariants"
        fileName = r"-human-c.elegans-variant-"
        fileType = FileType.TSV

        fd = FileReader(filePath, fileName, fileType)
        genesAndVariants = fd.read_all_genes_with_variants()

        de = DataExtracter().find_genes_variants_and_homologous(orthologousGenes, genesAndVariants)
        return de

    @staticmethod
    def getConditionsSetTest():
        frt = FileReaderTest()
        genesWithoutNonValidConditions = frt.readFileFilterConditionsTest()

        conditionsSet = DataExtracter().get_conditions_set(genesWithoutNonValidConditions, 0)
        print("set contains " + str(len(conditionsSet)) + " different conditions")
        return conditionsSet

    @staticmethod
    def find_genes_with_valid_condition_and_homolog_test():
        frt = FileReaderTest()
        genesWithoutNonValidConditions = frt.readFileFilterConditionsTest()
        humanGenes, homologousGenes = frt.readOrtholist()

        genesWithHomolougAndValidCondition = DataExtracter()\
            .find_genes_with_valid_condition_and_homolog(genesWithoutNonValidConditions, homologousGenes)
        return genesWithHomolougAndValidCondition

    @staticmethod
    def fromMultipleAccessionsToOneTest():
        frt = FileReaderTest()
        genesAndAccessions = frt.fromFileToDictWithPluralValuesTest()

        genesWithChosenAccessions = DataExtracter.from_multiple_accessions_to_one(genesAndAccessions)
        return genesWithChosenAccessions

    @staticmethod
    def getConservedDomainsTest(file_path, file_name, file_type):
        fr = FileReader(file_path,
                        file_name,
                        file_type)
        genes_list = fr.get_genes_list()

        de = DataExtracter()
        genesAndConservedDomains = de.get_conserved_domains(genes_list, "Conserved Domains (")
        return genesAndConservedDomains

    @staticmethod
    def checkReverseHomologyTest(num_of_files, file_path, file_name):
        de = DataExtracter()
        accession_and_ids, ids_and_hsp = de.get_hit_ids_for_accession_numbers(num_of_files, file_path, file_name)
        return accession_and_ids, ids_and_hsp

    @staticmethod
    def checkInitializer():
        de = DataExtracter()
        print(len(de.c_elegans_name_id_converter))
        print(de.c_elegans_name_id_converter["hmp-2"])
        print(len(de.human_name_id_converter))
        print(de.human_name_id_converter["CTNNB1"])

    # returns a file that consists of C.elegans genes id and the longest accession number
    @staticmethod
    def get_c_elegans_longest_accession_numbers(to_file: bool = False):
        genes_and_many_accessions = \
            FileReader(FileReader.research_path + r"\Data",
                       r"\all-genes-ids-and-accession-numbers.txt").from_file_to_dict_with_plural_values(key_index=0,
                                                                                                         value_index=1,
                                                                                                         delete_first=True)

        genes_and_accessions = DataExtracter().from_multiple_accessions_to_one(genes_and_many_accessions)

        print("number of genes id with accessions: ", str(len(genes_and_accessions)))
        if to_file:
            FileMaker().from_dict_to_file(genes_and_accessions, "c-elegans-gene-ids-and-accession-numbers")
        return genes_and_accessions

    # the function read the most updates accession file (120919) and uses a function to extract the sequence for each
    # accession number. it returns a dictionary of accession numbers as keys and sequences as values
    @staticmethod
    def get_c_elegans_sequences_for_accession_numbers(to_file: bool = False):
        accessions = list(FileReader(FileReader.research_path + r"\Data",
                                     r"\c-elegans-gene-ids-and-accession-numbers").from_file_to_dict().values())

        bp = BioPython()
        accessions_and_seqs = bp.make_accession_number_and_seq_dict(accessions, "ORIGIN", "translation=")

        print("Done creating a dictionary for accession numbers and sequences")
        TestFunctions("make_accession_number_and_seq_dict", dictionary=accessions_and_seqs).print_first_lines_in_dict(1)

        if to_file:
            fm = FileMaker()
            fm.from_dict_to_file(accessions_and_seqs, "c-elegans-accession-numbers-and-sequences")

# TESTING

det = DataExtracterTest()

# testFunc = TestFunctions("find_genes_variants_and_homologous", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.printFirstLines(5)
# testFunc.printRandomLines(5)
# testFunc.testFinisher()

# det.getConditionsSetTest()

# d: dict = det.find_genes_with_valid_condition_and_homolog_test()
# testFunc = TestFunctions("find_genes_with_valid_condition_and_homolog", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()


# d: dict = det.fromMultipleAccessionsToOneTest()
# testFunc = TestFunctions("from_multiple_accessions_to_one", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()

# d: dict = det.getConservedDomainsTest()
# testFunc = TestFunctions("get_conserved_domains", d)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()

# d1, d2 = det.checkReverseHomologyTest(1,
#                                       FileReader.research_path + r"\Test\Files\accessions-and-sequences",
#                                       r"\results-part-")
# testFunc = TestFunctions("get_hit_ids_for_accession_numbers - accessions and hit ids", d1)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()
#
# testFunc = TestFunctions("get_hit_ids_for_accession_numbers - hit ids and hsp values", d2)
# testFunc.checkSize()
# testFunc.checkBooleanResult()
# testFunc.print_first_lines_in_dict(5)
# testFunc.printRandomLinesInDict(5)
# testFunc.testFinisher()

# DataExtracter().addDescriptionFromWormBase(FileReader.research_path + r"\Test\Files\filterized-data-"
#                                            r"summarized-190619",
#                                            "filterized-data-with-all-summary-190619",
#                                            6,
#                                            "Not mentioned",
#                                            2)

# det.checkInitializer()
#
# det.get_c_elegans_longest_accession_numbers(True)
# det.get_c_elegans_sequences_for_accession_numbers(True)
