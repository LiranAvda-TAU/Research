from Code.Enum.FileType import FileType
from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileMaker import FileMaker
from Code.Files.FileReader import FileReader
from Code.Utils.BioPython import BioPython
from Test.Extraction.DataExtracterTest import DataExtracterTest
from Test.Files.FileReaderTest import FileReaderTest
from Test.TestFunctions import TestFunctions


class FileMakerTest:

    @staticmethod
    def fromSetToFileTest():
        det = DataExtracterTest()
        s = det.getConditionsSetTest()

        fm = FileMaker()
        fm.fromSetToFile(s, "conditions")

    @staticmethod
    def fromDictToFileTest():
        det = DataExtracterTest()
        d = det.findGenesWithValidConditionAndHomologTest()

        fm = FileMaker()
        fm.fromDictToUniqueFile(d, "genes, ortholougs and conditions")

    @staticmethod
    def fromTwoDictToFileTest():
        det = DataExtracterTest()
        data = det.findGenesWithValidConditionAndHomologTest()

        fd_humans = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data", r"\gene_id-gene_name-start-end.txt",
                               FileType.TSV)
        humanGenesLength = fd_humans.getGenesLength(1, 2, 3)
        fd_c_elegans = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data", r"\c_elegans_genes-start-end.txt",
                                  FileType.TSV)
        cElegansGenesLength = fd_c_elegans.getGenesLength(0, 1, 2)

        de = DataExtracter()
        orthologs, phenotypes = de.filter_genes_with_size_differences(data, humanGenesLength, cElegansGenesLength, p=10)
        print(str(len(orthologs)))
        print(str(len(phenotypes)))

        fm = FileMaker()
        fm.fromTwoDictToFile(orthologs, phenotypes, "genes-homolougs-phenotypes-filteredBySize")

    @staticmethod
    def fromDictToFileEasyTest():
        # frt = FileReaderTest()
        # print("extracting genes...\n")
        # genesList = frt.getCElegansGenesFromOutputTest()
        #
        # bp = BioPython()
        # print("getting genes' id...\n")
        # gisDict: dict = bp.getGenesId(genesList)
        # print("all genes' id are gotten, there are " + str(len(gisDict)) + " genes with id\n")
        #
        # testFunc = TestFunctions("getGenesId", dictionary=gisDict)
        # testFunc.printRandomLinesInDict(5)
        #
        # fm = FileMaker()
        # fm.fromDictToFile(gisDict, "relevant-c-elegans-genes-and-their-id")

        # det = DataExtracterTest()
        # genesAndAccessions = det.fromMultipleAccessionsToOneTest()
        #
        # fm = FileMaker()
        # fm.fromDictToFile(genesAndAccessions, "c-elegans-genes-and-longest-accession_number")

        fr = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Test\Files",
                                             r"\c-elegans-genes-and-longest-accession_number",
                                             FileType.TSV)
        genesAndAccessions = fr.fromFileToDict(0, 1)

        print("Done extracting file genesAndAccessions to a dictionary")
        TestFunctions("fromFileToDict", dictionary=genesAndAccessions).printFirstLinesInDict(1)

        accessions = list(genesAndAccessions.values())

        bp = BioPython()
        accessionsAndSeqs = bp.make_accession_number_and_seq_dict(accessions, "ORIGIN", "translation=")

        print("Done creating a dictionary for accession numbers and sequences")
        TestFunctions("make_accession_number_and_seq_dict", dictionary=accessionsAndSeqs).printFirstLinesInDict(1)

        fm = FileMaker()
        fm.fromDictToFile(accessionsAndSeqs, "accession-numbers-and-sequences")

    @staticmethod
    def fromOneFileToManyTest():

        fm = FileMaker()
        fm.fromOneFileToMany(r"C:\Users\Liran\PycharmProjects\Research\Test\Files\accession-numbers-and-sequences",
                             r"accessions-and-sequences\accessions-and-sequences-part-",
                             25)

    @staticmethod
    def fixTabbedFileTest():

        fm = FileMaker()
        fm.fixTabbedFile(
            r"C:\Users\Liran\PycharmProjects\Research\Executors\extra_genes_id_number_and_chosen_accessions-120619",
            r"C:\Users\Liran\PycharmProjects\Research\Executors\extra_genes_ids_number_and_chosen_accessions-120619")



### TESTING ###

# fmt = FileMakerTest()

# fmt.fromSetToFileTest()
# fmt.fromDictToFileTest()
# fmt.fromTwoDictToFileTest()
# FileMakerTest.fromDictToFileEasyTest()

# FileMakerTest.fromOneFileToManyTest()

# FileMakerTest.fixTabbedFileTest()

# FileMaker().separateOrthologsInfoInData(r"C:\Users\Liran\PycharmProjects\Research\Test\Files\filterized_data_and_gene_summary-190619",
#                                         "filterized-data-summarized-190619")
