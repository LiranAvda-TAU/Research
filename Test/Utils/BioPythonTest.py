from Code.Enum.FileType import FileType
from Code.Files.FileReader import FileReader
from Code.Utils.BioPython import BioPython
from Test.TestFunctions import TestFunctions


class BioPythonTest:

    @staticmethod
    def blastForAllTest():
        fr = FileReader(FileReader.research_path + r"\Test\Files",
                        r"\c-elegans-genes-and-longest-accession_number",
                        FileType.TSV)
        genesAndAccessions = fr.fromFileToDict(0, 1)

        TestFunctions("fromFileToDict", dictionary=genesAndAccessions).printFirstLinesInDict(1)

        bp = BioPython()
        bp.blast_for_all("blastp", "nr", genesAndAccessions)

    @staticmethod
    def blastpForAllTest():
        fr = FileReader(FileReader.research_path + r"\Test\Files",
                        r"\c-elegans-genes-and-longest-accession_number",
                        FileType.TSV)
        genesAndAccessions = fr.fromFileToDict(0, 1)

        TestFunctions("fromFileToDict", dictionary=genesAndAccessions).printFirstLinesInDict(1)

        bp = BioPython()
        bp.blastp_for_all_outdated("blastp", "nr", genesAndAccessions, end_word="ORIGIN", translation_word="translation=")

    @staticmethod
    def blastpForAll_improvedTest():
        fr = FileReader(FileReader.research_path + r"\Test\Files\accessions-and-sequences",
                        r"\accessions-and-sequences-part-18",
                        FileType.TSV)
        accessionsAndSequences = fr.fromFileToDict(0, 1)

        print("accessions and sequences dictionary is achieved, and has " + str(len(accessionsAndSequences)) + " items")
        TestFunctions("fromFileToDict", dictionary=accessionsAndSequences).printFirstLinesInDict(1)

        bp = BioPython()
        bp.blastp_by_accessions("blastp", "nr", accessionsAndSequences, end_word="ORIGIN",
                                translation_word="translation=")


### TESTING ###

# BioPythonTest.blastForAllTest()
# BioPythonTest.blastpForAll_improvedTest()
