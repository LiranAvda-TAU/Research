from Code.Enum.FileType import FileType
from Code.Files.FileReader import FileReader
from Code.Utils.BioPython import BioPython
from Test.TestFunctions import TestFunctions


class BioPythonTest:
    bp = BioPython()

    def blastForAllTest(self):
        fr = FileReader(FileReader.research_path + r"\Test\Files",
                        r"\c-elegans-genes-and-longest-accession_number",
                        FileType.TSV)
        genesAndAccessions = fr.from_file_to_dict(0, 1)

        TestFunctions("from_file_to_dict", dictionary=genesAndAccessions).print_first_lines_in_dict(1)

        self.bp.blast_for_all("blastp", "nr", genesAndAccessions)

    def blastpForAllTest(self):
        fr = FileReader(FileReader.research_path + r"\Test\Files",
                        r"\c-elegans-genes-and-longest-accession_number",
                        FileType.TSV)
        genesAndAccessions = fr.from_file_to_dict(0, 1)

        TestFunctions("from_file_to_dict", dictionary=genesAndAccessions).print_first_lines_in_dict(1)

        self.bp.blastp_for_all_outdated("blastp", "nr", genesAndAccessions, end_word="ORIGIN", translation_word="translation=")

    @staticmethod
    def blastpForAll_improvedTest(self):
        fr = FileReader(FileReader.research_path + r"\Test\Files\accessions-and-sequences",
                        r"\accessions-and-sequences-part-18",
                        FileType.TSV)
        accessionsAndSequences = fr.from_file_to_dict(0, 1)

        print("accessions and sequences dictionary is achieved, and has " + str(len(accessionsAndSequences)) + " items")
        TestFunctions("from_file_to_dict", dictionary=accessionsAndSequences).print_first_lines_in_dict(1)

        self.bp.blastp_by_accessions("blastp", "nr", accessionsAndSequences, end_word="ORIGIN",
                                     translation_word="translation=")

    # def pairwise_alignment_inspector_test1(self):
    #     human_aa_seq = self.bp.get_aa_seq_by_human_gene_name("TCP1")
    #     print("human seq:", human_aa_seq)
    #     c_elegans_aa_seq = self.bp.get_aa_seq_by_c_elegans_gene_name("cct-1")
    #     print("C.elegans seq:", c_elegans_aa_seq)
    #     original_aa = "V"
    #     variant_index = 74
    #     print(self.bp.pairwise_alignment_inspector(human_seq=human_aa_seq,
    #                                                c_elegans_seq=c_elegans_aa_seq,
    #                                                original_amino_acid=original_aa,
    #                                                variant_index=variant_index))

    # testing if answer_count works correctly - not conserved
    # def pairwise_alignment_inspector_test2(self):
    #     human_aa_seq = self.bp.get_aa_seq_by_human_gene_name("TCP1")
    #     print("human seq:", human_aa_seq)
    #     c_elegans_aa_seq = self.bp.get_aa_seq_by_c_elegans_gene_name("cct-1")
    #     print("C.elegans seq:", c_elegans_aa_seq)
    #     original_aa = "T"
    #     variant_index = 53
    #     print(self.bp.pairwise_alignment_inspector(human_seq=human_aa_seq,
    #                                                c_elegans_seq=c_elegans_aa_seq,
    #                                                original_amino_acid=original_aa,
    #                                                variant_index=variant_index))
    #
    # # similar
    # def pairwise_alignment_inspector_test3(self):
    #     human_aa_seq = self.bp.get_aa_seq_by_human_gene_name("TCP1")
    #     print("human seq:", human_aa_seq)
    #     c_elegans_aa_seq = self.bp.get_aa_seq_by_c_elegans_gene_name("cct-1")
    #     print("C.elegans seq:", c_elegans_aa_seq)
    #     original_aa = "I"
    #     variant_index = 49
    #     print(self.bp.pairwise_alignment_inspector(human_seq=human_aa_seq,
    #                                                c_elegans_seq=c_elegans_aa_seq,
    #                                                original_amino_acid=original_aa,
    #                                                variant_index=variant_index))
    #
    # # doesn't exist in human sequence
    # def pairwise_alignment_inspector_test4(self):
    #     human_aa_seq = self.bp.get_aa_seq_by_human_gene_name("TCP1")
    #     print("human seq:", human_aa_seq)
    #     c_elegans_aa_seq = self.bp.get_aa_seq_by_c_elegans_gene_name("cct-1")
    #     print("C.elegans seq:", c_elegans_aa_seq)
    #     original_aa = "I"
    #     variant_index = 83
    #     print(self.bp.pairwise_alignment_inspector(human_seq=human_aa_seq,
    #                                                c_elegans_seq=c_elegans_aa_seq,
    #                                                original_amino_acid=original_aa,
    #                                                variant_index=variant_index))

### TESTING ###

# BioPythonTest().blastForAllTest()
# # BioPythonTest().blastpForAll_improvedTest()
# BioPythonTest().pairwise_alignment_inspector_test1()
# BioPythonTest().pairwise_alignment_inspector_test2()
# BioPythonTest().pairwise_alignment_inspector_test3()
# BioPythonTest().pairwise_alignment_inspector_test4()
