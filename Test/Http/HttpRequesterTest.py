from Code.Http.HttpRequester import HttpRequester
from Code.Utils.BioMart import BioMart


class HttpRequesterTest:

    @staticmethod
    def makeRequestTest():
        url = "http://ortholist.shaye-lab.org/"
        geneName = "usp11"
        httpRequester = HttpRequester(url, geneName)
        httpRequester.make_request()

    @staticmethod
    def test_get_transcript():
        print("###TEST: test_get_transcript###")
        list_of_gene_ids = 'WBGene00000018', 'WBGene00018050', 'WBGene00018064', 'WBGene00018072', 'WBGene00018119',\
                           'WBGene00200942', 'WBGene00200958', 'WBGene00200968', 'WBGene00201208', 'WBGene00198430',\
                           'WBGene00271642', 'WBGene00019039', 'WBGene00018933', 'WBGene00000377'
        for gene_id in list_of_gene_ids:
            print("for gene id:", gene_id)
            HttpRequester.get_transcript(gene_id)

    @staticmethod
    def test_canonical_vs_longest_protein():
        list_of_human_genes_ids = ['ENSG00000166206', 'ENSG00000198888', 'ENSG00000107854', 'ENSG00000040531',
                                   'ENSG00000196826', 'ENSG00000253729', 'ENSG00000131374', 'ENSG00000285013',
                                   'ENSG00000177994', 'ENSG00000197275', 'ENSG00000157483', 'ENSG00000226490']
        for gene_id in list_of_human_genes_ids:
            print("canonical vs longest for:", gene_id)
            canonical = BioMart().get_swissprot_sequence_from_biomart(gene_id)
            print("canonical:", canonical)
            longest = HttpRequester().get_longest_human_protein_sequence_from_uniprot(gene_id)
            print("longest:", longest)
            if canonical != longest:
                print("different!")
            else:
                print("equals!")

# TESTING #

# HttpRequesterTest.makeRequestTest()
# HttpRequesterTest.test_get_transcript()
HttpRequesterTest.test_canonical_vs_longest_protein()
