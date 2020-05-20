from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.BioMart import BioMart


class BioMartTest:

    @staticmethod
    def test_get_swissprot_sequence_from_biomart():
        human_gene_ids = FileReader.get_human_genes_ids()
        bm = BioMart()
        for gene_id in human_gene_ids:
            uni_seq = bm.get_swissprot_sequence_from_biomart(gene_id)
            print("for gene id:", gene_id, "uniport_swissport seq:", uni_seq)

    @staticmethod
    def test_uniprot_ids():
        # checks if there are uniprotswissprot ids that won't be recognized as those in the gene's uniprot page
        human_gene_ids = FileReader.get_human_genes_ids()
        bm = BioMart()
        reviewed = 0
        for gene_id in human_gene_ids:
            reviewed_html = HttpRequester.get_uniprot_html(gene_id)
            uni_ids = bm.df.loc[bm.df['Gene stable ID'] == gene_id]['UniProtKB/Swiss-Prot ID']
            if uni_ids.size == 1:
                uni_id = uni_ids.iloc[0]
                if uni_id not in reviewed_html:
                    print("id:", uni_id, "for gene:", gene_id, "doesn't seen as reviewed")
                else:
                    reviewed += 1
                    if not reviewed % 1000:
                        print(reviewed, "are good")
### TESTING ###
if False:
    BioMartTest.test_get_swissprot_sequence_from_biomart()

if False:
    BioMartTest.test_uniprot_ids()
