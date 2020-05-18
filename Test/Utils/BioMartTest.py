from Code.Files.FileReader import FileReader
from Code.Utils.BioMart import BioMart

class BioMartTest:

    @staticmethod
    def test_get_uniport_swissport_id():
        human_gene_ids = FileReader.get_human_genes_ids()
        bm = BioMart()
        for gene_id in human_gene_ids:
            uni_id = bm.get_uniport_swissport_id(gene_id)
            print("for gene id:", gene_id, "uniport_swissport id:", uni_id)

### TESTING ###

# BioMartTest.test_get_uniport_swissport_id()