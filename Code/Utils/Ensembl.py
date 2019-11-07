import ensembl_rest


class Ensembl:

    @staticmethod
    def get_human_gene_id_by_gene_name(gene_name):
        try:
            record = ensembl_rest.symbol_lookup('homo sapiens', gene_name)
        except:
            print("No valid lookup found for symbol", gene_name)
            return None
        return record['id']

    @staticmethod
    def get_c_elegans_gene_id_by_gene_name(gene_name):
        try:
            record = ensembl_rest.symbol_lookup('caenorhabditis elegans', gene_name)
        except:
            print("No valid lookup found for symbol", gene_name)
            return None
        return record['id']

    @staticmethod
    def get_c_elegans_genes_ids_by_genes_names(genes_names: list):
        genes_ids = []
        for gene_name in genes_names:
            gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(gene_name)
            genes_ids.append(gene_id)
        return genes_ids

    @staticmethod
    def get_gene_name_by_gene_id(gene_id):
        try:
            record = ensembl_rest.lookup(gene_id)
        except:
            print("No valid lookup found for id", gene_id)
            return None
        return record['display_name']

    @staticmethod
    def get_nt_sequence_by_gene_id(gene_id):
        try:
            record = ensembl_rest.sequence_id(gene_id)
        except:
            print("Id", gene_id, "not found")
            return None
        return record['seq']

    @staticmethod
    def get_nt_sequence_by_gene_name(gene_name):
        try:
            return Ensembl.get_nt_sequence_by_gene_id(Ensembl.get_c_elegans_gene_id_by_gene_name(gene_name))
        except:
            return None

    @staticmethod
    def get_ncbi_id_by_gene_id(gene_id):
        try:
            record = ensembl_rest.xref_id(gene_id)
        except:
            print("Id:", gene_id, "not found")
            return None
        for record_cell in record:
            if "WBGene" not in record_cell["primary_id"] and "ENSG" not in record_cell["primary_id"]:
                return record_cell["primary_id"]
        return None

# print(Ensembl.get_nt_sequence_by_gene_name("cct-1"))