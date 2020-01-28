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

    # retrieves the gene WB id
    @staticmethod
    def get_c_elegans_gene_id_by_gene_name(gene_name):
        try:
            record = ensembl_rest.symbol_lookup('caenorhabditis elegans', gene_name)
        except:
            print("No valid lookup found for symbol", gene_name)
            return None
        return record['id']

    # receives a gene name and a species ("human" or "C.elegans") and returns the gene id
    @staticmethod
    def get_gene_id_by_gene_name(gene_name, species):
        if species == "Human":
            return Ensembl.get_human_gene_id_by_gene_name(gene_name)
        else:
            return Ensembl.get_c_elegans_gene_id_by_gene_name(gene_name)

    @staticmethod
    def get_c_elegans_genes_ids_by_genes_names(genes_names: list):
        genes_ids = []
        for gene_name in genes_names:
            gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(gene_name)
            if gene_id:
                genes_ids.append(gene_id)
        return genes_ids

    @staticmethod
    def get_genes_names_by_genes_ids(genes_ids):
        return [Ensembl.get_gene_name_by_gene_id(gene_id) for gene_id in genes_ids]

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

    # receives (1) list of genes names of human or c.elegans and (2) the string "human" or "c.elegans", and returns a
    # new list with the relevant genes ids.
    @staticmethod
    def convert_from_names_to_ids(genes_list, subject):
        new_list = []
        error = None
        for gene in genes_list:
            gene_id = None
            if subject == "Human":
                gene_id = Ensembl.get_human_gene_id_by_gene_name(gene)
            elif subject == "C.elegans":
                gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(gene)
            if gene_id:
                new_list.append(gene_id)
            else:
                error = "The name " + gene + " is not listed as a " + subject + " gene in our sources"

        return new_list, error

# print(Ensembl.get_nt_sequence_by_gene_name("cct-1"))