from pybiomart import Dataset
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.Ensembl import Ensembl


class BioMart:
    def __init__(self):
        self.dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        self.df = self.dataset.query(attributes=['ensembl_gene_id', 'uniprotswissprot']).dropna()

    def get_swissprot_sequence_from_biomart(self, gene_id):
        uni_ids = self.df.loc[self.df['Gene stable ID'] == gene_id]['UniProtKB/Swiss-Prot ID']
        if uni_ids.size == 1:
            return HttpRequester.get_protein_seq_by_uniprot_swissprot_id(uni_ids.iloc[0])
        elif uni_ids.size == 0:
            return None
        else:  # many ids
            lst = list(uni_ids)
            # check if uniprot reviewed all ids
            reviewed_html = HttpRequester.get_uniprot_html(gene_id)
            for uni_id in lst[:]:
                if uni_id not in reviewed_html:
                    lst.remove(uni_id)
            # choose by length
            chosen_seq = ''
            for uni_id in lst:
                optional_seq = HttpRequester.get_protein_seq_by_uniprot_swissprot_id(uni_id)
                if len(optional_seq) > len(chosen_seq):
                    chosen_seq = optional_seq
            return chosen_seq

            # access function

    def get_human_protein_from_uniprot_by_gene_id(self, gene_id):
        # first try canonical sequence
        canonical_protein_sequence = self.get_swissprot_sequence_from_biomart(gene_id)
        if canonical_protein_sequence:
            return canonical_protein_sequence
        else:
            # get longest sequence
            longest_seq = HttpRequester().get_longest_human_protein_sequence_from_uniprot(gene_id)
            if not longest_seq:
                return None
            return longest_seq

    def get_human_protein_from_uniprot_by_gene_name(self, human_gene_name):
        human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        if not human_gene_id:
            print("Human gene id for", human_gene_name, "cannot be found")
            return None
        else:
            return self.get_human_protein_from_uniprot_by_gene_id(human_gene_id)


