import requests
from pybiomart import Dataset

from Code.Http.HttpRequester import HttpRequester
from Code.Utils.Strings import Strings


class BioMart:
    def __init__(self):
        self.dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

    def get_uniport_swissport_id(self, gene_id):
        df = self.dataset.query(attributes=['ensembl_gene_id', 'uniprotswissprot']).dropna()
        uni_ids = df.loc[df['Gene stable ID'] == gene_id]['UniProtKB/Swiss-Prot ID']
        if uni_ids.size == 1:
            return uni_ids.iloc[0]
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
            chosen_uni_id = None
            chosen_seq = ''
            for uni_id in lst:
                try:
                    r = requests.get("https://www.uniprot.org/uniprot/" + uni_id + ".fasta")
                    response_body = r.text
                except Exception as e:
                    print("Problem in function get_uniport_swissport_id while trying to extract sequence for",
                          uni_id, "in gene", gene_id, ":", e)
                    return None
                optional_seq = Strings.from_fasta_seq_to_seq(response_body)
                if len(optional_seq) > len(chosen_seq):
                    chosen_seq = optional_seq
                    chosen_uni_id = uni_id
            return chosen_uni_id


