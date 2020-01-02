import requests
from requests.exceptions import HTTPError

from Code.Files.FileReader import FileReader
from Code.Utils.Strings import Strings


class HttpRequester:

    def __init__(self, url=""):
        self.url = url
        self.id_upi_converter = FileReader(FileReader.research_path + r"\Data",
                                           r"\human_gene_id_upi.txt").fromFileToDictWithPluralValues(0, 1, True)

    def makeRequest(self):
        # sending get request and saving the response as response object
        try:
            response = requests.get(url=self.url)
            # If the response was successful, no Exception will be raised
            response.raise_for_status()
        except HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
            exit()
        except Exception as err:
            print(f'Other error occurred: {err}')
            exit()
        else:
            # success
            return str(response.content)

    def get_human_protein_sequence_from_uniProt(self, gene_id):
        chosen_seq = ''
        if gene_id not in self.id_upi_converter:
            return None
        for upi in self.id_upi_converter[gene_id]:
            if upi == "":
                continue
            requestURL = "https://www.ebi.ac.uk/proteins/api/uniparc/upi/" + upi + "?rfTaxId=9606"
            try:
                r = requests.get(requestURL, headers={"Accept": "text/x-fasta"})
            except:
                print("Communication failure, please check your internet connection")
                exit()
            if not r.ok:
                r.raise_for_status()
                print("Something went wrong with get_human_protein_sequence_from_uniProt() while trying to extract sequence"
                      ", please check")
                return None

            responseBody = r.text
            optional_seq = Strings.fromFastaSeqToSeq(responseBody)
            if len(optional_seq) > len(chosen_seq):
                chosen_seq = optional_seq
        return chosen_seq

    def get_protein_sequence_from_ensembl(self, gene_id):
        requestURL = self.url + gene_id + "?type=protein;multiple_sequences=1"
        try:
            r = requests.get(requestURL, headers={"Accept": "text/x-fasta"})
        except:
            print("Communication failure, please check your internet connection")
            exit()
        if not r.ok:
            r.raise_for_status()
            print("Something went wrong with get_human_protein_sequence_from_uniProt() while trying to extract sequence"
                  ", please check")
            return None
        return r.text