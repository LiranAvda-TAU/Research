import requests
from requests.exceptions import HTTPError

from Code.Files.FileReader import FileReader
from Code.Utils.Strings import Strings


class HttpRequester:

    def __init__(self, url=""):
        self.url = url
        self.id_upi_converter = FileReader(FileReader.research_path + r"\Data",
                                           r"\human_gene_id_upi.txt").from_file_to_dict_with_plural_values(0, 1, True)
        # self.id_upa_converter = FileReader(FileReader.research_path + r"\Data",
        #                                    r"\human_gene_id_upa.txt").from_file_to_dict_with_plural_values(0, 1, True)

    def make_request(self):
        # sending get request and saving the response as response object
        try:
            response = requests.get(url=self.url)
            # If the response was successful, no Exception will be raised
            response.raise_for_status()
        except HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
            return None
        except Exception as err:
            print(f'Other error occurred: {err}')
            return None
        else:
            # success
            return str(response.content)

    @staticmethod
    def get_uniprot_html(gene_id):
        try:
            r = requests.get("https://www.uniprot.org/uniprot/?query=" + gene_id + "&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score")
        except Exception as e:
            print("Exception in get_uniprot_html:", e)
            print("Communication failure, please check your internet connection")
            return None
        if not r.ok:
            r.raise_for_status()
            print("Something went wrong with get_uniprot_html while trying to extract reviewed ids, please check")
            return None
        return r.text

    def get_longest_human_protein_sequence_from_uniprot(self, gene_id):
        chosen_seq = ''
        if gene_id not in self.id_upi_converter:
            return None
        for upi in self.id_upi_converter[gene_id]:
            if upi == "":
                continue
            request_url = "https://www.ebi.ac.uk/proteins/api/uniparc/upi/" + upi + "?rfTaxId=9606"
            try:
                r = requests.get(request_url, headers={"Accept": "text/x-fasta"})
            except Exception as e:
                print("Exception in get_longest_human_protein_sequence_from_uniprot:", e)
                print("Communication failure, please check your internet connection")
                return None
            if not r.ok:
                r.raise_for_status()
                print("Something went wrong with get_longest_human_protein_sequence_from_uniprot() while trying to extract sequence"
                      ", please check")
                return None

            response_body = r.text
            optional_seq = Strings.from_fasta_seq_to_seq(response_body)
            if len(optional_seq) > len(chosen_seq):
                chosen_seq = optional_seq
        return chosen_seq

    def get_protein_sequence_from_ensembl(self, gene_id):
        if not gene_id:
            return None
        request_url = self.url + gene_id + "?type=protein;multiple_sequences=1"
        try:
            r = requests.get(request_url, headers={"Accept": "text/x-fasta"})
        except Exception as e:
            print("Communication failure, please check your internet connection:", e)
            return None
        if not r.ok:
            r.raise_for_status()
            print("Something went wrong with get_longest_human_protein_sequence_from_uniprot() while trying to extract sequence"
                  ", please check")
            return None
        return r.text

    @staticmethod
    def get_transcript(gene_id):
        chosen_transcript = ''
        transcript_ids = HttpRequester(url="http://rest.wormbase.org/rest/widget/gene/").\
            get_list_of_transcript_ids(gene_id)
        if not transcript_ids or not len(transcript_ids):
            return None
        for transcript_id in transcript_ids:
            transcript = HttpRequester(url="http://rest.wormbase.org/rest/field/transcript/").\
                get_transcript_by_transcript_id(transcript_id)
            if transcript and len(transcript) > len(chosen_transcript):
                chosen_transcript = transcript
                print(transcript_id, "is better:", len(transcript))
        if not chosen_transcript:
            print("Couldn't find transcript")
            return None
        return chosen_transcript

    def get_transcript_by_transcript_id(self, transcript_id):
        if not transcript_id:
            return None
        request_url = self.url + transcript_id + "/unspliced_sequence_context"
        try:
            r = requests.get(request_url)
        except Exception as e:
            print("Communication failure in get_transcript_by_gene_id, please check your internet connection:", e)
            return None
        if not r.ok:
            r.raise_for_status()
            print("Something went wrong with get_transcript_by_gene_id while trying to extract transcript"
                  ", please check")
            return None
        try:
            sequence = r.json()['unspliced_sequence_context']['data']['positive_strand']['sequence']
        except Exception as e:
            print("Couldn't extract sequence from worm base:", e)
            return None
        return sequence

    def get_list_of_transcript_ids(self, gene_id):
        transcript_ids = []
        # url = "http://rest.wormbase.org/rest/widget/gene/"
        if not gene_id:
            return None
        request_url = self.url + gene_id + "/sequences"
        try:
            r = requests.get(request_url)
        except Exception as e:
            print("Communication failure in get_list_of_transcript_ids, please check your internet connection:", e)
            return None
        if not r.ok:
            r.raise_for_status()
            print("Something went wrong with get_list_of_transcript_ids while trying to extract transcript ids"
                  ", please check")
            return None
        try:
            records = r.json()['fields']['gene_models']['data']['table'][0]['model']
            if type(records) == list:
                for record in records:
                    transcript_ids.append(record['label'])
            if type(records) == dict:  # one dict, no list of dicts
                transcript_ids.append(records['label'])
        except Exception as e:
            print("Error in get_list_of_transcript_ids: couldn't extract transcript ids:", e)
            return None
        print("list of transcript ids:", *transcript_ids, sep="\n")
        return transcript_ids

    @staticmethod
    def get_protein_seq_by_uniprot_swissprot_id(uniprot_swissprot_id):
        try:
            r = requests.get("https://www.uniprot.org/uniprot/" + uniprot_swissprot_id + ".fasta")
            response_body = r.text
        except Exception as e:
            print("Problem in function get_protein_seq_by_uniprot_swissprot_id while trying to extract sequence for",
                  uniprot_swissprot_id, ":", e)
            return None
        seq = Strings.from_fasta_seq_to_seq(response_body)
        return seq if seq else None


