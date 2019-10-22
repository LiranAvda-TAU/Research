from Bio import Entrez
from Bio import pairwise2
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import Record
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

from Code.Enum.BlastType import BlastType
from Code.Enum.FileType import FileType
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.Ensembl import Ensembl
from Code.Utils.Strings import Strings

E_VALUE_THRESH = 0.04


class BioPython:
    def __init__(self):
        self.c_elegans_id_multiple_accessions = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                                           r"\\all-genes-ids-and-accession-numbers.txt",
                                                           FileType.TSV).fromFileToDictWithPluralValues(0, 1)

    @staticmethod
    def get_human_gene_sequence(geneName: str, converter: dict):
        if geneName in converter:
            geneId = converter[geneName]
        else:
            print("There is no record for this gene in the dictionary")
            exit()
        Entrez.email = "liranavda@gmail.com"
        handle = Entrez.esearch(db="gene", retmax=10, term=geneId)
        record = Entrez.read(handle)
        handle.close()
        recordId = record["IdList"][0]
        gene = Entrez.efetch(db="gene", id=recordId, rettype="gb", retmode="text")
        data = gene.read()
        return data

    # receives a list of genes' ids(WB) and returns a dict of a gene id(WB) and its id(number)
    @staticmethod
    def get_genes_id(c_elegans_genes: list):
        print("we have " + str(len(c_elegans_genes)) + " genes to extract id to. let's begin!\n")
        genes = {}
        index = 0
        for gene in c_elegans_genes:
            Entrez.email = "liranavda@gmail.com"
            gene_id = BioPython.get_gene_id(gene)
            if gene_id:
                genes[gene] = gene_id
                index += 1
                if not index % 100:
                    print("got to " + str(index) + "!\n")
            else:
                print("Gene id can't be found using Entrez: " + gene)
                continue
        return genes

    @staticmethod
    def get_gene_id(c_elegans_gene: str):
        Entrez.email = "liranavda@gmail.com"
        try:
            handle = Entrez.esearch(db="gene", retmax=1, term=c_elegans_gene)
            record = Entrez.read(handle)
            handle.close()
            if len(record["IdList"]) > 0:
                gene_id = record["IdList"][0]
                return gene_id
        except:
            return None

    @staticmethod
    def blast(program, db, accession_number):
        print("Accession Number: " + str(accession_number))
        try:
            result_handle = NCBIWWW.qblast(program=program, database=db, sequence=accession_number)
        except ValueError:
            result_handle = NCBIWWW.qblast(BlastType.fromNtoP(program), db, accession_number)
        record = NCBIXML.read(result_handle)
        result_handle.close()
        for alignment in record.alignments:
            if alignment.hit_def.startswith("Homo sapiens"):
                return True
        return False
        # print("****Alignment****")
        # print("sequence: ", alignment.title)
        # print("length: ", alignment.length)
        # for hsp in alignment.hsps:
        # if hsp.expect < E_VALUE_THRESH:
        #     print("****HSP****")
        #     print("e value:", hsp.expect)
        #     print(hsp.query[0:75] + "...")
        #     print(hsp.match[0:75] + "...")
        #     print(hsp.sbjct[0:75] + "...")

    @staticmethod
    def blast_with_seq(program, db, accession_number, sequence):
        state = False
        print("Accession Number: " + str(accession_number))
        result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence,
                                       entrez_query="Homo sapiens[organism]")
        record: Record = NCBIXML.read(result_handle)
        result_handle.close()
        for alignment in record.alignments:
            if BioPython.get_species(alignment.hit_def).startswith("Homo sapiens"):
                state = True
                print("****Alignment****")
                print(alignment.hit_id)
                print(alignment.hit_def)
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print("****HSP****")
                        print("e value:", hsp.expect)
                        print("hsp score:", hsp.score)
        return state

    # receives (1) an accession numbers and (2) the relevant sequencec, runs BLAST and returns a list of all
    # id of homo sapiens reported (first 50 by default)
    @staticmethod
    def pipelineBlastWithSeq(program, db, sequence):
        hit_ids = []
        result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence,
                                       entrez_query="Homo sapiens[organism]")
        record: Record = NCBIXML.read(result_handle)
        result_handle.close()
        for alignment in record.alignments:
            if BioPython.get_species(alignment.hit_def).startswith("Homo sapiens"):
                hit_id_record = alignment.hit_id
                hit_id = hit_id_record[hit_id_record.find("|") + 1:hit_id_record.find("|", hit_id_record.find("|") + 1)]
                extra_letter = hit_id_record[hit_id_record.find(hit_id) +
                                             len(hit_id) + 1:hit_id_record.find(hit_id) + len(hit_id) + 2]
                if extra_letter:
                    hit_id = hit_id + "_" + extra_letter
                hit_ids.append(hit_id)
        return hit_ids

    @staticmethod
    def blastForAll(blast_type, db, accessions_dict: dict):
        f = open("blast-for-all", "w+")
        index = 0
        keys = accessions_dict.keys()
        for key in keys:
            accession = accessions_dict[key]
            res = BioPython.blast(blast_type, db, accession)
            index += 1
            if not index % 100:
                print("got to " + str(index) + "!")
            f.write(key + "\t" + accession + "\t" + str(res))
            print(res)

    # receives a dictionary of gene ids and accession numbers
    @staticmethod
    def blastpForAll(blast_program, db, accessions_dict: dict, end_word, translation_word):
        index = 0
        gene_ids = accessions_dict.keys()
        for gene_id in gene_ids:
            accession = accessions_dict[gene_id]
            seq = BioPython.get_aa_seq_of_longest_isoform(accession_number=accession, end_word=end_word,
                                                          translation_word=translation_word)
            res = BioPython.blast_with_seq(blast_program, db, accession, seq)
            index += 1
            if not index % 100:
                print("got to " + str(index) + "!")
            print(res)

    # receives a dictionary of accession numbers and sequences
    @staticmethod
    def blastpForAll_improved(blast_program, db, accessions_dict: dict):
        accessions = accessions_dict.keys()
        for accession in accessions:
            seq = accessions_dict[accession]
            print("seq: " + seq)
            res = BioPython.blast_with_seq(blast_program, db, accession, seq)
            print(res)

    # receives (1) a sequence with white spaces, counts only the nucleotide letters, and returns the count
    @staticmethod
    def get_sequence_length(seq: str):
        return seq.count('t') + seq.count('a') + seq.count('c') + seq.count('g')

    @staticmethod
    def get_longest_isoform(accessions: set(), index_word, database: str = "nucleotide"):
        Entrez.email = "liranavda@gmail.com"
        longest_isoform = ""
        longest_isoform_length = 0
        if len(accessions) == 1:
            return accessions.pop()
        for accession_number in accessions:
            info = Entrez.efetch(db=database, id=accession_number, rettype="gb", retmode="text")
            record = info.read()
            info.close()
            if index_word in record:
                index = record.find(index_word)
                nt_seq = record[index + len(index_word):]
                seq_length = BioPython.get_sequence_length(nt_seq)
                if seq_length > longest_isoform_length:
                    longest_isoform = accession_number
                    longest_isoform_length = seq_length
            else:
                print(index_word + "does not exist in record for " + accession_number)
                exit()
        if not longest_isoform:  # longest isoform stayed the empty string
            print("Error: no isoform was found longer than the empty string")
            exit()
        return longest_isoform

    @staticmethod
    def get_aa_seq_of_longest_isoform(accession_number, end_word="ORIGIN", translation_word="translation=",
                                      database: str = "nucleotide"):
        Entrez.email = "liranavda@gmail.com"
        info = Entrez.efetch(db=database, id=accession_number, rettype="gb", retmode="text")
        record = info.read()
        info.close()
        # extracting aa sequence out of chosen isoform
        translation_index = record.find(translation_word)
        end_index = record.find(end_word)
        seq = record[translation_index + len(translation_word):end_index]
        extracted_seq = BioPython.get_amino_acid_seq(seq)
        return extracted_seq

    @staticmethod
    def get_amino_acid_seq(seq: str):
        new_seq = ""
        amino_acid_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                           'T', 'W', 'Y', 'V']
        for letter in seq:
            if letter in amino_acid_list:
                new_seq += letter
        return new_seq

    @staticmethod
    def get_species(alignment_hit: str):
        start_index = alignment_hit.find("[")
        end_index = alignment_hit.find("]")
        species = alignment_hit[start_index + 1:end_index]
        if species == "imported":
            return BioPython.get_species(alignment_hit[end_index + 1:])
        return species

    @staticmethod
    def make_accession_number_and_seq_dict(accessions, end_word, start_word):
        accessionsAnsSequences = {}
        for accession in accessions:
            seq = BioPython.get_aa_seq_of_longest_isoform(accession, end_word, start_word)
            accessionsAnsSequences[accession] = seq
        return accessionsAnsSequences

    @staticmethod
    def get_gene_name_from_protein_accession(hit_ids: [], search_term: str = "gene=\"", end_search_term: str = "\"\n"):
        print("Started working on extracting the gene name by the protein accession")
        accessions_and_names = {}
        Entrez.email = "liranavda@gmail.com"
        for hit_id in hit_ids:
            try:
                info = Entrez.efetch(db="protein", id=hit_id, rettype="gb", retmode="text")
                print("extraction was successful for: " + hit_id)
            except:
                print("extraction was not successful for: " + hit_id)
                continue
            record = info.read()
            info.close()
            # extracting aa sequence out of chosen isoform
            start_index = record.find(search_term)
            if start_index < 0:  # no gene name in record
                continue
            end_index = record.find(end_search_term, start_index)
            gene_name = record[start_index + len(search_term):end_index]
            accessions_and_names[hit_id] = gene_name
        return accessions_and_names

    @staticmethod
    def get_alignments(human_seq, c_elegans_seq):
        blosum_alignments = pairwise2.align.globaldx(human_seq, c_elegans_seq, matlist.blosum62)
        classic_alignments = pairwise2.align.globalxx(human_seq, c_elegans_seq)
        pam_alignments = pairwise2.align.globaldx(human_seq, c_elegans_seq, matlist.pam250)
        with_penalty_alignments = pairwise2.align.globalms(human_seq, c_elegans_seq, 2, -1, -.5, -.1)

        alignments = blosum_alignments + pam_alignments + classic_alignments + with_penalty_alignments
        return alignments

    # receives (1) the human sequence and (2) the C.elegans sequence, align them and check the maximum conservation out
    # of all the alignment, returns the score
    @staticmethod
    def get_conservation_score(human_seq, c_elegans_seq):
        alignments = BioPython.get_alignments(human_seq, c_elegans_seq)
        max_score = 0
        for alignment in alignments:
            score = BioPython.alignment_window_conservation_checker(alignment, 0, 0, whole_alignment=True)
            if score > max_score:
                max_score = score
        return max_score

    @staticmethod
    def pairwise_alignment_inspector(human_seq, c_elegans_seq, original_amino_acid, variant_index, window_size: int = 60):
        result = ''
        ultimate_result = ''
        count = 0
        c_elegans_aa_location = '-'
        window_conservation_score = 0
        alignments = BioPython.get_alignments(human_seq, c_elegans_seq)
        for alignment in alignments:
            human_alignment = alignment[0]
            c_elegans_alignment = alignment[1]
            alignment_index = Strings.getAminoAcidInLocationInAlignment(variant_index, human_alignment)
            status = human_alignment[alignment_index - 1], c_elegans_alignment[alignment_index - 1]

            if original_amino_acid == '-':  # termination
                return "Stop Codon", count, c_elegans_aa_location, window_conservation_score
            if human_alignment[alignment_index - 1] != original_amino_acid:
                result = "not in the human sequence, wanted: " + original_amino_acid + \
                         " and found: " + human_alignment[alignment_index - 1] + " alignment index " + str(alignment_index)
            elif human_alignment[alignment_index - 1] == c_elegans_alignment[
                        alignment_index - 1] == original_amino_acid:
                start = max(0, alignment_index - window_size)
                end = min(len(alignment[0]) - 1, alignment_index + window_size)
                window_conservation_score = BioPython.alignment_window_conservation_checker(alignment, start, end)
                count += 1
                ultimate_result = "conserved: " + original_amino_acid
                c_elegans_aa_location = BioPython.get_place_in_sequence(c_elegans_alignment, alignment_index)

            elif Strings.areAminoAcidsSimilar(human_alignment[alignment_index - 1],
                                              c_elegans_alignment[alignment_index - 1]):
                result = "not conserved, but similar: " + ", ".join(status)
            else:
                result = "not conserved: " + ", ".join(status)

        return (ultimate_result, count, c_elegans_aa_location, window_conservation_score) if ultimate_result else \
            (result, count, c_elegans_aa_location, window_conservation_score)

    @staticmethod
    def get_place_in_sequence(c_elegans_alignment, alignment_index):
        sequence_index = 0
        for char_index in range(alignment_index):
            if 'A' <= c_elegans_alignment[char_index] <= 'z':
                sequence_index += 1
        return sequence_index

    @staticmethod
    def alignment_window_conservation_checker(alignment, start, end, whole_alignment=False):

        lst = format_alignment(*alignment).split("\n")
        equals = 0
        if whole_alignment:
            # print(format_alignment(*alignment))
            segment_one = lst[0]
            segment_two = lst[2]
        else:
            segment_one = lst[0][start:end]
            segment_two = lst[2][start:end]
        for i in range(len(segment_one)):
            if segment_one[i] == segment_two[i]:
                equals += 1
        return equals / alignment[4]

    # receives (1) accession numbers for specific gene id, enact the function that chooses the longest isoform
    # of all accession number, and returns the chosen longest accession number
    @staticmethod
    def from_multiple_accessions_to_one_single_id(accessions):
        chosen_isoform = BioPython.get_longest_isoform(accessions=accessions, index_word="ORIGIN")
        return chosen_isoform

    @staticmethod
    def get_aa_seq_by_c_elegans_gene_name(c_elegans_gene_name):
        gene_id_WB = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
        return BioPython().get_aa_seq_by_c_elegans_gene_id_WB(gene_id_WB)

    def get_aa_seq_by_c_elegans_gene_id_WB(self, c_elegans_id_WB):
        c_elegans_id_number = Ensembl.get_ncbi_id_by_gene_id(c_elegans_id_WB)
        if not c_elegans_id_number:
            c_elegans_id_number = BioPython.get_gene_id(c_elegans_id_WB)
            if not c_elegans_id_number:
                print("Couldn't find C.elegans' id (number), moving on to the next gene")
                return None
        try:
            c_elegans_accession_numbers = set(self.c_elegans_id_multiple_accessions[c_elegans_id_number])
        except:
            print("cannot find accession numbers for gene", c_elegans_id_WB)
            return None
        try:
            c_elegans_accession_number = BioPython.from_multiple_accessions_to_one_single_id(
                c_elegans_accession_numbers)
        except:
            print("C.elegans gene's accession number cannot be extracted")
            return None
        try:
            c_elegans_seq = BioPython.get_aa_seq_of_longest_isoform(c_elegans_accession_number)
        except:
            print("C.elegans gene's sequence cannot be extracted")
            return None
        return c_elegans_seq

    @staticmethod
    def get_aa_seq_by_human_gene_name(human_gene_name):
        human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        if not human_gene_id:
            print("Human gene id for", human_gene_name, "cannot be found")
            return None
        human_seq = HttpRequester().get_human_protein_sequence_from_uniProt(human_gene_id)
        if not human_seq:
            print("human sequence cannot be extracted, perhaps upi cannot be found")
            return None
        return human_seq

# human_gene_name = "CAPZA1"
# c_elegans_gene_name = "cap-1"
# human_seq = BioPython.get_aa_seq_by_human_gene_name(human_gene_name)
# print("human seq:", human_seq)
# c_elegans_seq = BioPython().get_aa_seq_by_c_elegans_gene_name(c_elegans_gene_name)
# print("c-elegans seq", c_elegans_seq)
