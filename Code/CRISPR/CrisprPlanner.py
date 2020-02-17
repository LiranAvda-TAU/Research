import random
import re
from statistics import mean
import time

from Code.CRISPR.NamedTuples.CodonData import CodonData
from Code.CRISPR.NamedTuples.CodonMutation import CodonMutation
from Code.CRISPR.NamedTuples.PointMutation import PointMutation
from Code.CRISPR.NamedTuples.ReattachmentSection import ReattachmentSection
from Code.CRISPR.NamedTuples.RestrictionMutation import RestrictionMutation
from Code.CRISPR.NamedTuples.RestrictionSite import RestrictionSite
from Code.CRISPR.NamedTuples.SequenceSites import SequenceSites
from Code.CRISPR.Enum.AminoAcid import AminoAcid
from Code.CRISPR.Enum.DNASection import DNASection
from Code.CRISPR.Enum.MutationDirection import MutationDirection
from Code.CRISPR.Enum.RestrictionSiteType import RestrictionSiteType
from Code.Files.FileReader import FileReader
from Code.Utils.BioPython import BioPython


class CrisprPlanner:
    possibilities = {'A': 'C', 'C': 'T', 'T': 'G', 'G': 'A'}

    codon_dic = {'A': {'GCT', 'GCC', 'GCA', 'GCG'}, 'C': {'TGT', 'TGC'}, 'D': {'GAT', 'GAC'}, 'E': {'GAA', 'GAG'},
                 'F': {'TTT', 'TTC'}, 'G': {'GGT', 'GGC', 'GGA', 'GGG'}, 'H': {'CAT', 'CAC'},
                 'I': {'ATT', 'ATC', 'ATA'}, 'K': {'AAA', 'AAG'}, 'L': {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'},
                 'M': {'ATG'}, 'N': {'AAT', 'AAC'}, 'P': {'CCT', 'CCC', 'CCA', 'CCG'}, 'Q': {'CAA', 'CAG'},
                 'R': {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'}, 'S': {'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'},
                 'T': {'ACT', 'ACC', 'ACA', 'ACG'}, 'V': {'GTT', 'GTC', 'GTA', 'GTG'}, 'W': {'TGG'},
                 'Y': {'TAT', 'TAC'}, 'STOP': {'TAA', 'TAG', 'TGA'}}

    amino_acid_dic = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TGT': 'C', 'TGC': 'C', 'GAT': 'D', 'GAC': 'D',
                      'GAA': 'E', 'GAG': 'E', 'TTT': 'F', 'TTC': 'F', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                      'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K', 'TTA': 'L',
                      'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATG': 'M', 'AAT': 'N', 'AAC': 'N',
                      'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
                      'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                      'AGT': 'S', 'AGC': 'S', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GTT': 'V', 'GTC': 'V',
                      'GTA': 'V', 'GTG': 'V', 'TGG': 'W', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
                      'TGA': 'STOP'}

    variables = {'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
                 'H': ('A', 'C', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                 'N': ('A', 'C', 'G', 'T'), 'R': ('A', 'G'), 'S': ('C', 'G'),
                 'V': ('A', 'C', 'G'), 'W': ('A', 'T'), 'Y': ('C', 'T')}

    codon_usage = {'UUU': 23.3, 'UCU': 16.7, 'UAU': 17.5, 'UGU': 11.2, 'UUC': 23.9, 'UCC': 10.6, 'UAC': 13.7,
                   'UGC': 9.1,
                   'UUA': 9.8, 'UCA': 20.6, 'UAA': 1.6, 'UGA': 1.4, 'UUG': 20.0, 'UCG': 12.2, 'UAG': 0.6, 'UGG': 11.1,
                   'CUU': 21.2, 'CCU': 8.8, 'CAU': 14.1, 'CGU': 11.2, 'CUC': 14.8, 'CCC': 4.4, 'CAC': 9.2, 'CGC': 5.1,
                   'CUA': 7.9, 'CCA': 26.1, 'CAA': 27.4, 'CGA': 12.1, 'CUG': 12.1, 'CCG': 9.7, 'CAG': 14.4, 'CGG': 4.7,
                   'AUU': 32.2, 'ACU': 18.9, 'AAU': 30.2, 'AGU': 12.1, 'AUC': 18.9, 'ACC': 10.4, 'AAC': 18.3,
                   'AGC': 8.4,
                   'AUA': 9.5, 'ACA': 20.0, 'AAA': 37.5, 'AGA': 15.4, 'AUG': 26.1, 'ACG': 8.9, 'AAG': 25.8, 'AGG': 4.0,
                   'GUU': 24.1, 'GCU': 22.4, 'GAU': 35.8, 'GGU': 10.9, 'GUC': 13.6, 'GCC': 12.6, 'GAC': 17.1,
                   'GGC': 6.7,
                   'GUA': 9.8, 'GCA': 19.8, 'GAA': 40.8, 'GGA': 31.7, 'GUG': 14.3, 'GCG': 8.2, 'GAG': 24.5, 'GGG': 4.4}

    restriction_sites = ['ACCGGT', 'GTGCAC', 'ATTAAT', 'GGATCC', 'AGATCT', 'ATCGAT', 'GAATTC', 'GATATC', 'AGCGCT',
                         'GGCGCT', 'AGCGCC', 'GGCGCC', 'AAGCTT', 'GGTGA', 'GGTACC', 'TGGCCA', 'CCATGG', 'CATATG',
                         'GCTAGC', 'GCGGCCGC', 'CTGCAG', 'CAGCTG', 'GAGCTC', 'CCTGCAGG', 'TACGTA', 'ACTAGT',
                         'TCTAGA', 'CTCGAG', 'CCCGGG', 'AATT', 'GTAC', 'TCGA', 'GATC', 'CCGG']

    def __init__(self, gene_name, aa_mutation_site, sense_strand, amino_acid_sequence: str = ""):
        self.gene_name = gene_name
        self.sense_strand = sense_strand
        self.anti_sense_strand = CrisprPlanner.get_anti_sense_strand(sense_strand)
        self.amino_acid_sequence = amino_acid_sequence if amino_acid_sequence else \
            BioPython().get_aa_seq_by_c_elegans_gene_name(gene_name)
        self.amino_acid_mutation_site = aa_mutation_site
        self.sense_mutation_site = -1
        self.anti_sense_mutation_site = -1
        self.mutated_strand = None
        self.mutation_direction = None
        self.ssODN_direction = None
        self.mutated_sites = []
        self.restriction_enzymes = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Code\CRISPR",
                                              r"\parsed_restriction_enzymes.txt").get_restriction_enzymes_list()

    def get_relevant_strand(self, direction):
        if direction > 0:
            return self.sense_strand
        else:
            return self.anti_sense_strand

    def plan_my_crispr(self,
                       from_aa: AminoAcid,
                       to_aa: AminoAcid,
                       check_consistency: bool = False,
                       window_size: int = 30,
                       PAM_size: int = 3):
        print("restriction enzymes:", self.restriction_enzymes)
        print("Gene's sense sequence:", self.sense_strand)
        if not self.amino_acid_sequence:
            # amino acid couldn't be extracted
            self.amino_acid_sequence = input("Insert amino acid sequence for gene " + self.gene_name + ":")
        print("Gene's amino acid sequence: ", self.amino_acid_sequence)
        self.sense_mutation_site = self.find_site_in_nt_seq(amino_acid_site=self.amino_acid_mutation_site,
                                                            check_sequence_consistency=check_consistency)
        print("mutation site in sense strand is: " + str(self.sense_mutation_site), ", and the codon is:",
              self.sense_strand[self.sense_mutation_site:self.sense_mutation_site + 3])
        self.anti_sense_mutation_site = self.get_mutation_site_for_anti_sense(self.sense_mutation_site)
        print("mutation site in anti-sense strand: " + str(self.anti_sense_mutation_site))

        # get optional sequences for crRNA
        sense_strand_sequences = CrisprPlanner.get_list_of_potential_sequences(self.sense_strand,
                                                                               self.sense_mutation_site,
                                                                               window_size,
                                                                               PAM_size)
        anti_sense_strand_sequences = CrisprPlanner. \
            get_list_of_potential_sequences(self.anti_sense_strand,
                                            self.anti_sense_mutation_site,
                                            window_size,
                                            PAM_size)
        print("Sense crRNAs: ")
        for sequence in sense_strand_sequences:
            print(sequence)
            pam_site_start = sequence[1][1]
            pam_site_end = sequence[1][1] + 3
            print("PAM sequence: " + self.sense_strand[pam_site_start:pam_site_end])
        print("Anti Sense crRNAs: ")
        for sequence in anti_sense_strand_sequences:
            print(sequence)
            pam_site_start = sequence[1][1]
            pam_site_end = sequence[1][1] + 3
            print("PAM sequence: " + self.anti_sense_strand[pam_site_start:pam_site_end])

        chosen_crrna, strand_direction = CrisprPlanner.choose_best_crrna(sense_strand_sequences,
                                                                         anti_sense_strand_sequences)
        # now we have our cr_rna
        print("chosen crRNA: " + str(chosen_crrna))
        strand_str = "anti-" if strand_direction < 0 else ""
        print("The relevant strand is " + strand_str + "sense")

        pam_site_start = chosen_crrna[1][1]
        pam_site_end = chosen_crrna[1][1] + 2
        pam_sites = (pam_site_start, pam_site_end)

        self.ssODN_direction, self.mutation_direction = CrisprPlanner. \
            choose_ssODN_strand(self.sense_mutation_site if strand_direction > 0 else self.anti_sense_mutation_site,
                                strand_direction,
                                pam_site_start)
        print("Mutation direction is:", str(self.mutation_direction), "thus ssODN direction is:",
              str(self.ssODN_direction))
        if self.ssODN_direction != strand_direction:
            print("The ssODN is on the opposite strand!")
            # ssODN is not on the strand from which we got the crRNA, thus we need to find its pam sites again
            pam_sites = (self.find_index_in_parallel_strand(pam_site_end),
                         self.find_index_in_parallel_strand(pam_site_start))
            print("pam sites on the ssODN strand:", pam_sites)

        # to be mutated
        self.mutated_strand = self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand
        ssODN_codon_mutation_site = self.sense_mutation_site if self.ssODN_direction > 0 else self.anti_sense_mutation_site

        # getting the zone where all mutations will be located: after to before DSB
        # sanity check
        print("mutation codon is at", ssODN_codon_mutation_site)

        first_original_strand = self.mutated_strand
        first_original_mutated_sites = self.mutated_sites
        # first try - when at last we try to insert restriction site
        insertion_result = self.insert_mutations(ssODN_codon_mutation_site, from_aa, to_aa, pam_sites, RestrictionSiteType.INSERTED)
        if not insertion_result:
            self.mutated_strand = first_original_strand
            self.mutated_sites = first_original_mutated_sites
            deletion_result = self.insert_mutations(ssODN_codon_mutation_site, from_aa, to_aa, pam_sites, RestrictionSiteType.REMOVED)
            if not deletion_result:
                print("Could not insert or remove restriction site")
                return False
        return True

    # this function is responsible for inserting all needed mutations to the string: to change amino acid, to prevent
    # reattachment and to add\remove restriction site
    def insert_mutations(self, ssODN_codon_mutation_site, from_aa, to_aa, pam_sites, restriction_site_type):
        # 1. mutation to change codon to change the amino acid
        possible_codon_mutations = self.get_possible_mutations_demands(from_aa, to_aa, self.sense_mutation_site)
        second_original_strand = self.mutated_strand
        second_original_mutated_sites = self.mutated_sites
        for possible_codon_mutation in possible_codon_mutations:
            print("Now working on codon mutation:", possible_codon_mutation)
            self.mutated_strand = second_original_strand
            self.mutated_sites = second_original_mutated_sites
            self.apply_codon_mutation(ssODN_codon_mutation_site,
                                      possible_codon_mutation.dict_of_mutations)

            # making sure that the strand is good so far
            print("mutated strand after adding codon mutation:", self.mutated_strand)
            print("sites changed:", self.mutated_sites)

            # 2. mutations to change nt to prevent re-attachments - in PAM or crRNA sequence
            # if mutation is DOWNSTREAM the section to be changed is PAM site, if UPSTREAM - crRNA.
            if self.mutation_direction == MutationDirection.DOWNSTREAM:
                section_to_mutate = ReattachmentSection(1, DNASection.PAM_SITE, SequenceSites(pam_sites[0], pam_sites[1]))
            else:
                section_to_mutate = ReattachmentSection(4, DNASection.CR_RNA,
                                                        SequenceSites(pam_sites[0] - 20, pam_sites[0] - 6))

            print("section to mutate to prevent re-attachment:", section_to_mutate.section_sites, "in",
                  section_to_mutate.section_type)

            possible_mutations, number_of_mutants = self.get_possible_codon_mutations(section_to_mutate,
                                                                                      ssODN_codon_mutation_site)
            # subsets with indexes that represent mutations that sum up to enough mutations overall
            valid_index_subsets = CrisprPlanner.get_valid_subsets(number_of_mutants, possible_mutations)
            third_original_strand = self.mutated_strand
            third_original_mutated_sites = self.mutated_sites
            print("all valid index subsets:", valid_index_subsets)
            for index_subset in valid_index_subsets:
                self.mutated_strand = third_original_strand
                self.mutated_sites = third_original_mutated_sites
                point_mutations = self.get_point_mutations(possible_mutations, index_subset) # also mutates the sequence
                print("now working on index subset:", index_subset, "and point mutations:", point_mutations)
                print("mutated strand after adding anti-reattachment mutations:", self.mutated_strand)
                self.mutated_sites.extend(point_mutations)
                print("all mutated sites:", self.mutated_sites)

                # 3. mutations to add/remove restriction sites
                print("Add or Remove restriction sites:")
                result = self.add_remove_restriction_sites(pam_sites,
                                                           ssODN_codon_mutation_site,
                                                           restriction_site_type)
                if result:
                    return True
        return False

    # receives (1) number of mutations needed, (2) the list of all possible anti-reattachment section mutations, and
    # (3) power set of all mutations indexes, and leave out only the ones with sufficient amount of mutations. if there,
    # leaves only subsets with maximum number of mutations
    @staticmethod
    def get_valid_subsets(number_of_mutations, possible_mutations):
        index_power_set = CrisprPlanner.get_power_set(list(range(len(possible_mutations))))
        valid_subsets = []
        underscored_subsets = []
        max_underscored = 0
        for subset in index_power_set:
            mutations_in_subset = 0
            sum_usage = 0
            for index in subset:
                if mutations_in_subset >= number_of_mutations:
                    break
                mutations_in_subset += possible_mutations[index][1].number_of_mutations
                sum_usage += possible_mutations[index][1].usage
            else:
                if mutations_in_subset < number_of_mutations:  # not enough mutations
                    if mutations_in_subset >= max_underscored:
                        if mutations_in_subset > max_underscored:
                            underscored_subsets = []
                            max_underscored = mutations_in_subset
                        underscored_subsets.append((subset, mutations_in_subset, sum_usage))
                else:
                    valid_subsets.append((subset, mutations_in_subset, sum_usage))
        if valid_subsets:
            # sort subsets according to number of mutstions max usage
            print("valid subsets:", valid_subsets)
            valid_subsets.sort(key=lambda mutation: (mutation[1], -mutation[2]))
            print("valid subsets after sort:", valid_subsets)
            return list(map(lambda item: item[0], valid_subsets))
        else:
            print("WARNING: the maximum amount of mutations is", max_underscored, ", there are no subsets with",
                  number_of_mutations, "mutations")
            # sort according to number of mutations first, and codons usage second
            underscored_subsets.sort(key=lambda mutation: (mutation[1], -mutation[2]))
            return list(map(lambda item: item[0], underscored_subsets))





    # receives (1) the section that needs to be mutated to prevent reattachment, (2) the indexes of the amino acid's
    # mutation codon on the ssODN strand,(3) the mutation direction, and (4) the sites that have
    #  changed while mutating the amino acid, computes which indexes can be changed in order to insert the section with
    # silent mutations and returns the possibilities for silent mutations
    def get_possible_codon_mutations(self,
                                     section_to_mutate: ReattachmentSection,
                                     ssODN_mutation_codon):
        print("***GET MUTATIONS TO PREVENT REATTACHMENT***")
        # format of codons_data: ('GCG', (1, 3))
        codons_data = CrisprPlanner.get_relevant_codons(section_to_mutate, ssODN_mutation_codon, self.mutated_strand)
        possible_mutations = []
        for codon_data in codons_data:
            # all codons_data that translate to the same amino acid as the relevant codon does
            same_aa_codons = list(CrisprPlanner.get_similar_codons(codon_data.codon))
            print("for codon:", codon_data, "same aa codons:", same_aa_codons)
            if section_to_mutate.section_type == DNASection.PAM_SITE:
                # removes PAM mutations of sequence NGA and NGG
                CrisprPlanner.check_PAM_mutation(codon_data, same_aa_codons, section_to_mutate, self.mutation_direction)
            # removes codons that run over the nucleotide change to change amino acid
            CrisprPlanner.check_already_mutated_sites(codon_data, same_aa_codons, self.mutated_sites)

            optional_codon_mutations = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(codon_data.codon, same_aa_codons)
            for optional_mutation in optional_codon_mutations:
                possible_mutations.append((codon_data, optional_mutation))
        # sort mutations according to which is closest the the mean of former point mutations
        possible_mutations.sort(key=lambda x: abs(x[0].codon_sites.start + mean(x[1].dict_of_mutations.keys()) -
                                                  mean(list(map(lambda item: item.index, self.mutated_sites)))))
        print("possible mutations:", possible_mutations)
        number_of_mutants = CrisprPlanner.get_number_of_mutants(section_to_mutate, self.mutated_sites)
        return possible_mutations, number_of_mutants

    # receives (1) all possible mutations to prevent reattachment and (2) list of indexes corresponding to the mutations
    # list and returns the point mutations in format of namedTuples: PointMutations(index, new nucleotide)
    def get_point_mutations(self, possible_mutations, indexes_subset):
        # possible mutation format: (CodonData(codon='TTT', codon_sites=SequenceSites(start=1174, end=1176)), OptionalCodonChange('TTC', 1, {2: 'C'}, 23.9))
        point_mutations = []
        for index in indexes_subset:
            possible_mutation = possible_mutations[index]
            dict_of_changes = possible_mutation[1].dict_of_mutations
            for key in dict_of_changes:
                point_mutations.append(PointMutation(possible_mutation[0].codon_sites.start + key, dict_of_changes[key]))
                self.mutated_strand = CrisprPlanner.change_char_in_string(self.mutated_strand,
                                                                          possible_mutation[0].codon_sites.start + key,
                                                                          dict_of_changes[key])
        return point_mutations

    @staticmethod
    def get_power_set(lst, power_set=[[]]):
        if lst:
            new_power_set = []
            new_item = lst[0]
            for item in power_set:
                new_power_set.append(item)
                new_power_set.append(item + [new_item])
            return CrisprPlanner.get_power_set(lst[1:], new_power_set)
        else:
            return power_set

    # Takes into consideration whether mutations to change amino acid have also modified the reattachment sequence, and
    # thus we need fewer mutant to introduce to the reattachment sequence
    @staticmethod
    def get_number_of_mutants(section_to_mutate, mutated_sites):
        number_of_mutants = section_to_mutate.number_of_mutations
        extra = 0 if section_to_mutate.section_type == DNASection.PAM_SITE else 2
        for mutated_site in mutated_sites:
            if section_to_mutate.section_sites.start <= mutated_site.index <= section_to_mutate.section_sites[1]+extra:
                number_of_mutants -= 1
                print("already a mutant in nucleotide:", mutated_site)
        if number_of_mutants < 0:
            number_of_mutants = 0
        return number_of_mutants

    # receives (1) codon_data the contains the codon and its sites, (2) possible codons that the former codons could be
    # mutated into but keep the amino acid, and (3) sites that have changed and the new mutations shouldn't run over.
    # The function goes through each of the mutated site (input 3) and checks if it is inside the range of the new
    # codon sites. If so, it checks whether in the possible codons list there are codons that run over that site (change
    #  nucleotide), and if there are, removes them from the possible codons list.
    @staticmethod
    def check_already_mutated_sites(codon_data, possible_codons, mutated_sites):
        print("checking if some codons need to be removed to prevent them from changing mutations")
        for mutated_site in mutated_sites:
            if codon_data.codon_sites.start <= mutated_site.index <= codon_data.codon_sites.end:
                for possible_codon in possible_codons[:]:
                    # mutated_codon_details format = (1, {2: 'G'}, 14.3)
                    mutated_codon_details = CrisprPlanner.how_to_get_b_codons_from_a(codon_data.codon, [possible_codon])
                    for key in mutated_codon_details.dict_of_mutations:
                        # for each nucleotide mutation, key is the nucleotide place in the codon: 0,1,2
                        if codon_data.codon_sites[0] + key == mutated_site.index:
                            print(mutated_codon_details, codon_data.codon_sites, "removed because", mutated_site.index,
                                  "has already been mutated")
                            possible_codons.remove(possible_codon)
                            break

    # receives (1) list of mutations to check, (2) list of possible codons in that site, (3) section to mutate (PAM site
    #  or crRNA and the details), and (4) mutation direction (upstream or downstream) and checks if one of the PAM
    # mutations is substituting NGG site with an NGA or NGG site, and if so, removes it from the valid list of codons
    # returned.
    @staticmethod
    def check_PAM_mutation(codon_data: CodonData, possible_codons, section_to_mutate, mutation_direction):

        if codon_data.codon_sites.start < section_to_mutate.section_sites.start:
            # before beginning of PAM site
            if codon_data.codon_sites.end == section_to_mutate.section_sites.end:
                # only contains first nucleotide
                if mutation_direction == MutationDirection.DOWNSTREAM:
                    # cannot be T to prevent NGA
                    for possible_codon in possible_codons[:]:
                        if possible_codon[0] == 'T':
                            print("possible codon removed:", possible_codon, "in", codon_data)
                            possible_codons.remove(possible_codon)
                else:
                    # mutation only changes the N part of NGG PAM site - not useful
                    possible_codons.clear()
            elif codon_data.codon_sites.end > section_to_mutate.section_sites.start:
                # contains first two amino acids
                if mutation_direction == MutationDirection.DOWNSTREAM:
                    for possible_codon in possible_codons[:]:
                        if possible_codon[1:] == 'TC':
                            print("possible codon removed:", possible_codon, "in", codon_data)
                            possible_codons.remove(possible_codon)
        elif codon_data.codon_sites == section_to_mutate.section_sites:
            for possible_codon in possible_codons[:]:
                if mutation_direction == MutationDirection.UPSTREAM:
                    # the site is NGG
                    if possible_codon[1:] == 'GA':
                        print("possible codon removed:", possible_codon, "in", codon_data)
                        possible_codons.remove(possible_codon)
                else:
                    # the site is CCN
                    if possible_codon[:2] == 'TC':
                        print("possible codon removed:", possible_codon, "in", codon_data)
                        possible_codons.remove(possible_codon)
        else:
            if codon_data.codon_sites.start < section_to_mutate.section_sites.end:
                # contains two last amino acids
                if mutation_direction == MutationDirection.UPSTREAM:
                    for possible_codon in possible_codons[:]:
                        if possible_codon[:2] == 'GA':
                            print("possible codon removed:", possible_codon, "in", codon_data)
                            possible_codons.remove(possible_codon)
            elif codon_data.codon_sites.start == section_to_mutate.section_sites.end:
                # contains only last PAM nucleotide
                if mutation_direction == MutationDirection.UPSTREAM:
                    for possible_codon in possible_codons[:]:
                        if possible_codon[0] == 'A':
                            print("possible codon removed:", possible_codon, "in", codon_data)
                            possible_codons.remove(possible_codon)
                else:
                    print("all mutations of this sort must be erased since mutation in that codon will leave out the NGG")
                    possible_codons.clear()

    # receives (1) a codon and returns all other codons that translated to the same amino acid
    @staticmethod
    def get_similar_codons(codon):
        amino_acid = CrisprPlanner.amino_acid_dic[codon]
        relevant_codons = CrisprPlanner.codon_dic[amino_acid].copy()
        relevant_codons.remove(codon)
        return relevant_codons

    # receives (1) the section to be mutated (PAM site or crRNA indexes), (2) the aa mutation codon indexes of the ssODN
    # strand, and (3) the mutated strand (sense or anti-sense) and returns a list of all codons and their indexes that
    # are part of the section to be mutated indexes
    # format of the results: ('GCG', (1, 3))
    @staticmethod
    def get_relevant_codons(section_to_mutate: ReattachmentSection, ssODN_mutation_codon, mutated_strand):
        codons_data = []
        remainder = abs(section_to_mutate.section_sites.start - ssODN_mutation_codon) % 3
        if section_to_mutate.section_sites.start >= ssODN_mutation_codon:
            start = section_to_mutate.section_sites.start - remainder
        else:
            if remainder > 0:
                start = section_to_mutate.section_sites.start - (3 - remainder)
            else:
                start = section_to_mutate.section_sites.start
        end = section_to_mutate.section_sites.end
        # if section is crRNA, we don't want to make changes in the DSB
        while start <= end:
            codon = mutated_strand[start:start + 3]
            codons_data.append(CodonData(codon, SequenceSites(start, start + 2)))
            start += 3
        print(str(len(codons_data)), "relevant codons that we can mutate:", codons_data)
        return codons_data

    # This function receives (1) an amino acid sequence site and (2) value that indicates whether to check the
    # correlation between amino acid sequence and nucleotide sequence (codon to amino acid) or not, and returns the
    # site of the mutation in the nucleotide sequence, meaning the site or the corresponding codon
    def find_site_in_nt_seq(self, amino_acid_site, check_sequence_consistency: bool = True):
        sense_without_utr = self.remove_utr()
        print("sense without utr:", sense_without_utr)
        if check_sequence_consistency:
            # consistency check
            nt_index = 0
            aa_index = 0
            while aa_index < len(self.amino_acid_sequence):
                CrisprPlanner.check_sequence_consistency(sense_without_utr[nt_index:nt_index + 3],
                                                         self.amino_acid_sequence[aa_index])
                nt_index += 3
                aa_index += 1

        nt_site_in_sense_without_utr = (amino_acid_site - 1) * 3
        utrs = 0
        seq_index = 0
        exon_index = 0
        while seq_index < len(self.sense_strand):
            if self.sense_strand[seq_index].islower():
                # region is intron or utr
                utrs += 1
            else:
                # region is exon
                exon_index += 1
                if exon_index == nt_site_in_sense_without_utr:
                    return exon_index + utrs
            seq_index += 1

    # receives (1) a codon string and (2) the relevant amino and checks if the codon translates to this amino acid
    @staticmethod
    def check_sequence_consistency(input_codon, amino_acid):
        codon_table = {'GCT': 'A', 'TGT': 'C', 'GAT': 'D', 'GAA': 'E', 'TTT': 'F', 'GGT': 'G', 'CAT': 'H',
                       'ATT': 'I', 'AAA': 'K', 'TTA': 'L', 'ATG': 'M', 'AAT': 'N', 'CCT': 'P', 'CAA': 'Q',
                       'CGT': 'R', 'TCT': 'S', 'ACT': 'T', 'GTT': 'V', 'TGG': 'W', 'TAT': 'Y', 'TTC': 'F',
                       'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATC': 'I', 'ATA': 'I',
                       'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCC': 'P',
                       'CCA': 'P', 'CCG': 'P', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCC': 'A', 'GCA': 'A',
                       'GCG': 'A', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'CAC': 'H', 'CAG': 'Q', 'AAC': 'N',
                       'AAG': 'K', 'GAC': 'D', 'GAG': 'E', 'TGC': 'C', 'TGA': 'STOP', 'CGC': 'R', 'CGA': 'R',
                       'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGC': 'G', 'GGA': 'G',
                       'GGG': 'G'}
        if codon_table[input_codon] != amino_acid:
            print("Instead of " + codon_table[input_codon], ", the amino acid in the sequence is " + amino_acid)

    # receives a (1) strand of nt, (2) the mutation site, (3) the window size around the mutation site where DSB point
    # can be located and (4) size of PAM sequence, checks in the strand for 20 bp sequences around the mutation size
    # that end in PAM sequence ('NGG') and returns a list of those sequences, as a tuple of the sequences and its start
    # and end indexes.
    @staticmethod
    def get_list_of_potential_sequences(strand, mutation_site, window_size, PAM_size, cr_rna_size: int = 20):
        potential_cr_rnas = []
        start_point = max(0, mutation_site - window_size + PAM_size)
        end_point = min(mutation_site + window_size + PAM_size, len(strand) - 1)
        print("start point: " + str(start_point) + ", end point: " + str(end_point))
        for i in range(start_point, end_point):
            if strand[i + 1:i + 3] == "GG":
                try:
                    potential_cr_rna = strand[i - cr_rna_size:i]
                    potential_cr_rnas.append((potential_cr_rna, (i - cr_rna_size, i)))
                except:
                    print("Not enough space for a whole crRNA, i = " + str(i))
                    continue
        return potential_cr_rnas

    # receives (1) the mutation position in the sense strand, and returns the position of this mutation in the
    # anti-sense
    def get_mutation_site_for_anti_sense(self, sense_mutation_site):
        return len(self.anti_sense_strand) - sense_mutation_site - 3

    # receives (1) the sense strand and returns its complementary, from 5' to 3'
    @staticmethod
    def get_anti_sense_strand(given_sense_strand):
        anti_sense_strand = ''
        pair_nucleotides = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
        for letter in given_sense_strand:
            anti_sense_strand += pair_nucleotides[letter]
        return anti_sense_strand[::-1]

    # receives (1) a list of sense strands and returns their complementary strands, from 5' to 3'
    @staticmethod
    def get_anti_sense_strands(given_sense_strands):
        complementary_strands = []
        for sense_strand in given_sense_strands:
            anti_sense_strand = CrisprPlanner.get_anti_sense_strand(sense_strand)
            complementary_strands.append(anti_sense_strand)
        return complementary_strands

    # this function gets a nucleotide sequence and adds to a dictionary all restriction sites in the sequence
    @staticmethod
    def find_restriction_sites(nucleotide_seq, start_point, methylated: bool = False):
        restriction_sites_in_seq = []
        for restriction_site in CrisprPlanner.restriction_sites:
            indexes = [m.start() for m in re.finditer('(?=' + restriction_site + ')', nucleotide_seq)]
            for index in indexes:
                if restriction_site == "GGTGA":
                    if index + 8 >= len(nucleotide_seq):
                        continue
                elif restriction_site == "GATC":
                    if not methylated:
                        continue
                restriction_sites_in_seq.append(RestrictionSite(SequenceSites(
                    start_point + index,
                    start_point + index + len(restriction_site)), restriction_site))
        print("restriction sites:", restriction_sites_in_seq)
        return restriction_sites_in_seq

    @staticmethod
    def change_char_in_string(seq, position, new_char):
        lst = list(seq)
        lst[position] = new_char
        return "".join(lst)

    @staticmethod
    def choose_best_crrna(sense_options: list, anti_sense_options: list):
        num = random.uniform(0, 1)
        if num > 0.5:
            return random.choice(sense_options), 1
        else:
            return random.choice(anti_sense_options), -1

    @staticmethod
    def choose_ssODN_strand(mutation_site, strand_direction, pam_site_start):
        DSB = pam_site_start - 3
        # default assignment
        mutation_direction = MutationDirection.UPSTREAM
        print("mutation site: " + str(mutation_site) + ", pam start: " + str(pam_site_start) + ", DSB: " + str(DSB))
        if mutation_site > DSB:
            mutation_direction = MutationDirection.DOWNSTREAM
        elif mutation_site < DSB:
            mutation_direction = MutationDirection.UPSTREAM
        if mutation_direction == MutationDirection.DOWNSTREAM:
            return -1 * strand_direction, mutation_direction
        else:
            return strand_direction, mutation_direction

    def find_index_in_parallel_strand(self, index):
        return len(self.sense_strand) - index - 1

    # this function will return the minimal nucleotides you need to change (and the range of options for change) in
    # order to get from codon a to amino acid b
    @staticmethod
    def how_to_get_b_from_a(a_codon, b_amino_acid: AminoAcid):
        b_codons = CrisprPlanner.codon_dic[b_amino_acid.value]
        return CrisprPlanner.how_to_get_b_codons_from_a(a_codon, b_codons)

    # this function will return the minimal nucleotides you need to change (and the range of options for change) in
    # order to get from codon a to possible codons b. result in format: ('GCG', (1, {0: 'G'}, 8.2))
    @staticmethod
    def how_to_get_b_codons_from_a(a_codon, b_codons):
        min_distance = 4
        options = {}
        for b_codon in b_codons:
            distance, demands = CrisprPlanner.codon_distance(a_codon, b_codon)
            if distance < min_distance:
                options = {}
                min_distance = distance
                options[b_codon] = (distance, demands, CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(b_codon)])
            elif distance == min_distance:
                options[b_codon] = (distance, demands, CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(b_codon)])
        if options:
            # result is a tuple of codon, details
            return CrisprPlanner.get_codon_with_higher_usage(options)
        return None

    # this function will return a list of all possible codon changes, sorted by the number of point mutations needed and
    # by the usage value of the codon, in order to get from codon a to possible codons b.
    # result in format of list of OptionlCodonChange
    @staticmethod
    def get_list_of_how_to_get_b_codons_from_a(a_codon, b_codons):
        options = []
        for b_codon in b_codons:
            distance, demands = CrisprPlanner.codon_distance(a_codon, b_codon)
            options.append(CodonMutation(b_codon,
                                         distance,
                                         demands,
                                         CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(b_codon)]))
        if options:
            # result is a list of sorted OptionalCodonChange
            options.sort(key=lambda option: (option.number_of_mutations, -option.usage))
        return options

    @staticmethod
    def codon_distance(a_codon, b_codon):
        distance = 0
        demands = {}
        for i in range(3):
            if a_codon[i] != b_codon[i]:
                distance += 1
                demands[i] = b_codon[i]
        return distance, demands

    @staticmethod
    def change_t_to_u(codon):
        new_codon = ""
        for nt in codon:
            if nt == 'T':
                new_codon += 'U'
            else:
                new_codon += nt
        return new_codon

    # gets a dictionary of options, such that the keys are tuples of the original codon (doesn't change) and new, and
    # the values are 3-tuples which consists of the number of mutations, the instructions (dictionary of codon indexes and
    # wanted nucleotides, and the usage percentage, and returns a MutatedCodonDetails that holds the said info + the key
    # codon with the highest codon usage
    @staticmethod
    def get_codon_with_higher_usage(codon_options):
        max_codon = ""
        max_usage = 0
        for codon in codon_options:
            value = codon_options[codon]
            usage = value[2]
            if usage > max_usage:
                max_usage = usage
                max_codon = codon
        return CodonMutation(max_codon, codon_options[max_codon][0], codon_options[max_codon][1], max_usage)

    # calculated the demands to go through the codon from the mutation site to one of the codons translating to the
    # wanted amino acid
    def get_possible_mutations_demands(self, from_amino_acid: AminoAcid, to_amino_acid: AminoAcid, mutation_site):
        sense_codon = self.sense_strand[mutation_site:mutation_site + 3]
        codon_options = self.codon_dic[from_amino_acid.value]
        if sense_codon not in codon_options:
            print(f'Codon in mutation site {mutation_site}: {sense_codon} does not translate into {from_amino_acid}')
            exit()
        a_codon = sense_codon if self.ssODN_direction > 0 else CrisprPlanner.get_anti_sense_strand(sense_codon)
        b_codons = CrisprPlanner.codon_dic[to_amino_acid.value]
        if self.ssODN_direction > 0:
            b_codons_info = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(a_codon, list(b_codons))
        else:
            b_codons_info = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(a_codon,
                                                                                 CrisprPlanner.get_anti_sense_strands(b_codons))
        print(len(b_codons_info), "possible codons change:", b_codons_info)
        print("The best option is to change", a_codon, "to", b_codons_info[0].codon, "with demands:", str(b_codons_info[0].dict_of_mutations))

        return b_codons_info

    # receives (1) the mutation site in the strand (internal indexing), and (2) the demands to make the change with
    # indexes in respect to the codon, not the strand (so, 0, 1, 2...) and returns the
    # mutated strand and a list of the changes in format of PointMutation tuples: (index, new nucleotide)
    def apply_codon_mutation(self, mutation_site, demands):
        for index in demands:
            nt = demands[index]
            position = mutation_site + int(index)
            self.mutated_sites.append(PointMutation(position, nt))
            self.mutated_strand = CrisprPlanner.change_char_in_string(self.mutated_strand, position, nt)

    @staticmethod
    def check_if_pam_site_mutated(mutated_sites: list, pam_sites: tuple):
        pam_site_mutated = False
        for mutated_site in mutated_sites:
            if pam_sites[0] <= mutated_site <= pam_sites[1]:
                pam_site_mutated = True
                break

        if pam_site_mutated:
            print("pam site mutated")
        else:
            print("pam site not mutated")
        return pam_site_mutated

    # receives (1) the mutation direction (UPSTREAM or DOWNSTREAM) and (2) the pam sites and according to the scheme
    # that shows that for UPSRTREAM mutations, the homology arm that can undergo changes is from the DSB and for
    # DOWNSTREAM mutations, the homology arm that can undergo changes is up until the DSB, returns the right zone of
    # 13 nt. The range is open in the end (doesn't include the end nucleotide)
    @staticmethod
    def get_mutations_zone(mutation_direction, pam_sites):

        if mutation_direction == MutationDirection.UPSTREAM:
            # untouched homology arm is until DSB
            DSB = pam_sites[0] - 3
        else:
            # untouched homology arm is after DSB
            DSB = pam_sites[1] + 4
        return SequenceSites(max(0, DSB - 13), DSB)

    @staticmethod
    def is_within_mutation_zone(mutation_zone, place):
        if mutation_zone[0] < place < mutation_zone[1]:
            return True
        else:
            return False

    # extracts the part of original and mutated sequences and checks to see whether one of them has a restriction site
    # the other one doesn't. returns the distinctive restriction sites if the lists don't match, and a None value
    # if they do
    def get_distinctive_restriction_sites(self, mutated_sites, mutated_strand):
        # first we'll find the boundaries in which to search
        min_site = min(list(map(lambda item: item.index, mutated_sites)))
        max_site = max(list(map(lambda item: item.index, mutated_sites)))
        max_restriction_site = max(list(map(lambda item: len(item), self.restriction_sites)))

        original_sequence = self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand
        original_section = original_sequence[min_site - max_restriction_site:max_site + max_restriction_site]
        mutated_section = mutated_strand[min_site - max_restriction_site:max_site + max_restriction_site]

        restrictions_sites_in_original_seq = self.find_restriction_sites(original_section,
                                                                         min_site - max_restriction_site)
        restrictions_sites_in_mutated_seq = self.find_restriction_sites(mutated_section,
                                                                        min_site - max_restriction_site)

        if self.do_restriction_lists_dictionaries_match(restrictions_sites_in_original_seq,
                                                        restrictions_sites_in_mutated_seq):
            return None
        else:
            return self.find_extra_restriction_sites(restrictions_sites_in_original_seq, restrictions_sites_in_mutated_seq)

    def get_restriction_sites(self, mutated_strand, mutation_zone):
        max_restriction_site = max(list(map(lambda item: len(item), self.restriction_sites)))

        # minus 3 to avoid changing codon after the DSB
        mutated_section = mutated_strand[mutation_zone.start - max_restriction_site:mutation_zone.end - 3 + max_restriction_site]
        restrictions_sites_in_mutated_seq = self.find_restriction_sites(mutated_section,
                                                                        mutation_zone.start - max_restriction_site)
        return restrictions_sites_in_mutated_seq

    def add_remove_restriction_sites(self,
                                     pam_sites,
                                     ssODN_mutation_codon,
                                     restriction_site_type):

        mutation_zone = self.get_mutations_zone(self.mutation_direction, pam_sites)
        if restriction_site_type == RestrictionSiteType.INSERTED:
            # first check if the mutated strand so far (with codon mutations) has different number of restriction sites
            # than the original one
            restriction_sites = self.get_distinctive_restriction_sites(self.mutated_sites, self.mutated_strand)
            if restriction_sites:  # value is not None
                print("dictionaries don't match!")  # different restriction sites
                # TODO
                return True

            else:  # dictionaries match
                print("let's modify some sites! first, trying to add...")
                print("mutation zone:", mutation_zone)
                restriction_mutation = self.get_new_restriction_site(mutation_zone,
                                                                     ssODN_mutation_codon,
                                                                     max_mutations=2)
                if restriction_mutation:
                    print("restriction mutation added!", restriction_mutation)
                    self.mutated_strand = restriction_mutation.mutated_strand
                    return True
                else:
                    return False

        else:  # trying to remove a restriction site
                print("now let's try to remove...")
                if self.choose_restriction_site_to_remove(mutation_zone, ssODN_mutation_codon):
                    return True
                else:
                    return False

    # receives the mutation zone and extracts all restriction sites there. after filtering and sorting, chooses a
    # restriction site to mutate and does so if possible
    def choose_restriction_site_to_remove(self, mutation_zone, ssODN_mutation_codon):
        rest_sites = self.get_restriction_sites(self.mutated_strand, mutation_zone)
        if rest_sites:
            print("No restriction sits in the relevant section")
        print("restriction sites:", rest_sites)
        CrisprPlanner.filter_out_restriction_sites_with_no_space(rest_sites, self.mutated_strand)
        if rest_sites:
            # if there are still restriction sites left, we  will now sort them
            self.sort_restriction_sites(rest_sites)
            print("rest sites after sorting:", rest_sites)
            chosen_index = CrisprPlanner.get_chosen_index(rest_sites)
            restriction_site_removed = False
            sites_tested = 0
            chosen_rest_site = None
            while not restriction_site_removed and sites_tested < len(rest_sites):
                chosen_rest_site = rest_sites[chosen_index]
                restriction_site_removed = self.remove_restriction_site(rest_sites[chosen_rest_site], ssODN_mutation_codon)
                chosen_index = (chosen_index + 1) % len(rest_sites)
                sites_tested += 1
            if restriction_site_removed:  # succeeded in finding and removing a restriction site
                print("Restriction site removed! chosen restriction site:", chosen_rest_site)
                return True
            else:
                print("No restriction site to remove was found")
                return False
        else:
            print("No restriction sites left after filtering out")
            return False

    # receives (1)the chosen restriction site and (2) the start of the ssODN mutation codon tries to insert silent mutations to remove it
    def remove_restriction_site(self, rest_site: RestrictionSite, ssODN_mutation_codon):
        # it's not a reattachment section, but I need this namedTuple for the get_relevant_codons method
        section_to_mutate = ReattachmentSection(1, DNASection.MutationZone, rest_site.index)
        possible_mutations, number_of_mutants = self.get_possible_codon_mutations(self, section_to_mutate,
                                                                                  ssODN_mutation_codon)
        fourth_original_strand = self.mutated_strand
        fourth_original_sites = self.mutated_sites
        for index in range(len(possible_mutations)):
            self.mutated_strand = fourth_original_strand
            self.mutated_sites = fourth_original_sites
            point_mutations = self.get_point_mutations(possible_mutations, [index])  # also mutates the sequence
            self.mutated_sites.append(point_mutations)
            print("Restriction site removed with:", point_mutations)
            return True
        return False

    @staticmethod
    def get_chosen_index(rest_sites: list):
        chosen_index = input("Dear user, we believe the first restriction site is the best candidate, but please"
                             "tell us what the index of your choise is (starting from 0). Otherwise, we will go"
                             "with the first one on the list above")
        try:
            chosen_index = int(chosen_index)
        except ValueError:
            chosen_index = 0
        else:
            if chosen_index < 0 or chosen_index > len(rest_sites) - 1:
                chosen_index = 0
        return chosen_index

    # receives a list of restriction sites (RestrictionSite) and sorts them according to what restriction sites the lab
    # has, its distance from same type restriction sites and frequency along the sequence
    def sort_restriction_sites(self, rest_sites):
        # first, sort according to existing restriction sites in the lab, second, sort according to the distance to next
        # same type restriction site, and third
        existing_restriction_sites = input("Of all the restriction sites enlishted above: \n" + self.restriction_sites
                                            + "\nwhich are the ones that are possession of your lab? write down a string"
                                              "of indexes, such as: 0 5 7 12")
        existing_restriction_sites_list = existing_restriction_sites.split(" ")
        rest_sites.sort(key=lambda r_site: (self.restriction_sites.index(r_site.site) in existing_restriction_sites_list,
                                            self.check_for_distance(r_site, vicinity=250),
                                            self.check_rareness(r_site, 250)),
                        reverse=True)

    def check_rareness(self, restriction_site, distance):
        site_length = len(restriction_site.site)
        mutated_section = self.mutated_strand[max(0, restriction_site.index.start-distance-site_length):
                                              min(len(self.mutated_strand)-1, restriction_site.index.end+distance+site_length)]
        occurrences = [m.start() for m in re.finditer('(?=' + restriction_site.site + ')', mutated_section)]
        return 1-len(occurrences)/len(mutated_section)

    # receives a list of restriction sites and the strand that they are taken from, and looks for other same restriction
    # sites in the vicinity of certain nucleotides. if the vicinity is too close, remove the restriction site from the
    # list
    @staticmethod
    def filter_out_restriction_sites_with_no_space(rest_sites, mutated_strand, vicinity: int = 50):
        for restriction_site in rest_sites[:]:
            if CrisprPlanner.check_for_distance(restriction_site, mutated_strand, vicinity) < vicinity:
                rest_sites.remove(restriction_site)

    # receives a restriction site and a strand, and checks in what distance there is another restriction site from the
    # same type
    def check_for_distance(self, restriction_site, vicinity: int = 50):
        distance = vicinity + 1
        site_length = len(restriction_site.site)
        mutated_section = self.mutated_strand[max(0, restriction_site.index.start-vicinity-site_length):
                                              min(len(self.mutated_strand)-1,
                                                  restriction_site.index.end+vicinity+site_length)]
        occurrences = [m.start() for m in re.finditer('(?=' + restriction_site.site + ')', mutated_section)]
        for occurrence in occurrences:
            if 0 < abs(occurrence-restriction_site.index.start) < distance:
                distance = occurrence
        return distance

    @staticmethod
    def subtract_two_lists(lst1, lst2):
        return list(set(lst1)-set(lst2))

    def get_mean_of_subtracted_lists(self, item: RestrictionMutation, shorter_list):
        subtracted_list = self.subtract_two_lists(item.mutated_sites, shorter_list)
        if not subtracted_list:
            # if empty
            subtracted_list = shorter_list
        return mean(list(map(lambda site: site.index, subtracted_list)))

    def get_new_restriction_site(self, mutation_zone, ssODN_mutation_codon,
                                 max_mutations: int = 2):
        possible_results = self. get_possible_restriction_sites(self.mutated_strand, mutation_zone, self.mutated_sites,
                                                                ssODN_mutation_codon, max_mutations)
        print("All possible restriction sites mutations:")
        for possible_result in possible_results:
            print(possible_result)
        mean_former_mutations_index = mean(list(map(lambda site: site.index, self.mutated_sites)))
        if possible_results:
            possible_results.sort(key=lambda item: (item.number_of_mutations,
                                                    abs(mean_former_mutations_index-self.get_mean_of_subtracted_lists(item, self.mutated_sites))))
            return possible_results[0]
        else:
            return None

    # this function extracts all possible codons in the mutation zone section and tries to insert silent mutation in
    # them, recursively, and collects all options that add a new restriction site to the strand. number of max_mutations
    # is chosen
    def get_possible_restriction_sites(self, mutated_strand, mutation_zone, mutated_sites, ssODN_mutation_codon,
                                       max_mutations, mutations_so_far: int = 0):
        possible_results = []
        section_to_mutate = ReattachmentSection(-1,
                                                DNASection.MUTATION_ZONE,
                                                SequenceSites(mutation_zone.start, mutation_zone.end - 3))
        codons_data = CrisprPlanner.get_relevant_codons(section_to_mutate, ssODN_mutation_codon, mutated_strand)
        for codon_data in codons_data:
            print("for codon:", codon_data.codon, "at", codon_data.codon_sites)
            # all codons_data that translate to the same amino acid as the relevant codon does
            same_aa_codons = list(CrisprPlanner.get_similar_codons(codon_data.codon))
            print("same codons:", same_aa_codons)
            CrisprPlanner.check_already_mutated_sites(codon_data, same_aa_codons, mutated_sites)
            for same_aa_codon in same_aa_codons:
                num_of_mutations = mutations_so_far
                mutated_codon_details = CrisprPlanner.how_to_get_b_codons_from_a(codon_data.codon, [same_aa_codon])
                print("adding mutation:", mutated_codon_details, codon_data.codon_sites)
                dict_of_changes = mutated_codon_details.dict_of_mutations
                current_mutated_strand = mutated_strand
                temporary_point_mutations = []
                for key in dict_of_changes:
                    current_mutated_strand = CrisprPlanner.change_char_in_string(current_mutated_strand,
                                                                                 codon_data.codon_sites.start + key,
                                                                                 dict_of_changes[key])
                    num_of_mutations += 1
                    temporary_point_mutations.append(PointMutation(codon_data.codon_sites.start + key, dict_of_changes[key]))
                restriction_sites = self.get_distinctive_restriction_sites(mutated_sites + temporary_point_mutations,
                                                                           current_mutated_strand)
                print("mutations so far:", num_of_mutations)
                if restriction_sites:
                    print(restriction_sites, codon_data, mutated_codon_details)
                    possible_results.append(RestrictionMutation(current_mutated_strand,
                                                                restriction_sites,
                                                                num_of_mutations,
                                                                mutated_sites + temporary_point_mutations))
                if num_of_mutations < max_mutations:
                    print("going recursive!", num_of_mutations, max_mutations)
                    possible_results.extend(self.get_possible_restriction_sites(current_mutated_strand, mutation_zone,
                                                                                 mutated_sites + temporary_point_mutations,
                                                                                ssODN_mutation_codon,
                                                                                max_mutations, num_of_mutations))
        return possible_results

    # gets two dictionaries of the format of: {(3, 9): 'TCTAGA'}) and checks if both dictionaries have the same restriction sites
    @staticmethod
    def do_restriction_lists_dictionaries_match(sites_in_original_dic, sites_in_mutated_dic):
        if len(sites_in_original_dic) != len(sites_in_mutated_dic):
            return False
        else:
            for restriction_site in sites_in_original_dic:
                filtered_list = list(filter(lambda item: item.index == restriction_site.index, sites_in_mutated_dic))
                # if list is empty, no restriction sites in that index
                if not filtered_list:
                    return False
                else:
                    if restriction_site.site != filtered_list[0].site:
                        return False
        return True

    # receives the original sequence (with only codon mutations) and the optional sequence with more mutations that has
    # more\less restrictions sites than the original, and returns the added restriction site and indexes, and whether
    # the change is added restriction site or removed. format: (2, {'TCGA': (19, 23), 'CCGG': (0, 4)})
    @staticmethod
    def find_extra_restriction_sites(original_seq_restriction_sites, mutant_seq_restriction_sites):
        extra_restrictions_sites = {}

        # restriction site was removed
        for restriction_site in original_seq_restriction_sites:
            filtered_list = list(filter(lambda item: item.index == restriction_site.index, mutant_seq_restriction_sites))
            if not filtered_list or \
                            restriction_site.site != filtered_list[0].site:
                extra_restrictions_sites[restriction_site] = RestrictionSiteType.REMOVED

        # restriction site was added
        for restriction_site in mutant_seq_restriction_sites:
            filtered_list = list(
                filter(lambda item: item.index == restriction_site.index, original_seq_restriction_sites))
            if not filtered_list or \
                            restriction_site.site != filtered_list[0].site:
                extra_restrictions_sites[restriction_site] = RestrictionSiteType.INSERTED

        print(extra_restrictions_sites)
        return extra_restrictions_sites

    def remove_utr(self):
        seq_without_utr = ''
        for letter in self.sense_strand:
            if 'A' <= letter <= 'Z':
                seq_without_utr += letter
            else:
                continue
        return seq_without_utr

    def find_sites(self, seq):
        print("seq:", seq)
        for i in range(len(seq)):
            for j in range(i + 1, len(seq) + 1):
                for enzyme in self.restriction_enzymes:
                    for derivative in enzyme.derivatives:
                        if seq[i:j] == derivative.replace("/", ""):
                            print(enzyme.site, i, j)


work1 = False
if work1:
    repo_1_sense_strand = "ATGGACTTTCAGAACAGAGCTGGAGGAAAAACGGGAAGCGGAGGAGTGGCTTCGGCCGCCGATGCTGGTGTTGATCGACGGGAACGGCTCCGCCAGTTGGCTCTAGAGACAATTGATCTTCAAAAGGATCCGTATTTCATGCGAAATCACATTGGAACGTACGAATGCAAGCTGTGTCTTACTCTTCACAACAATGAAGGATCTTATTTGGCACATACACAAGGAAAGAAGCATCAAGCGAATCTTGCACGGCGTGCCGCTAAAGAACAATCTGAACAACCATTTCTACCAGCTCCACAGAAAGCTGCAGTTGAAACTAAAAAGTTTGTGAAAATCGGACGTCCTGGATACAAGGTAACAAAAGAACGTGATCCAGGAGCTGGCCAGCAAGCACTTCTCTTCCAAATTGATTATCCGGAGATTGCTGACGGTATTGCGCCACGTCATCGATTTATGTCTGCTTATGAGCAAAAGATTCAGCCTCCAGACAAGAGATGGCAATACCTCTTGTTTGCTGCTGAGCCGTATGAAACGATTGGATTCAAAATTCCATCAAGgtgaggctttacaacattttagcacttttctatctcatagttacgattaaaaaaattgtatataccaagtaattttttccagAGAAGTTGACAAATCTGAAAAATTTTGGACGATGTGGAACAAAGACACGAAGCAATTCTTCTTACAAGTCGCATTCAAATTGGAACGACTCGATGATCAGCCGTACTATTGAtactctatgtttttatctttttgatttcaaaattcaaaacaattttttcgtgtttttcgatgatctaacaataaattattttcctttttttt"
    # sense_strand = "gcattgtaaggagaagccgggtaattaatacgataggcgccgttacaaaccgccaactggtgatcattattctctgaaaATGGAGCCGCGGACAGACGGAGCAGAATGCGGTGTCCAGgtattaattttccccgcctagattttccaatttcatattgttttcagGTATTTTGTCGTATTCGGCCGCTCAATAAGACCGAGGAGAAGAATGCGGACCGTTTCCTGCCCAAATTCCCTTCCGAGGACAGTATATCGCTTGGGgtgagtaatacaaaggggtcatagggaacaattatgtcaacagggacgggaagcacgggggatgcaggtgtgtcaattctctcacatgacacattcatctgtttgaaaagtacacgaaaagtgcaaagttgaatatatatatatatatcgattgatttgttggaatttttcagGGAAAAGTATACGTGTTCGATAAAGTGTTCAAGCCGAACACCACGCAAGAGCAAGTGTACAAAGGAGCCGCTTATCACATCGTACAGGATGTATTATCCGGTTATAATGGAACAGTTTTTGCATATGGACAAACATCTTCCGGAAAAACACATACAATGGAGgtaggaattatgaaaaccttgataattacgtagaatgcgacaaagacaatcaagttgtaatatcaacagtgcaaatctttactgattaatgaaaagaaaagtttgagaactaattttcagcagttatttccgaaatcgaatgcccgaaagatttttgataatttttacgtttaaaatattcggcgctcgatgaaattaacaatataatttaattattttcatattttttacagGGAGTAATCGGTGATAATGGCTTGTCGGGAATCATTCCACGTATCGTTGCTGACATCTTTAACCACATTTATAGTATGGACGAGAATCTTCAATTTCACATCAAAGTGTCCTATTATGAAATTTACAACGAGAAGATTCGAGATTTATTAGACCCCGAGAAGGTCAATTTGTCCATTCATGAAGATAAAAATCGAGTGCCATACGTGAAGGGAGCCACCGAACGGTTTGTTGGAGGACCCGATGAGGTTCTTCAGGCAATCGAAGATGGAAAATCCAACAGAATGGTTGCAGTTACGAgtgagtaaacttaaaattaaacaaattaacatgtgaacgaaatttcagACATGAACGAACATTCTTCTCGATCTCATTCCGTCTTCTTGATTACTGTGAAACAAGAACATCAGACAACAAAGAAACAGCTCACCGGAAAGCTTTATCTTGTTGATTTGGCTGGTTCTGAGAAAGTGAGCAAAACTGGAGCTCAAGGAACAGTTTTAGAAGAAGCCAAAAACATCAACAAGTCACTTACTGCACTCGGAATAGTTATTTCAGCATTGGCTGAAGGAACTgtgagttgtttaaattatgaccttcttaaaacgaatatttatttcagAAATCTCATGTTCCATATCGTGATTCCAAACTGACTCGTATTCTTCAAGAATCTCTAGGAGGAAATTCCCGTACTACAGTTATTATTTGTGCTTCTCCGTCACATTTCAACGAAGCTGAAACTAAATCCACACTTTTGTTCGGAGCACGTGCGAAGACTATCAAGAATGTTGTACAAATCAACGAAGAGCTCACAGCAGAAGAATGGAAACGGCGATATGAGAAAGAAAAAGAGAAGAATACTCGATTGGCCGCCCTTCTCCAGGCAGCGGCTTTGGAACTTTCACGCTGGCGTGCTGGAGAATCAGTGTCTGAGGTTGAATGGGTCAATCTATCAGATTCTGCTCAAATGGCTGTGTCGGAAGTTTCTGGTGGGTCGACTCCACTCATGGAACGTTCGATTGCTCCAGCTCCTCCAATGCTAACTTCTACAACTGGCCCGATCACTGACGAAGAGAAGAAGAAGTACGAAGAGGAACGTGTCAAACTGTATCAGCAACTCGACGAGAAAGATGATGAGATTCAAAAAGTTTCGCAAGAGCTTGAGAAGCTTAGACAACAAGTTCTTCTCCAAGAAGAAGCTTTGGGAACTATGCGTGAAAACGAGGAGCTGATCCGTGAAGAGAACAACCGATTCCAAAAAGAAGCTGAAGACAAGCAGCAAGAAGGAAAGGAAATGATGACAGCTCTGGAAGAGATTGCTGTCAACTTGGATGTTCGACAAGCAGAATGCGAAAAATTGAAGAGAGAGTTGGAAGTTGTTCAAGAAGATAACCAGAGTTTGGAAGATCGAATGAACCAAGCAACATCACTCCTCAATGCTCATCTTGACGAATGTGGTCCAAAAATCCGTCATTTCAAAGAAGGAATCTACAATGTTATTCGTGAATTCAACATTGCTGACATTGCCTCTCAAAATGATCAACTTCCTGATCACGATCTTCTGAACCATGTCAGAATCGGAGTTTCAAAACTCTTCTCAGAATACTCTGCTGCGAAAGAGAGCAGTACAGCTGCCGAGCATGATGCTGAAGCGAAACTTGCAGCTGATGTTGCTCGTGTTGAATCTGGTCAAGACGCGGGTAGAATGAAACAATTGCTGGTGAAGGATCAGGCGGCAAAGGAGATCAAGCCACTAACAGATCGTGTCAATATGGAGCTTACAACGTTGAAGAATTTGAAAAAGGAGTTCATGAGAGTACTTGTTGCTCGATGCCAAGCCAATCAAGACACCGAGGGAGAAGATTCTCTCAGTGGACCAGCTCAAAAGCAACGAATTCAGTTCTTGGAGAACAATTTGGACAAGTTGACGAAGGTTCACAAGCAGgtttgtcgtttttattctcattttgattatcttaaaacttgaattttcagCTTGTTCGCGACAATGCCGATTTGCGCGTTGAACTGCCAAAGATGGAAGCTCGTCTTCGTGGTCGTGAAGATCGCATCAAAATATTAGAAACTGCTCTTCGTGATTCGAAGCAACGTAGTCAAGCAGAACGAAAGAAGTATCAACAAGAAGTTGAACGAATCAAGGAAGCTGTTCGACAACGTAACATGCGACGAATGAATGCTCCACAAATTGTGAAGCCAATCCGTCCAGGACAAGTGTATACGTCTCCGTCAGCAGGAATGTCACAAGGAGCTCCAAATGGCTCAAACGgtgtgtttagtcagacatctacaccttcaacatctcgcaatcagataccatcaaaaatgactatttcacagttgattgcagaaatttaagatttttttaaaaattctttagtgctcatgtaatttttcacaagtaattatactatgaattagaattagagtgagtgtttctttttcttcctaccgtattatcaaatttaacagtcttttgtccgtccatttttcactaatcaaagtttttcagCATAAtgtctcccaacaacaacatcaactcatcgtcttctttgatccaatcaatacactgaagactgacattcaaatgcttctctatctctcttcttttcccggctttgtgatatactttcgatgggcttttctgtttattttaaaatctagtaacttatacaattacgcggcttctggaagtttcaacaaaaatatcttcatttggttggttgtgtctccccatttcgttccttggcttctcgtcttccatgtagaatacaaaacttcaaaagctaaaagtatttaaagcttccctccacccccacccaaattgcctttttccgcctttttgttctaatagtctgtttctatacgattttcctgtttcagttttactaatctgacacgaggttttgtctggttcttccccccgtcacccaccaacactcctatgattgttttttgcatgcgtttgagtgtctttaaagcttgcttgctaaatccccctatcattcttcataagaaatcaacttgtttcgtttctgcacaattcggcccccaaatccccgcacatcccaattg"
    repo_1_aa_sequence = "MDFQNRAGGKTGSGGVASAADAGVDRRERLRQLALETIDLQKDPYFMRNHIGTYECKLCLTLHNNEGSYLAHTQGKKHQANLARRAAKEQSEQPFLPAPQKAAVETKKFVKIGRPGYKVTKERDPGAGQQALLFQIDYPEIADGIAPRHRFMSAYEQKIQPPDKRWQYLLFAAEPYETIGFKIPSREVDKSEKFWTMWNKDTKQFFLQVAFKLERLDDQPYY"
    # amino_acid_sequence = "MEPRTDGAECGVQVFCRIRPLNKTEEKNADRFLPKFPSEDSISLGGKVYVFDKVFKPNTTQEQVYKGAAYHIVQDVLSGYNGTVFAYGQTSSGKTHTMEGVIGDNGLSGIIPRIVADIFNHIYSMDENLQFHIKVSYYEIYNEKIRDLLDPEKVNLSIHEDKNRVPYVKGATERFVGGPDEVLQAIEDGKSNRMVAVTNMNEHSSRSHSVFLITVKQEHQTTKKQLTGKLYLVDLAGSEKVSKTGAQGTVLEEAKNINKSLTALGIVISALAEGTKSHVPYRDSKLTRILQESLGGNSRTTVIICASPSHFNEAETKSTLLFGARAKTIKNVVQINEELTAEEWKRRYEKEKEKNTRLAALLQAAALELSRWRAGESVSEVEWVNLSDSAQMAVSEVSGGSTPLMERSIAPAPPMLTSTTGPITDEEKKKYEEERVKLYQQLDEKDDEIQKVSQELEKLRQQVLLQEEALGTMRENEELIREENNRFQKEAEDKQQEGKEMMTALEEIAVNLDVRQAECEKLKRELEVVQEDNQSLEDRMNQATSLLNAHLDECGPKIRHFKEGIYNVIREFNIADIASQNDQLPDHDLLNHVRIGVSKLFSEYSAAKESSTAAEHDAEAKLAADVARVESGQDAGRMKQLLVKDQAAKEIKPLTDRVNMELTTLKNLKKEFMRVLVARCQANQDTEGEDSLSGPAQKQRIQFLENNLDKLTKVHKQLVRDNADLRVELPKMEARLRGREDRIKILETALRDSKQRSQAERKKYQQEVERIKEAVRQRNMRRMNAPQIVKPIRPGQVYTSPSAGMSQGAPNGSNA"
    cp = CrisprPlanner("repo-1",
                       aa_mutation_site=31,
                       sense_strand=repo_1_sense_strand,
                       amino_acid_sequence=repo_1_aa_sequence)
    cp.plan_my_crispr(from_aa=AminoAcid.ARGININE, to_aa=AminoAcid.GLUTAMINE)
    # cp.plan_my_crispr(aa_mutation_site=26)

# work2 = False
# if work2:
#     # seq = "GCCCCAACAAT"
#     # demands = {4: 'C', 5: 'A', 6: 'A'}
#     # pam_sites = (6, 8)
#     mutants_dic = {}
#     CrisprPlanner.modify_seq_to_change_restriction_sites("TCCAACA", {2: "C", 3: "A", 4: "A"},
#                                                          mutants=mutants_dic, max_mutations=2, mutations_so_far=0)
#     print(mutants_dic)
#     print(str(len(mutants_dic)))

work3 = False
if work3:
    mutated_seq = "CCGGACAA"
    CrisprPlanner.find_restriction_sites("TCCGCCAG")

work4 = False
if work4:
    a = "AAA"
    b = AminoAcid.ASPARAGINE
    codon_details = CrisprPlanner.how_to_get_b_from_a(a, b)
    print(codon_details)

work5 = False
if work5:
    human_gene_name = 'cct-1'
    amino_acid_mutation_site = 287
    nt_seq = "gtaATGGCATCAGCTGGAGATTCCATTCTTGCCCTCACCGGTAAAAGAACTACTGGACAAGGCATCAGATCTCAGAATGgtaacaccgaaagctcaatataagtatacattaattaattgcagTCACCGCGGCAGTTGCGATCGCCAATATTGTGAAGTCATCTCTTGGCCCTGTCGGACTTGATAAAATGCTTGTCGATGATGTTGGAGATGTCATTGTCACAAATGACGGAGCCACAATTCTGAAACAACTCGAGGTTGAGCATCCGGCTGGAAAAGTGCTTGTAGAACTTGCACAGCTGCAAGACGAGGAGGTCGGAGATGGAACTACTTCTGTCGTTATTGTGGCGGCTGAGCTCTTGAAGAGAGCCGATGAGCTTGTGAAACAAAAAGTTCATCCGACGACTATTATCAATGGTTACCGTCTCGCGTGCAAGGAAGCCGTCAAGTACATTAGTGAAAACATCTCATTCACTTCCGACTCGATTGGTAGACAATCAGTTGTCAACGCTGCCAAAACTTCCATGAGCAGTAAGATTATCGGACCgtgagtttggtgttgtctatgcttcaagaaaattgatttttcagAGACGCCGATTTCTTCGGAGAGCTGGTTGTTGATGCCGCGGAAGCTGTTCGTGTGGAAAATAACGGGAAAGTCACTTATCCTATCAATGCAGTCAATGTTCTGAAGGCCCACGGAAAGAGCGCTCGCGAATCAGTTTTGGTGAAAGGATATGCACTCAATTGCACAGTTGCCAGTCAGGCCATGCCACTTCGTGTTCAAAATGCCAAGATCGCATGTCTCGATTTCTCTTTGATGAAGGCTAAGATGCACCTCGGTATTTCAGTCGTTGTTGAAGATCCAGCCAAGCTTGAGGCTATTCGCAGAGAgtgagttgaaactattcgtttctttttaagctatggaattttcagAGAATTCGATATTACCAAACGCCGCATTGATAAAATTTTGAAAGCCGGAGCCAACGTTGTTCTTACAACTGGAGGTATCGATGATTTGTGCTTGAAGCAATTTGTCGAATCTGGAGCTATGGCTGTTCGTCGATGCAAGAAATCAGACTTGAAGAGAATTGCCAAAGCTACTGGAGCCACATTGACTGTTTCCTTGGCTACTTTGGAAGGAGATGAAGCTTTCGATGCCTCGCTTCTTGGACATGCCGATGAAATTGTTCAAGAAAGAATTAGTGACGACGAGCTCATTCTCATCAAGGGACCGAAATCTCGTACTGCCAGCAGCATTATCCTCCGTGGAGCGAACGATGTGATGCTCGATGAAATGGAGAGATCGGTTCACGACTCACTCTGTGTTGTTCGTAGAGTTCTGGAAAGCAAGAAACTTGTGGCTGGAGGAGGTGCTGTTGAGACTTCTCTCAGTCTTTTCCTTGAAACTTATGCACAAACCTTGTCTTCTCGCGAGCAGCTTGCTGTTGCTGAATTCGCTTCAGCGCTTCTCATCATTCCGAAGGTTTTGGCAAGCAATGCTGCAAGAGATTCTACTGATTTAGTGACAAAACTCCGCGCGTACCACTCCAAAGCTCAATTGATCCCACAACTTCAACACCTCAAGTGgtaagtgaaaatgttttttttaaagagtaggttattacatgttagcttaatgtaataaaattaaaataatttatttcaaaaaatttcgttttgtgcttagaaaaagcgtctaattcatgttttctgaatttgagtcagtttattcactctttttttagGGCTGGTTTGGATCTCGAAGAAGGCACGATCCGCGATAACAAGGAGGCTGGAATTTTGGAGCCAGCTCTTAGTAAGGTCAAGTCTCTGAAGTTCGCCACTGAGGCAGCCATTACGATATTGCGTATTGATGACCTCATCAAACTTGACAAGCAAGAGCCACTTGGAGGAGATGATTGCCACGCTTAAattttcccgtttaccccgtttatatatccctgttttccgcgtgcttctcacataattccgatctgctgctccttatcccaaattctcatgttcagcttttgttttcttcttttgatgatactttattgaacgaaatgttgtaagttttaatgttttgatttcaaagttgtttgtattcgtttttcattattcaaacaatgaagaagctttgccac"
    CrisprPlanner(gene_name=human_gene_name, aa_mutation_site=amino_acid_mutation_site, sense_strand=nt_seq). \
        plan_my_crispr(from_aa=AminoAcid.ASPARAGINE,
                       to_aa=AminoAcid.SERINE,
                       check_consistency=True)
work6 = False
if work6:
    print(CrisprPlanner.how_to_get_b_from_a("ACG", AminoAcid.ALANINE))

work7 = True
if work7:
    start = time.time()
    fin = open(r"C:\Users\RZBlab\PycharmProjects\Research-RZB\Code\CRISPR\restriction_enzymes.txt")
    for line in fin:
        lst = line.rstrip("\n").split("\t")
        enzyme_site = lst[0]
        enzyme_name = lst[1]
        parsed_sites_rec = FileReader.find_enzymes_derivatives(enzyme_site, rec=True)
    fin.close()
    end = time.time()
    print(end - start)
    start = time.time()
    fin = open(r"C:\Users\RZBlab\PycharmProjects\Research-RZB\Code\CRISPR\restriction_enzymes.txt")
    for line in fin:
        lst = line.rstrip("\n").split("\t")
        enzyme_site = lst[0]
        enzyme_name = lst[1]
        parsed_sites_iter = FileReader.find_enzymes_derivatives(enzyme_site, rec=False)
    fin.close()
    end = time.time()
    print(end - start)


work8 = False
if work8:
    human_gene_name = 'cct-1'
    amino_acid_mutation_site = 287
    nt_seq = "gtaATGGCATCAGCTGGAGATTCCATTCTTGCCCTCACCGGTAAAAGAACTACTGGACAAGGCATCAGATCTCAGAATGgtaacaccgaaagctcaatataagtatacattaattaattgcagTCACCGCGGCAGTTGCGATCGCCAATATTGTGAAGTCATCTCTTGGCCCTGTCGGACTTGATAAAATGCTTGTCGATGATGTTGGAGATGTCATTGTCACAAATGACGGAGCCACAATTCTGAAACAACTCGAGGTTGAGCATCCGGCTGGAAAAGTGCTTGTAGAACTTGCACAGCTGCAAGACGAGGAGGTCGGAGATGGAACTACTTCTGTCGTTATTGTGGCGGCTGAGCTCTTGAAGAGAGCCGATGAGCTTGTGAAACAAAAAGTTCATCCGACGACTATTATCAATGGTTACCGTCTCGCGTGCAAGGAAGCCGTCAAGTACATTAGTGAAAACATCTCATTCACTTCCGACTCGATTGGTAGACAATCAGTTGTCAACGCTGCCAAAACTTCCATGAGCAGTAAGATTATCGGACCgtgagtttggtgttgtctatgcttcaagaaaattgatttttcagAGACGCCGATTTCTTCGGAGAGCTGGTTGTTGATGCCGCGGAAGCTGTTCGTGTGGAAAATAACGGGAAAGTCACTTATCCTATCAATGCAGTCAATGTTCTGAAGGCCCACGGAAAGAGCGCTCGCGAATCAGTTTTGGTGAAAGGATATGCACTCAATTGCACAGTTGCCAGTCAGGCCATGCCACTTCGTGTTCAAAATGCCAAGATCGCATGTCTCGATTTCTCTTTGATGAAGGCTAAGATGCACCTCGGTATTTCAGTCGTTGTTGAAGATCCAGCCAAGCTTGAGGCTATTCGCAGAGAgtgagttgaaactattcgtttctttttaagctatggaattttcagAGAATTCGATATTACCAAACGCCGCATTGATAAAATTTTGAAAGCCGGAGCCAACGTTGTTCTTACAACTGGAGGTATCGATGATTTGTGCTTGAAGCAATTTGTCGAATCTGGAGCTATGGCTGTTCGTCGATGCAAGAAATCAGACTTGAAGAGAATTGCCAAAGCTACTGGAGCCACATTGACTGTTTCCTTGGCTACTTTGGAAGGAGATGAAGCTTTCGATGCCTCGCTTCTTGGACATGCCGATGAAATTGTTCAAGAAAGAATTAGTGACGACGAGCTCATTCTCATCAAGGGACCGAAATCTCGTACTGCCAGCAGCATTATCCTCCGTGGAGCGAACGATGTGATGCTCGATGAAATGGAGAGATCGGTTCACGACTCACTCTGTGTTGTTCGTAGAGTTCTGGAAAGCAAGAAACTTGTGGCTGGAGGAGGTGCTGTTGAGACTTCTCTCAGTCTTTTCCTTGAAACTTATGCACAAACCTTGTCTTCTCGCGAGCAGCTTGCTGTTGCTGAATTCGCTTCAGCGCTTCTCATCATTCCGAAGGTTTTGGCAAGCAATGCTGCAAGAGATTCTACTGATTTAGTGACAAAACTCCGCGCGTACCACTCCAAAGCTCAATTGATCCCACAACTTCAACACCTCAAGTGgtaagtgaaaatgttttttttaaagagtaggttattacatgttagcttaatgtaataaaattaaaataatttatttcaaaaaatttcgttttgtgcttagaaaaagcgtctaattcatgttttctgaatttgagtcagtttattcactctttttttagGGCTGGTTTGGATCTCGAAGAAGGCACGATCCGCGATAACAAGGAGGCTGGAATTTTGGAGCCAGCTCTTAGTAAGGTCAAGTCTCTGAAGTTCGCCACTGAGGCAGCCATTACGATATTGCGTATTGATGACCTCATCAAACTTGACAAGCAAGAGCCACTTGGAGGAGATGATTGCCACGCTTAAattttcccgtttaccccgtttatatatccctgttttccgcgtgcttctcacataattccgatctgctgctccttatcccaaattctcatgttcagcttttgttttcttcttttgatgatactttattgaacgaaatgttgtaagttttaatgttttgatttcaaagttgtttgtattcgtttttcattattcaaacaatgaagaagctttgccac"
    restriction_sites = CrisprPlanner(gene_name=human_gene_name,
                                      aa_mutation_site=amino_acid_mutation_site,
                                      sense_strand=nt_seq).restriction_enzymes
    print(restriction_sites)

work9 = False
if work9:
    human_gene_name = 'cct-1'
    amino_acid_mutation_site = 287
    nt_seq = "gtaATGGCATCAGCTGGAGATTCCATTCTTGCCCTCACCGGTAAAAGAACTACTGGACAAGGCATCAGATCTCAGAATGgtaacaccgaaagctcaatataagtatacattaattaattgcagTCACCGCGGCAGTTGCGATCGCCAATATTGTGAAGTCATCTCTTGGCCCTGTCGGACTTGATAAAATGCTTGTCGATGATGTTGGAGATGTCATTGTCACAAATGACGGAGCCACAATTCTGAAACAACTCGAGGTTGAGCATCCGGCTGGAAAAGTGCTTGTAGAACTTGCACAGCTGCAAGACGAGGAGGTCGGAGATGGAACTACTTCTGTCGTTATTGTGGCGGCTGAGCTCTTGAAGAGAGCCGATGAGCTTGTGAAACAAAAAGTTCATCCGACGACTATTATCAATGGTTACCGTCTCGCGTGCAAGGAAGCCGTCAAGTACATTAGTGAAAACATCTCATTCACTTCCGACTCGATTGGTAGACAATCAGTTGTCAACGCTGCCAAAACTTCCATGAGCAGTAAGATTATCGGACCgtgagtttggtgttgtctatgcttcaagaaaattgatttttcagAGACGCCGATTTCTTCGGAGAGCTGGTTGTTGATGCCGCGGAAGCTGTTCGTGTGGAAAATAACGGGAAAGTCACTTATCCTATCAATGCAGTCAATGTTCTGAAGGCCCACGGAAAGAGCGCTCGCGAATCAGTTTTGGTGAAAGGATATGCACTCAATTGCACAGTTGCCAGTCAGGCCATGCCACTTCGTGTTCAAAATGCCAAGATCGCATGTCTCGATTTCTCTTTGATGAAGGCTAAGATGCACCTCGGTATTTCAGTCGTTGTTGAAGATCCAGCCAAGCTTGAGGCTATTCGCAGAGAgtgagttgaaactattcgtttctttttaagctatggaattttcagAGAATTCGATATTACCAAACGCCGCATTGATAAAATTTTGAAAGCCGGAGCCAACGTTGTTCTTACAACTGGAGGTATCGATGATTTGTGCTTGAAGCAATTTGTCGAATCTGGAGCTATGGCTGTTCGTCGATGCAAGAAATCAGACTTGAAGAGAATTGCCAAAGCTACTGGAGCCACATTGACTGTTTCCTTGGCTACTTTGGAAGGAGATGAAGCTTTCGATGCCTCGCTTCTTGGACATGCCGATGAAATTGTTCAAGAAAGAATTAGTGACGACGAGCTCATTCTCATCAAGGGACCGAAATCTCGTACTGCCAGCAGCATTATCCTCCGTGGAGCGAACGATGTGATGCTCGATGAAATGGAGAGATCGGTTCACGACTCACTCTGTGTTGTTCGTAGAGTTCTGGAAAGCAAGAAACTTGTGGCTGGAGGAGGTGCTGTTGAGACTTCTCTCAGTCTTTTCCTTGAAACTTATGCACAAACCTTGTCTTCTCGCGAGCAGCTTGCTGTTGCTGAATTCGCTTCAGCGCTTCTCATCATTCCGAAGGTTTTGGCAAGCAATGCTGCAAGAGATTCTACTGATTTAGTGACAAAACTCCGCGCGTACCACTCCAAAGCTCAATTGATCCCACAACTTCAACACCTCAAGTGgtaagtgaaaatgttttttttaaagagtaggttattacatgttagcttaatgtaataaaattaaaataatttatttcaaaaaatttcgttttgtgcttagaaaaagcgtctaattcatgttttctgaatttgagtcagtttattcactctttttttagGGCTGGTTTGGATCTCGAAGAAGGCACGATCCGCGATAACAAGGAGGCTGGAATTTTGGAGCCAGCTCTTAGTAAGGTCAAGTCTCTGAAGTTCGCCACTGAGGCAGCCATTACGATATTGCGTATTGATGACCTCATCAAACTTGACAAGCAAGAGCCACTTGGAGGAGATGATTGCCACGCTTAAattttcccgtttaccccgtttatatatccctgttttccgcgtgcttctcacataattccgatctgctgctccttatcccaaattctcatgttcagcttttgttttcttcttttgatgatactttattgaacgaaatgttgtaagttttaatgttttgatttcaaagttgtttgtattcgtttttcattattcaaacaatgaagaagctttgccac"
    crisprPlanner = CrisprPlanner(gene_name=human_gene_name,
                                      aa_mutation_site=amino_acid_mutation_site,
                                      sense_strand=nt_seq)
    print("finished with the restriction sites list")
    crisprPlanner.find_sites("TTGCA")
    crisprPlanner.find_sites("TCAGA")
    crisprPlanner.find_sites("GATTC")

