import random
import re
from statistics import mean

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

    def __init__(self, gene_name, aa_mutation_site, sense_strand, amino_acid_sequence: str = "",
                 existing_enzymes=None):
        self.gene_name = gene_name
        self.sense_strand = sense_strand
        self.anti_sense_strand = CrisprPlanner.get_complementary_sequence(sense_strand)
        self.amino_acid_sequence = amino_acid_sequence if amino_acid_sequence else \
            BioPython().get_aa_seq_by_c_elegans_gene_name(gene_name)
        self.amino_acid_mutation_site = aa_mutation_site
        self.sense_mutation_site = -1
        self.anti_sense_mutation_site = -1
        self.ssODN_mutation_codon_start = -1
        self.mutated_strand = None
        self.mutation_direction = None
        self.DSB = None
        self.pam_sites = None
        self.ssODN_direction = None
        self.mutated_sites = []
        self.restriction_enzymes = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Code\CRISPR",
                                              r"\parsed_restriction_enzymes.txt").get_parsed_restriction_enzymes_list()
        self.relevant_restriction_enzymes = self.get_relevant_restriction_enzymes(existing_enzymes) if existing_enzymes else self.restriction_enzymes

    def get_relevant_strand(self, direction):
        if direction > 0:
            return self.sense_strand
        else:
            return self.anti_sense_strand

    def plan_my_crispr(self,
                       from_aa: AminoAcid,
                       to_aa: AminoAcid,
                       check_consistency: bool = True,
                       window_size: int = 30,
                       PAM_size: int = 3):

        self.initiate_crispr(check_consistency)

        chosen_crrna, crrna_strand = self.get_crrna(window_size, PAM_size)

        if not chosen_crrna:
            return None
        # now we have our cr_rna
        print("chosen crRNA: " + str(chosen_crrna))
        strand_str = "anti-" if crrna_strand < 0 else ""
        print("The relevant strand is " + strand_str + "sense")

        self.complete_fields(crrna_strand, to_aa, chosen_crrna)

        first_strand = self.mutated_strand
        first_mutated_sites = self.mutated_sites
        # first try - when at last we try to insert restriction site
        insertion_result = self.insert_mutations(self.ssODN_mutation_codon_start, from_aa, to_aa,
                                                 RestrictionSiteType.INSERTED)
        if not insertion_result:
            self.mutated_strand = first_strand
            self.mutated_sites = first_mutated_sites
            deletion_result = self.insert_mutations(self.ssODN_mutation_codon_start, from_aa, to_aa,
                                                    RestrictionSiteType.REMOVED)
            if not deletion_result:
                print("Could not insert or remove restriction site")
                return False
        return True

    # initiates the CRISPR planner. fills in the amino_acid_sequence if not already filled in, and finds the nt site in
    # sense strand and anti-sense strand
    def initiate_crispr(self, check_consistency: bool):
        print("Gene's sense sequence:", self.sense_strand)
        if not self.amino_acid_sequence:
            # amino acid couldn't be extracted
            self.amino_acid_sequence = input("Insert amino acid sequence for gene " + self.gene_name + ":")
            if not self.amino_acid_sequence:
                print("Exception in initiate_crispr: no amino acid sequence was delivered or could be extracted")
                return None
        print("Gene's amino acid sequence: ", self.amino_acid_sequence)
        self.sense_mutation_site = self.find_site_in_nt_seq(amino_acid_site=self.amino_acid_mutation_site,
                                                            check_sequence_consistency=check_consistency)
        print("mutation site in sense strand is: " + str(self.sense_mutation_site), ", and the codon is:",
              self.sense_strand[self.sense_mutation_site:self.sense_mutation_site + 3])
        self.anti_sense_mutation_site = self.get_mutation_site_for_anti_sense(self.sense_mutation_site)
        print("mutation site in anti-sense strand: " + str(self.anti_sense_mutation_site))

    # start_crrna_index = for testing purposes
    def get_crrna(self, window_size, PAM_size, start_crrna_index=None):
        # get optional sequences for crRNA
        sense_strand_sequences = self.get_list_of_potential_sequences(self.sense_strand,
                                                                      self.sense_mutation_site,
                                                                      window_size,
                                                                      PAM_size)
        anti_sense_strand_sequences = self.get_list_of_potential_sequences(self.anti_sense_strand,
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

        return self.choose_best_crrna(sense_strand_sequences, anti_sense_strand_sequences, start_crrna_index)

    def complete_fields(self, crrna_strand_direction, to_aa, chosen_crrna):
        self.pam_sites = SequenceSites(start=chosen_crrna[1][1] + 1, end=chosen_crrna[1][1] + 3)

        self.ssODN_direction, self.mutation_direction = self.choose_ssODN_strand(crrna_strand_direction, to_aa)
        print("Mutation direction is:", str(self.mutation_direction), "thus ssODN direction is:",
              str(self.ssODN_direction))
        if self.ssODN_direction != crrna_strand_direction:
            print("The ssODN is on the opposite strand!")
            # ssODN is not on the strand from which we got the crRNA, thus we need to find its pam sites again
            self.pam_sites = SequenceSites(self.find_index_in_parallel_strand(self.pam_sites.end),
                                           self.find_index_in_parallel_strand(self.pam_sites.start))
            print("pam sites on the ssODN strand:", self.pam_sites)
        self.DSB = self.get_dsb(self.pam_sites, self.mutation_direction)

        # to be mutated
        self.mutated_strand = self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand
        self.ssODN_mutation_codon_start = self.sense_mutation_site if self.ssODN_direction > 0 else self.anti_sense_mutation_site

        # getting the zone where all mutations will be located: after to before DSB
        # sanity check
        print("mutation codon is at", self.ssODN_mutation_codon_start)

    # this function is responsible for inserting all needed mutations to the string: to change amino acid, to prevent
    # reattachment and to add\remove restriction site
    def insert_mutations(self, ssODN_codon_mutation_site, from_aa, to_aa, restriction_site_type):
        # 1. mutation to change codon to change the amino acid
        possible_codon_mutations = self.get_possible_mutations_demands(from_aa, to_aa, self.sense_mutation_site)
        second_strand = self.mutated_strand
        second_mutated_sites = self.mutated_sites
        for possible_codon_mutation in possible_codon_mutations:
            print("Now working on codon mutation:", possible_codon_mutation)
            self.mutated_strand = second_strand
            self.mutated_sites = second_mutated_sites
            self.apply_codon_mutation(ssODN_codon_mutation_site,
                                      possible_codon_mutation.dict_of_mutations)

            # making sure that the strand is good so far
            print("mutated strand after adding codon mutation:", self.mutated_strand)
            print("sites changed:", self.mutated_sites)

            # 2. mutations to change nt to prevent re-attachments - in PAM or crRNA sequence
            # if mutation is DOWNSTREAM the section to be changed is PAM site, if UPSTREAM - crRNA.
            if self.mutation_direction == MutationDirection.DOWNSTREAM:
                section_to_mutate = ReattachmentSection(1,
                                                        DNASection.PAM_SITE,
                                                        self.pam_sites)
            else:
                section_to_mutate = ReattachmentSection(4,
                                                        DNASection.CR_RNA,
                                                        SequenceSites(self.pam_sites.start - 20, self.pam_sites.start - 1 - 3))

            print("section to mutate to prevent re-attachment:", section_to_mutate.section_sites, "(close segment) in",
                  section_to_mutate.section_type)

            possible_mutations, number_of_mutants = self.get_possible_codon_mutations(section_to_mutate,
                                                                                      ssODN_codon_mutation_site)
            # subsets with indexes that represent mutations that sum up to enough mutations overall
            valid_index_subsets = CrisprPlanner.get_valid_subsets(number_of_mutants, possible_mutations, self.pam_sites,
                                                                  self.mutation_direction)
            third_strand = self.mutated_strand
            third_mutated_sites = self.mutated_sites
            print("all valid index subsets:", valid_index_subsets)
            for index_subset in valid_index_subsets:
                self.mutated_strand = third_strand
                self.mutated_sites = third_mutated_sites
                point_mutations = self.get_point_mutations(possible_mutations,
                                                           index_subset)  # also mutates the sequence
                print("now working on index subset:", index_subset, "and point mutations:", point_mutations)
                print("mutated strand after adding anti-reattachment mutations:", self.mutated_strand)
                self.mutated_sites.extend(point_mutations)
                print("all mutated sites:", self.mutated_sites)

                # 3. mutations to add/remove restriction sites
                print("Add or Remove restriction sites:")
                result = self.add_remove_restriction_sites(restriction_site_type)
                if result:
                    print("Restriction site chosen:", result)
                    return True
        return False

    # receives (1) number of mutations needed, (2) the list of all possible anti-reattachment section mutations, and
    # (3) power set of all mutations indexes, and leave out only the ones with sufficient amount of mutations. if there
    # are not any, leaves only subsets with maximum number of mutations
    @staticmethod
    def get_valid_subsets(number_of_mutations, possible_mutations, pam_sites, mutation_direction):
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
                    # check if mutations don't override each other by checking codon sites
                    if not CrisprPlanner.are_codons_overriden(subset, possible_mutations):
                        valid_subsets.append((subset, mutations_in_subset, sum_usage))
        if valid_subsets != [([], 0, 0)]:
            # sort subsets according to number of mutations max usage
            print("valid subsets:", valid_subsets)
            # sort according to minimal number of mutations, minimal distance from the DSB and maximal usage percentage
            valid_subsets.sort(key=lambda mutation: (mutation[1],
                                                     CrisprPlanner.get_group_distance_from_dsb(mutation[0],
                                                                                               possible_mutations,
                                                                                               pam_sites,
                                                                                               mutation_direction),
                                                     CrisprPlanner.get_avg_distance_from_dsb(mutation[0],
                                                                                             possible_mutations,
                                                                                             pam_sites,
                                                                                             mutation_direction),
                                                     - mutation[2]))
            print("valid subsets after sort:", valid_subsets)
            return list(map(lambda item: item[0], valid_subsets))
        else:
            print("WARNING: the maximum amount of mutations is", max_underscored, ", there are no subsets with",
                  number_of_mutations, "mutations")
            # sort according to number of mutations first, and codons usage second
            underscored_subsets.sort(key=lambda mutation: (mutation[1], -mutation[2]))
            return list(map(lambda item: item[0], underscored_subsets))

    # return the distance of the average point in subset from DSB
    @staticmethod
    def get_avg_distance_from_dsb(subset, possible_mutations, pam_sites, mutation_direction):
        DSB = CrisprPlanner.get_dsb(pam_sites, mutation_direction)
        mean_point = 0
        count = 0
        for index in subset:
            start = possible_mutations[index][0].codon_sites.start
            for mutation_inner_index in possible_mutations[index][1].dict_of_mutations:
                mean_point += start + mutation_inner_index
                count += 1
        mean_point /= count
        return abs(mean_point-DSB)

    # returns the distance of the furthest point from the DSB
    @staticmethod
    def get_group_distance_from_dsb(subset, possible_mutations, pam_sites, mutation_direction):
        DSB = CrisprPlanner.get_dsb(pam_sites, mutation_direction)
        max_distance = 0
        for index in subset:
            start = possible_mutations[index][0].codon_sites.start
            average_mutated_index = mean([start+item for item in possible_mutations[index][1].dict_of_mutations])
            if abs(average_mutated_index-DSB) > max_distance:
                max_distance = abs(average_mutated_index-DSB)
        return max_distance

    @staticmethod
    def get_dsb(pam_sites, mutation_direction):
        if mutation_direction == MutationDirection.DOWNSTREAM:
            DSB = pam_sites.end + 4
        else:
            DSB = pam_sites.start - 3
        return DSB

    # checks if the subset points to mutations that try to mutate the same codons
    @staticmethod
    def are_codons_overriden(subset, possible_mutations):
        codons = set()
        for index in subset:
            mutation = possible_mutations[index]
            codon_data = mutation[0]
            if codon_data in codons:
                return True
            codons.add(codon_data)
        return False

    # receives (1) the section that needs to be mutated to prevent reattachment, (2) the indexes of the amino acid's
    # mutation codon on the ssODN strand,(3) the mutation direction, and (4) the sites that have
    #  changed while mutating the amino acid, computes which indexes can be changed in order to insert the section with
    # silent mutations and returns the possibilities for silent mutations
    def get_possible_codon_mutations(self,
                                     section_to_mutate: ReattachmentSection):
        # format of codons_data: ('GCG', (1, 3))
        codons_data = CrisprPlanner.get_relevant_codons(section_to_mutate, self.ssODN_mutation_codon_start, self.mutated_strand)
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
            # removes codons in which the difference is in part that in outside of the mutation zone
            CrisprPlanner.check_outside_codon_mutations(codon_data, same_aa_codons, section_to_mutate)

            optional_codon_mutations = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(codon_data.codon,
                                                                                            same_aa_codons)
            for optional_mutation in optional_codon_mutations:
                possible_mutations.append((codon_data, optional_mutation))
        # sort mutations according to which is closest the the mean of former point mutations
        possible_mutations.sort(key=lambda x: abs(x[0].codon_sites.start+mean(x[1].dict_of_mutations.keys())-self.DSB))
        print("possible mutations:", *possible_mutations, sep="\n")
        number_of_mutants = CrisprPlanner.get_number_of_mutants(section_to_mutate, self.mutated_sites)
        return possible_mutations, number_of_mutants

    # receives (1) codon data (2) list of codons that are translated into same amino acid as the codon, and (3) the
    # section to mutate, and removes all codons from the same_aa_codons_ list in which the difference in nucleotides
    # is outside of the section to mutate. Also if the section is crRNA we leave out all codons that write on the DSB.
    @staticmethod
    def check_outside_codon_mutations(codon_data, same_aa_codons, section_to_mutate: ReattachmentSection):
        codon_start = codon_data.codon_sites.start
        for same_aa_codon in same_aa_codons[:]:
            _, demands = CrisprPlanner.codon_distance(codon_data.codon, same_aa_codon)
            print(same_aa_codon, ":", demands)
            valid_demands = False
            dsb_safe = True
            for index in demands:
                if section_to_mutate.section_sites.start <= codon_start + index <= section_to_mutate.section_sites.end:
                    valid_demands = True
                if section_to_mutate.section_type == DNASection.CR_RNA:
                    if section_to_mutate.section_sites.end < codon_start + index:
                        dsb_safe = False
            if not valid_demands:
                print("codon:", same_aa_codon, "is removed because of outside mutations")
                same_aa_codons.remove(same_aa_codon)
            if same_aa_codon in same_aa_codons and not dsb_safe:
                print("codon:", same_aa_codon, "is removed because of writing over DSB")
                same_aa_codons.remove(same_aa_codon)

    # receives (1) all possible mutations to prevent reattachment and (2) list of indexes corresponding to the mutations
    # list and returns the point mutations in format of namedTuples: PointMutations(index, new nucleotide)
    def get_point_mutations(self, possible_mutations, indexes_subset):
        # possible mutation format: (CodonData(codon='TTT', codon_sites=SequenceSites(start=1174, end=1176)), OptionalCodonChange('TTC', 1, {2: 'C'}, 23.9))
        point_mutations = []
        for index in indexes_subset:
            possible_mutation = possible_mutations[index]
            dict_of_changes = possible_mutation[1].dict_of_mutations
            for key in dict_of_changes:
                point_mutations.append(
                    PointMutation(possible_mutation[0].codon_sites.start + key, dict_of_changes[key]))
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
            if section_to_mutate.section_sites.start <= mutated_site.index <= section_to_mutate.section_sites[1] + extra:
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
                    print(
                        "all mutations of this sort must be erased since mutation in that codon will leave out the NGG")
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
    def get_relevant_codons(section_to_mutate: ReattachmentSection, ssODN_mutation_codon_start, mutated_strand):
        codons_data = []
        remainder = abs(section_to_mutate.section_sites.start - ssODN_mutation_codon_start) % 3
        if section_to_mutate.section_sites.start >= ssODN_mutation_codon_start:
            start_point = section_to_mutate.section_sites.start - remainder
        else:
            if remainder > 0:
                start_point = section_to_mutate.section_sites.start - (3 - remainder)
            else:
                start_point = section_to_mutate.section_sites.start
        if section_to_mutate.section_type == DNASection.CR_RNA:
            end_point = section_to_mutate.section_sites.end - 2  # so codons won't cross the DSB
        else:
            end_point = section_to_mutate.section_sites.end  # codons can cross the end of section
        while start_point <= end_point:  # far enough so a codon won't cross the end:
            codon = mutated_strand[start_point:start_point + 3]
            codons_data.append(CodonData(codon, SequenceSites(start_point, start_point + 2)))
            start_point += 3
        print(str(len(codons_data)), "relevant codons that we can mutate:", codons_data)
        return codons_data

    # This function receives (1) an amino acid sequence site and (2) value that indicates whether to check the
    # correlation between amino acid sequence and nucleotide sequence (codon to amino acid) or not, and returns the
    # site of the mutation in the nucleotide sequence, meaning the site or the corresponding codon
    def find_site_in_nt_seq(self, amino_acid_site, check_sequence_consistency: bool = True):
        no_utr_strand = self.remove_utr()
        print("sense without utr:", no_utr_strand)
        if check_sequence_consistency:
            # consistency check
            nt_index = 0
            aa_index = 0
            while aa_index < len(self.amino_acid_sequence):
                CrisprPlanner.check_sequence_consistency(no_utr_strand[nt_index:nt_index + 3],
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
        for i in range(start_point, end_point):
            if strand[i + 1:i + 3] == "GG":
                try:
                    potential_cr_rna = strand[i - cr_rna_size:i]
                    potential_cr_rnas.append((potential_cr_rna, SequenceSites(i - cr_rna_size, i-1)))
                except Exception as e:
                    print("Not enough space for a whole crRNA, i = " + str(i))
                    print("Exception in get_list_of_potential_sequences:", e)
                    continue
        return potential_cr_rnas

    # receives (1) the mutation position in the sense strand, and returns the position of this mutation in the
    # anti-sense
    def get_mutation_site_for_anti_sense(self, sense_mutation_site):
        return len(self.anti_sense_strand) - sense_mutation_site - 3

    # receives (1) a sequence and returns its complementary, from 5' to 3'
    @staticmethod
    def get_complementary_sequence(given_sense_strand):
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
            anti_sense_strand = CrisprPlanner.get_complementary_sequence(sense_strand)
            complementary_strands.append(anti_sense_strand)
        return complementary_strands

    # this function gets a nucleotide sequence and adds to a dictionary all restriction sites in the sequence
    def find_restriction_sites(self, nucleotide_seq, start_point, methylated: bool = False):
        restriction_sites_in_seq = []
        for restriction_enzyme in self.restriction_enzymes:
            for restriction_site in restriction_enzyme.derivatives:
                opposite_restriction_site = self.get_complementary_sequence(restriction_site)
                sites = [restriction_site, opposite_restriction_site]
                for site in sites:
                    indexes = [m.start() for m in re.finditer('(?=' + site + ')', nucleotide_seq)]
                    for index in indexes:
                        if restriction_site == "GGTGA":
                            if index + 8 >= len(nucleotide_seq):
                                continue
                        seq_sites = SequenceSites(start_point + index, start_point + index + len(site) - 1)  # closed section
                        restriction_sites_in_seq.append(RestrictionSite(seq_sites, restriction_enzyme))
        print("Restriction Sites:", *restriction_sites_in_seq, sep='\n')
        return restriction_sites_in_seq

    @staticmethod
    def change_char_in_string(seq, position, new_char):
        lst = list(seq)
        lst[position] = new_char
        return "".join(lst)

    @staticmethod
    def choose_best_crrna(sense_options: list, anti_sense_options: list, start_crrna_index):
        if not start_crrna_index:
            num = random.uniform(0, 1)
            if num > 0.5:
                return random.choice(sense_options), 1
            else:
                return random.choice(anti_sense_options), -1
        else:  # TESTING
            for crrna in sense_options:
                if crrna[1].start == start_crrna_index:
                    return crrna, 1
            for crrna in anti_sense_options:
                if crrna[1].start == start_crrna_index:
                    return crrna, -1
        return None

    def choose_ssODN_strand(self, strand_direction, to_aa: AminoAcid):

        def check_inside_dsb():
            downstream = 0
            upstream = 0
            codon_mutation = CrisprPlanner.how_to_get_b_from_a(crrna_strand[mutation_site:mutation_site+3], to_aa)
            for key in codon_mutation.dict_of_mutations:
                if mutation_site + key >= DSB:
                    downstream += 1
                else:
                    upstream += 1
            if downstream > upstream:
                return MutationDirection.DOWNSTREAM
            return MutationDirection.UPSTREAM

        mutation_site = self.sense_mutation_site if strand_direction > 0 else self.anti_sense_mutation_site
        crrna_strand = self.sense_strand if strand_direction > 0 else self.anti_sense_strand

        DSB = self.pam_sites.start - 3
        # default assignment
        mutation_direction = MutationDirection.UPSTREAM
        print("mutation site: " + str(mutation_site) + ", pam start: " + str(self.pam_sites.start) + ", DSB: " + str(DSB))
        if mutation_site <= DSB <= mutation_site + 2:
            # DSB is inside the codon we need to mutate, further investigation is needed
            mutation_direction = check_inside_dsb()
        elif mutation_site >= DSB:
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
        codon_mutation = CodonMutation(max_codon, codon_options[max_codon][0], codon_options[max_codon][1], max_usage)
        # print(codon_mutation)
        return codon_mutation

    # calculated the demands to go through the codon from the mutation site to one of the codons translating to the
    # wanted amino acid
    def get_possible_mutations_demands(self, from_amino_acid: AminoAcid, to_amino_acid: AminoAcid, mutation_site):
        sense_codon = self.sense_strand[mutation_site:mutation_site + 3]
        codon_options = self.codon_dic[from_amino_acid.value]
        if sense_codon not in codon_options:
            print(f'Codon in mutation site {mutation_site}: {sense_codon} does not translate into {from_amino_acid}')
            return []
        a_codon = sense_codon if self.ssODN_direction > 0 else CrisprPlanner.get_complementary_sequence(sense_codon)
        b_codons = CrisprPlanner.codon_dic[to_amino_acid.value]
        if self.ssODN_direction > 0:
            b_codons_info = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(a_codon, list(b_codons))
        else:
            b_codons_info = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(a_codon,
                                                                                 CrisprPlanner.get_anti_sense_strands(
                                                                                     b_codons))
        print(len(b_codons_info), "possible codons change:", b_codons_info)
        print("The best option is to change", a_codon, "to", b_codons_info[0].codon, "with demands:",
              str(b_codons_info[0].dict_of_mutations))

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
    def check_if_pam_site_mutated(mutated_sites: list, pam_sites: SequenceSites):
        pam_site_mutated = False
        for mutated_site in mutated_sites:
            if pam_sites.start <= mutated_site <= pam_sites.end:
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
        DSB = CrisprPlanner.get_dsb(pam_sites, mutation_direction)
        return SequenceSites(max(0, DSB - 13), DSB-1)  # closed segment

    @staticmethod
    def is_within_mutation_zone(mutation_zone, place):
        if mutation_zone[0] < place < mutation_zone[1]:
            return True
        else:
            return False

    # extracts the part of original and mutated sequences and checks to see whether one of them has a restriction site
    # the other one doesn't. returns the distinctive restriction sites if the lists don't match, and a None value
    # if they do
    def get_distinctive_restriction_sites(self):
        # first we'll find the boundaries in which to search
        lowest_site = min(list(map(lambda item: item.index, self.mutated_sites)))
        highest_site = max(list(map(lambda item: item.index, self.mutated_sites)))
        # max_restriction_site = max(list(map(lambda item: len(item.site), self.restriction_enzymes)))
        max_restriction_site = 30  # prior analyze

        original_sequence = self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand
        original_section = original_sequence[lowest_site - max_restriction_site:highest_site + max_restriction_site]
        mutated_section = self.mutated_strand[lowest_site - max_restriction_site:highest_site + max_restriction_site]

        restrictions_sites_in_original_seq = self.find_restriction_sites(original_section,
                                                                         lowest_site - max_restriction_site)
        restrictions_sites_in_mutated_seq = self.find_restriction_sites(mutated_section,
                                                                        lowest_site - max_restriction_site)

        if self.do_restriction_lists_dictionaries_match(restrictions_sites_in_original_seq,
                                                        restrictions_sites_in_mutated_seq):
            return None
        else:
            return self.find_extra_restriction_sites(restrictions_sites_in_original_seq,
                                                     restrictions_sites_in_mutated_seq)

    def get_restriction_sites(self, mutation_zone):
        max_len_site = max(list(map(lambda item: len(item.site), self.restriction_enzymes)))
        print("max len site:", max_len_site)

        start = max(0, mutation_zone.start-max_len_site)
        end = min(mutation_zone.end+1+max_len_site, len(self.mutated_strand)-1)
        mutated_section = self.mutated_strand[start:end]
        restrictions_sites_in_mutated_seq = self.find_restriction_sites(mutated_section,
                                                                        start)
        return restrictions_sites_in_mutated_seq

    def add_remove_restriction_sites(self,
                                     restriction_site_type):

        mutation_zone = self.get_mutations_zone(self.mutation_direction, self.pam_sites)
        print("Mutation zone:", mutation_zone)
        if restriction_site_type == RestrictionSiteType.INSERTED:
            # first check if the mutated strand so far (with codon mutations) has different number of restriction sites
            # than the original one
            print("Checking for already distinctive rest sites...")
            restriction_sites_dict = self.get_distinctive_restriction_sites()
            if restriction_sites_dict:  # value is not None
                print("dictionaries don't match!")  # different restriction sites
                # sorts and filters restriction sites
                rest_sites_list = self.check_mutated_restriction_sites(restriction_sites_dict)
                if rest_sites_list:
                    return rest_sites_list[0]

            # dictionaries match
            print("let's modify some sites! first, trying to add...")
            print("mutation zone:", mutation_zone)
            inserted_restriction_mutation = self.get_new_restriction_site(mutation_zone)
            if inserted_restriction_mutation:
                print("restriction mutation added!", inserted_restriction_mutation)
                self.mutated_strand = inserted_restriction_mutation.mutated_strand
                return inserted_restriction_mutation
            else:
                return None

        else:  # trying to remove a restriction site
            print("now let's try to remove...")
            return self.choose_restriction_site_to_remove(mutation_zone)

    # receives a dict of RestrictionSite as keys, and RestrictionSiteType as values of restriction sites that have been
    # mutated into insertion\deletion in prior steps. filters out irrelevant sites (too close to same sites) and sorts
    # them
    def check_mutated_restriction_sites(self, restriction_sites_dic: dict):
        rest_sites_list = list(restriction_sites_dic.keys())
        for rest_site in restriction_sites_dic:
            if self.check_distance(rest_site, 100) < 100:
                rest_sites_list.remove(rest_site)
        rest_sites_list.sort(key=lambda r_site: (r_site.enzyme in self.relevant_restriction_enzymes,
                                                 restriction_sites_dic[r_site] == RestrictionSiteType.INSERTED,
                                                 -(abs((r_site.index.start + r_site.index.end) / 2 - self.DSB)),
                                                 self.check_distance(r_site, vicinity=250),
                                                 - self.check_rareness(r_site, 250)),
                             reverse=True)
        print("Found restriction sits after sort:", *rest_sites_list, sep="\n")
        return rest_sites_list

    # receives the mutation zone and extracts all restriction sites there. after filtering and sorting, chooses a
    # restriction site to mutate and does so if possible
    def choose_restriction_site_to_remove(self, mutation_zone):
        rest_sites = self.get_restriction_sites(mutation_zone)
        if rest_sites:
            print("No restriction sits in the relevant section")
            return None
        print("Restriction Sites:", rest_sites)
        CrisprPlanner.filter_out_restriction_sites_with_no_space(rest_sites, self.mutated_strand)
        if rest_sites:
            # if there are still restriction sites left, we  will now sort them
            self.sort_restriction_sites(rest_sites)
            print("rest sites after sorting:", rest_sites)
            for rest_site in rest_sites:
                restriction_site_removed = self.remove_restriction_site(rest_site, mutation_zone,
                                                                        len(rest_sites))
                if restriction_site_removed:  # succeeded in finding and removing a restriction site
                    print("Restriction site removed! chosen restriction site:", rest_site)
                    return rest_site
            print("No restriction site to remove was found")
            return None
        else:
            print("No restriction sites left after filtering out")
            return None

    # receives (1)the chosen restriction site and (2) the start of the ssODN mutation codon tries to insert silent
    # mutations to remove it
    def remove_restriction_site(self, rest_site: RestrictionSite, mutation_zone, num_rest_sites):
        # it's not a reattachment section, but I need this namedTuple for the get_relevant_codons method
        section_to_mutate = ReattachmentSection(1, DNASection.MutationZone, rest_site.index)
        possible_mutations, _ = self.get_possible_codon_mutations(section_to_mutate)
        fourth_strand = self.mutated_strand
        fourth_sites = self.mutated_sites
        for index in range(len(possible_mutations)):
            self.mutated_strand = fourth_strand
            self.mutated_sites = fourth_sites
            point_mutations = self.get_point_mutations(possible_mutations, [index])  # also mutates the sequence
            if self.is_restriction_site_well_removed(point_mutations, rest_site, mutation_zone, num_rest_sites):
                self.mutated_sites.append(point_mutations)
                print("Restriction site removed with:", point_mutations)
                return True
        return False

    # receives (1) point mutations and the restriction site we wish to remove, mutate a copy of the relevant strand
    # according to the point mutations, take the part in which the restriction site is, and check if after the mutations
    # the restriction site is removed and can't function as a derivative of another restriction site
    def is_restriction_site_well_removed(self, point_mutations, rest_site: RestrictionSite, mutation_zone,
                                         num_rest_sites):
        for point_mutation in point_mutations:
            self.mutated_strand = self.change_char_in_string(self.mutated_strand,
                                                             point_mutation.index,
                                                             point_mutation.new_nucleotide)
            print("strand after restriction site removal mutation:", self.mutated_strand)
        start_index, end_index = rest_site.index.start, rest_site.index.end
        new_rest_site = self.mutated_strand[start_index:end_index + 1]
        if FileReader.is_derivative(rest_site.enzyme.site, new_rest_site):
            print("new rest site", new_rest_site, "is a derivative for", rest_site.enzyme.name, rest_site.enzyme.site)
            return False
        else:
            print("new rest site", new_rest_site, "is not a derivative for", rest_site.enzyme.name)
            # check if a new restriction site has been made
            new_rest_sites = self.find_sites(self.mutated_strand[mutation_zone.start:mutation_zone.end + 1])
            if len(new_rest_sites) == num_rest_sites:  # new site was created
                print("new rest site was created though")
                return False
            return True

    # as for the moment, not used
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
        rest_sites.sort(key=lambda r_site: (r_site in self.relevant_restriction_enzymes,
                                            -(abs((r_site.index.start + r_site.index.end)/2-self.DSB)),
                                            self.check_distance(r_site, vicinity=250),
                                            -self.check_rareness(r_site, 250)),
                        reverse=True)

    def check_rareness(self, restriction_site, distance=250):
        occurrences = []
        site_length = len(restriction_site.enzyme.site)
        mutated_section = self.mutated_strand[max(0, restriction_site.index.start - distance - site_length):
                                              min(len(self.mutated_strand) - 1,
                                                  restriction_site.index.end + distance + site_length)]
        for derivative in restriction_site.enzyme.derivatives:
            occurrences += [m.start() for m in re.finditer('(?=' + derivative + ')', mutated_section)]
            occurrences += [m.start() for m in re.finditer('(?=' + self.get_complementary_sequence(derivative) + ')',
                                                          mutated_section)]
        return len(occurrences)

    # receives a list of restriction sites and the strand that they are taken from, and looks for other same restriction
    # sites in the vicinity of certain nucleotides. if the vicinity is too close, remove the restriction site from the
    # list
    @staticmethod
    def filter_out_restriction_sites_with_no_space(rest_sites, mutated_strand, vicinity: int = 100):
        for restriction_site in rest_sites[:]:
            if CrisprPlanner.check_distance(restriction_site, mutated_strand, vicinity) < vicinity:
                rest_sites.remove(restriction_site)

    # receives a restriction site and a strand, and checks in what distance there is another restriction site from the
    # same type
    def check_distance(self, restriction_site: RestrictionSite, vicinity: int = 100):
        occurrences = []
        distance = vicinity + 1
        site_length = len(restriction_site.enzyme.site)
        start_point = max(0, restriction_site.index.start - vicinity - site_length)
        end_point = min(len(self.mutated_strand) - 1, restriction_site.index.end + vicinity + site_length)
        mutated_section = self.mutated_strand[start_point:end_point]
        print("current site being tested:", restriction_site)
        for derivative in restriction_site.enzyme.derivatives:
            occurrences += [m.start() for m in re.finditer('(?=' + derivative + ')', mutated_section)]
            occurrences += [m.start() for m in re.finditer('(?=' + self.get_complementary_sequence(derivative) + ')', mutated_section)]
            print("derivative", derivative, "complementary:", self.get_complementary_sequence(derivative))
        print("occurrences:", occurrences, ", start:", start_point)
        for occurrence in occurrences:
            if 0 < abs(start_point + occurrence - restriction_site.index.start) < distance:
                distance = abs(start_point + occurrence - restriction_site.index.start)
                print(distance, start_point+occurrence, mutated_section[occurrence:occurrence+9])
        return distance

    def get_mean_of_subtracted_lists(self, item: RestrictionMutation):
        def subtract_two_lists(lst1, lst2):
            return list(set(lst1) - set(lst2))
        only_new_mutated_sites = subtract_two_lists(item.mutated_sites, self.mutated_sites)
        if not only_new_mutated_sites:
            # if empty
            only_new_mutated_sites = self.mutated_sites
        return mean(list(map(lambda site: site.index, only_new_mutated_sites)))

    # gets a list of RestrictionMutation objects, sorts them and returns the optimal one
    def get_new_restriction_site(self, mutation_zone):
        possible_restriction_mutations = self.get_possible_restriction_mutations(mutation_zone)
        print("All possible restriction sites mutations:", *possible_restriction_mutations, sep="\n")

        if possible_restriction_mutations:
            possible_restriction_mutations.sort(key=lambda r_mutation: (r_mutation.number_of_mutations,
                                                          abs(self.DSB -
                                                              self.get_mean_of_subtracted_lists(r_mutation,
                                                                                                self.mutated_sites)),
                                                          -self.check_distance(r_mutation.restriction_site,
                                                                               vicinity=250),
                                                          self.check_rareness(r_mutation.restriction_site)))
            return possible_restriction_mutations[0]
        else:
            return None

    # try to replace get_possible_restriction_sites
    # aims to find an inserted restriction site and return a list of restriction mutation objects, each contains data to
    # achieve an inserted restriction site
    def get_possible_restriction_mutations(self, mutation_zone):
        possible_restriction_mutations = []
        # it's not a reattachment section, but I need this namedTuple for the get_relevant_codons method
        section_to_mutate = ReattachmentSection(1, DNASection.MUTATION_ZONE, mutation_zone)
        possible_mutations, _ = self.get_possible_codon_mutations(section_to_mutate)
        index_power_set = CrisprPlanner.get_power_set(list(range(len(possible_mutations))))

        original_strand = self.mutated_strand
        original_mutated_sites = self.mutated_sites
        print("The one and only original!!!!", original_mutated_sites)
        print("All index subsets to insert restriction site:", index_power_set)
        for index_subset in index_power_set:
            # check that mutation in subset don't override each other
            if CrisprPlanner.are_codons_overriden(index_subset, possible_mutations) and index_subset:
                print("subset", index_subset, "is ruled out because of overriden mutations")
            else:
                self.mutated_strand = original_strand
                self.mutated_sites = original_mutated_sites[:]
                print("original mutated sites:", original_mutated_sites)
                point_mutations = self.get_point_mutations(possible_mutations, index_subset)  # also mutates the sequence
                print("now working on index subset:", index_subset, "and point mutations:", point_mutations)
                self.mutated_sites.extend(point_mutations)
                print("all mutated sites:", self.mutated_sites)
                restriction_sites = self.get_distinctive_restriction_sites()
                if restriction_sites:
                    for restriction_site in restriction_sites:
                        if self.check_distance(restriction_site) > 100 and \
                                        restriction_sites[restriction_site] == RestrictionSiteType.INSERTED:
                            possible_restriction_mutations.append(RestrictionMutation(restriction_site,
                                                                                      len(index_subset),
                                                                                      self.mutated_sites))
        return possible_restriction_mutations

    # this function extracts all possible codons in the mutation zone section and tries to insert silent mutation in
    # them, recursively, and collects all options that add a new restriction site to the strand. number of max_mutations
    # is chosen
    # out of use function
    def get_possible_restriction_sites(self, mutated_strand, mutation_zone, mutated_sites, ssODN_mutation_codon,
                                       mutations_so_far: int = 0):
        possible_results = []
        section_to_mutate = ReattachmentSection(-1,
                                                DNASection.MUTATION_ZONE,
                                                SequenceSites(mutation_zone.start, mutation_zone.end))
        codons_data = CrisprPlanner.get_relevant_codons(section_to_mutate, ssODN_mutation_codon, mutated_strand)
        print("mutated sites so far:", mutated_sites)
        print("all possible codons for restriction mutations:", [data.codon for data in codons_data])
        for codon_data in codons_data:
            print("for codon:", codon_data.codon, "at", codon_data.codon_sites)
            # all codons_data that translate to the same amino acid as the relevant codon does
            same_aa_codons = list(CrisprPlanner.get_similar_codons(codon_data.codon))
            print("same codons:", same_aa_codons)
            CrisprPlanner.check_already_mutated_sites(codon_data, same_aa_codons, mutated_sites)
            # removes codons in which the difference is in part that in outside of the mutation zone
            CrisprPlanner.check_outside_codon_mutations(codon_data, same_aa_codons, section_to_mutate)

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
                    temporary_point_mutations.append(
                        PointMutation(codon_data.codon_sites.start + key, dict_of_changes[key]))
                    num_of_mutations += 1
                restriction_sites = self.get_distinctive_restriction_sites(mutated_sites + temporary_point_mutations,
                                                                           current_mutated_strand)
                if restriction_sites:
                    print("restriction sites found:", restriction_sites, codon_data, mutated_codon_details)
                    for restriction_site in restriction_sites:
                        if self.check_distance(restriction_site) > 100 and \
                                        restriction_sites[restriction_site] == RestrictionSiteType.INSERTED:
                            possible_results.append(RestrictionMutation(restriction_site,
                                                                        num_of_mutations,
                                                                        mutated_sites + temporary_point_mutations))
                print("going recursive!", num_of_mutations)
                possible_results.extend(self.get_possible_restriction_sites(current_mutated_strand, mutation_zone,
                                                                            mutated_sites + temporary_point_mutations,
                                                                            ssODN_mutation_codon,
                                                                            num_of_mutations))
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
            if restriction_site not in mutant_seq_restriction_sites:
                extra_restrictions_sites[restriction_site] = RestrictionSiteType.REMOVED

        # restriction site was added
        for restriction_site in mutant_seq_restriction_sites:
            if restriction_site not in original_seq_restriction_sites:
                extra_restrictions_sites[restriction_site] = RestrictionSiteType.INSERTED

        print("Found restriction sites:", *extra_restrictions_sites.items(), sep="\n")
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
        sites = []
        print("Finding sites in seq:", seq)
        for i in range(len(seq)):
            for j in range(i + 1, len(seq) + 1):
                for enzyme in self.restriction_enzymes:
                    for derivative in enzyme.derivatives:
                        if seq[i:j] == derivative:
                            print(enzyme.name, enzyme.site, derivative, i, j)
                            sites.append((enzyme.name, enzyme.site, derivative, i, j))
        return sites

    def get_relevant_restriction_enzymes(self, existing_enzymes_names: list):
        relevant_enzymes = list(filter(lambda r_enzyme: r_enzyme.name in existing_enzymes_names, self.restriction_enzymes))
        print("relevant enzymes:", relevant_enzymes)
        return relevant_enzymes

