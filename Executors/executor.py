from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType
from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileMaker import FileMaker
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.BioPython import BioPython
from Code.Utils.Ensembl import Ensembl
from Code.Utils.Strings import Strings
from Test.Extraction.DataExtracterTest import DataExtracterTest
from Test.TestFunctions import TestFunctions

import ast


class executor:

    # obtained the data of genes after filtering out genes by condition, reads files of human and C.elegans'
    # genes' coding domains length, and passes the info to a DE method, who computed the ratio between each gene
    # and its ortholog, and keeps only those who passes the bar. The compatible data is written to a new file.
    @staticmethod
    def filterGenesBySize(human_length_file_path, human_length_file_name: str, c_elegans_length_file_path: str,
                          c_elegans_length_file_name: str, new_file_name: str):
        det = DataExtracterTest()
        data = det.findGenesWithValidConditionAndHomologTest()
        print("obtained data...")

        fd_humans = FileReader(human_length_file_path, human_length_file_name, FileType.TSV)
        humanGenesLength = fd_humans.fromFileToDict(1, 2, True)
        print("number of items in human genes length dictionary is: " + str(len(humanGenesLength)))
        fd_c_elegans = FileReader(c_elegans_length_file_path, c_elegans_length_file_name, FileType.TSV)
        cElegansGenesLength = fd_c_elegans.get_genes_cd_length(0, 2, True)
        print("number of items in C.elegans genes length dictionary is: " + str(len(cElegansGenesLength)))

        de = DataExtracter()
        orthologs, phenotypes = de.filter_genes_with_size_differences(data, humanGenesLength, cElegansGenesLength, p=10)
        print("genes and orthologs dict: " + str(len(orthologs)))
        print("genes and phenotypes dict: " + str(len(phenotypes)))

        fm = FileMaker()
        fm.fromTwoDictToFile(orthologs, phenotypes, new_file_name)

    # creating a file to restore the c.elegans gene id (WB...) and the number of its conserved domains
    @staticmethod
    def c_elegans_conserved_domains_file():
        c_elegans_genes_id_WB = FileReader(FileReader.research_path + r"\Executors",
                                           r"\data-230619-fixed-ratio-with-C-elegans-phenotypes",
                                           FileType.TSV).get_c_elegans_genes_from_output(1)
        # c_elegans_genes_id_wb_number_dict = BioPython.get_genes_id(c_elegans_genes_id_WB)
        # FileMaker().fromDictToUniqueFile(c_elegans_genes_id_wb_number_dict, "relevant-c-elegans-genes-and-their-id-110619")

        de = DataExtracter()
        c_elegans_genes_and_conserved_domains = de.get_conserved_domains(c_elegans_genes_id_WB)
        print("The length of the genes_id_domain is " + str(len(c_elegans_genes_and_conserved_domains)))

        fm = FileMaker()
        fm.fromDictToFile(c_elegans_genes_and_conserved_domains, "c-elegans-genes-and-conserved-domains-230619")

    # creating a file to restore the c.elegans gene id (WB...) and the number of its conserved domains
    @staticmethod
    def human_conserved_domains_file():
        fr2 = FileReader(FileReader.research_path + r"\Executors",
                         r"\genes-orthologs-phenotypes-filteredBySize-100619",
                         FileType.TSV)
        genes_names = fr2.get_genes_list()

        de = DataExtracter()
        genes_id = de.from_human_genes_names_to_human_genes_id(genes_names)
        print("The length of the genes ids list is " + str(len(genes_id)))
        print("human genes ids for example: " + genes_id[0] + ", " + genes_id[1])

        human_genes_and_conserved_domains = de.get_conserved_domains(genes_id)

        fm = FileMaker()
        fm.fromDictToFile(human_genes_and_conserved_domains, "human-genes-and-conserved-domains-110619")

    @staticmethod
    def filter_by_conserved_domains_ratio(orthologs_dic, domains_range):
        # first we need the C.elegans genes id-number of domains dictionary
        c_elegans_genes_domains_dic = FileReader(FileReader.research_path + r"\Data",
                                                 r"\c-elegans-genes-and-conserved-domains-230619",
                                                 FileType.TSV).fromFileToDict(0, 1)
        # second we need the human genes id-number of domains dictionary
        human_genes_domains_dic = FileReader(FileReader.research_path + r"\Data",
                                             r"\human-genes-and-conserved-domains-230619",
                                             FileType.TSV).fromFileToDict(0, 1)
        relevant_orthologs_dic = DataExtracter.filter_by_conserved_domains(c_elegans_genes_domains_dic,
                                                                           human_genes_domains_dic,
                                                                           orthologs_dic,
                                                                           domains_range)
        return relevant_orthologs_dic

    # using the (1) C.elegans genes conserved domains file, (2) human genes conserved domains file, (3) most
    # recent data file, and (4) human gene id-gene name file, extracts the info needed, conducts the calculation
    # needed and creates a file of data that includes the ratio of conserved domains
    @staticmethod
    def addDomainsScoreInfoToFile():
        # first we need the C.elegans genes-number of domains dictionary
        fr1 = FileReader(FileReader.research_path + r"\Executors",
                         r"\c-elegans-genes-and-conserved-domains-110619",
                         FileType.TSV)
        c_elegans_genes_domains_dic = fr1.fromFileToDict(0, 1)
        tf1 = TestFunctions("fromFileToDict to read c.elegans genes and domains",
                            dictionary=c_elegans_genes_domains_dic)
        tf1.checkSize()
        tf1.printFirstLinesInDict(2)

        # second we need the human genes-number of domains dictionary
        fr2 = FileReader(FileReader.research_path + r"\Executors",
                         r"\human-genes-and-conserved-domains-110619",
                         FileType.TSV)
        human_genes_domains_dic = fr2.fromFileToDict(0, 1)
        tf2 = TestFunctions("fromFileToDict to read human genes and domains",
                            dictionary=human_genes_domains_dic)
        tf2.checkSize()
        tf2.printFirstLinesInDict(2)

        # then we'll need the dictionary to convert human gene name to its c.elegans orthologs
        fr3 = FileReader(FileReader.research_path + r"\Executors",
                         r"\genes-orthologs-phenotypes-filteredBySize-090619",
                         FileType.TSV)
        orthologs = fr3.fromFileToDict(0, 1)
        tf3 = TestFunctions("fromFileToDict to read human genes and orthologs",
                            dictionary=orthologs)
        tf3.checkSize()
        tf3.printFirstLinesInDict(2)

        # and last we need a dictionary to convert human gene name to human gene id
        fr4 = FileReader(FileReader.research_path + r"\Data",
                         r"\human_gene_id-gene_name.txt",
                         FileType.TSV)
        human_genes_names_ids_dict = fr4.fromFileToDict(1, 0, True)
        tf4 = TestFunctions("fromFileToDict to read human genes names and ids",
                            dictionary=human_genes_names_ids_dict)
        tf4.checkSize()
        tf4.printFirstLinesInDict(2)

        human_genes_id_orthologs_dic = DataExtracter.add_conserved_domains_info(c_elegans_genes_domains_dic,
                                                                                human_genes_domains_dic,
                                                                                orthologs,
                                                                                human_genes_names_ids_dict)
        tf5 = TestFunctions("add_conserved_domains_info",
                            dictionary=human_genes_id_orthologs_dic)
        tf5.checkSize()
        tf5.printFirstLinesInDict(7)

        # now lets write all this back to a file, but for that we need two things:
        # first, a converter from gene id to gene name:
        human_genes_ids_names_dict = fr4.fromFileToDict(0, 1, True)
        tf6 = TestFunctions("fromFileToDict to read human genes ids and names",
                            dictionary=human_genes_ids_names_dict)
        tf6.checkSize()
        tf6.printFirstLinesInDict(2)

        # and second, a dictionary of human genes' names and conditions
        fr5 = FileReader(FileReader.research_path + r"\Executors",
                         r"\genes-orthologs-phenotypes-filteredBySize-100619",
                         FileType.TSV)
        genes_names_conditions = fr5.fromFileToDict(0, 2)
        tf7 = TestFunctions("fromFileToDict for dictionary of human genes' names and conditions",
                            dictionary=genes_names_conditions)
        tf7.checkSize()
        tf7.printFirstLinesInDict(2)

        fm = FileMaker()
        fm.fromTwoDictOfDifferentKeysToOneFile(human_genes_id_orthologs_dic,
                                               genes_names_conditions,
                                               human_genes_ids_names_dict,
                                               "genes_id-genes_names-orthologs-conditions-110619")

    @staticmethod
    def getOnlyGenesWithHumanOrtholog():

        print("*** Function \"getOnlyGenesWithHumanOrtholog\" has started ***")

        accession_numbers_and_hit_ids, hit_ids_and_hsps = DataExtracter.get_hit_ids_for_accession_numbers(
            31,
            FileReader.research_path + r"\Test\Files\accessions-and-sequences",
            r"\results-part-")
        print("accessions numbers and hit ids and hsp have been obtained")

        # hit_ids_and_accessions = DataExtracter.from_plural_valued_dict_to_plural_valued_reversed_dict(accession_numbers_and_hit_ids)
        # hit_ids = hit_ids_and_accessions.keys()
        # print("There are " + str(len(hit_ids_and_accessions)) + " hit ids")
        # hit_ids_and_gene_names = BioPython.get_gene_name_from_protein_accession(hit_ids)
        #
        # print("Gene names for the hit ids have been obtained!")
        # tf1 = TestFunctions("get_gene_name_from_protein_accession",
        #                     dictionary=hit_ids_and_gene_names)
        # tf1.checkSize()
        # tf1.printFirstLinesInDict(2)

        # FileMaker().fromDictToFile(hit_ids_and_gene_names, "blast_hit_ids_and-c-elegans-gene_names0-30")
        # FileMaker().fromDictToFile(hit_ids_and_hsps, "blast_hit_ids_and_hsp_scores0-30")

        hit_ids_and_gene_names = FileReader(FileReader.research_path + r"\Executors",
                                            r"\blast-hit-ids-and-c-elegans-gene-names0-30",
                                            FileType.TSV).fromFileToDict(0, 1)

        tf1 = TestFunctions("get_gene_name_from_protein_accession",
                            dictionary=hit_ids_and_gene_names)
        tf1.checkSize()
        tf1.printFirstLinesInDict(2)

        gene_ids_and_accession_numbers = FileReader(FileReader.research_path + r"\Executors",
                                                    r"\extra_genes_ids_number_and_chosen_accessions-120619",
                                                    FileType.TSV).fromFileToDict(0, 1, False)
        tf2 = TestFunctions("fromFileToDict to get genes id and the accession numbers",
                            dictionary=gene_ids_and_accession_numbers)
        tf2.checkSize()
        tf2.printFirstLinesInDict(2)

        gene_names_WB_and_gene_ids = FileReader(FileReader.research_path + r"\Executors",
                                                r"\fixed-relevant-c-elegans-genes-and-their-id-110619",
                                                FileType.TSV).fromFileToDict(0, 1)

        tf3 = TestFunctions("fromFileToDict to get genes names and the gene ids",
                            dictionary=gene_names_WB_and_gene_ids)
        tf3.checkSize()
        tf3.printFirstLinesInDict(2)

        # last one - gene names to their homoloug

        gene_ids_WB_and_human_names = FileReader(FileReader.research_path + r"\Executors",
                                                 r"\genes_id-genes_names-orthologs-conditions-110619",
                                                 FileType.TSV).getCelegansIdToHumanNameDict()

        tf4 = TestFunctions("getCelegansIdToHumanNameDict",
                            dictionary=gene_ids_WB_and_human_names)
        tf4.checkSize()
        tf4.printFirstLinesInDict(2)

        gene_ids_WB_and_true_homolougs = {}

        WBGenes = gene_ids_WB_and_human_names.keys()
        print("There are " + str(len(WBGenes)) + " WBGenes")
        for WBGene in WBGenes:
            gene_name = gene_ids_WB_and_human_names[WBGene]
            print("human gene name for " + WBGene + ": " + gene_name)
            try:
                gene_id = gene_names_WB_and_gene_ids[WBGene]
                print("c-elegans gene id for " + gene_name + ": " + gene_id)
            except:
                print("couldn't find the gene id for " + gene_name)
                continue
            try:
                accession_number = gene_ids_and_accession_numbers[gene_id]
            except:
                print("couldn't find the accession number for " + gene_id)
                continue
            print("c-elegans accession number for " + gene_id + ": " + accession_number)
            try:
                hit_ids = accession_numbers_and_hit_ids[accession_number]
                print("hit ids for " + accession_number + ": " + str(hit_ids))
            except KeyError:
                print("couldn't find hit ids for " + accession_number)
                continue

            for hit_id in hit_ids:
                try:
                    blasted_gene_name = hit_ids_and_gene_names[hit_id]
                except KeyError:
                    print("couldn't find " + hit_id + " in the hit_ids_and_gene_names_dict")
                    continue
                print(blasted_gene_name, gene_name)
                if blasted_gene_name == gene_name:
                    gene_ids_WB_and_true_homolougs[WBGene] = gene_name

        tf5 = TestFunctions("unplaced function to find true homologs",
                            dictionary=gene_ids_WB_and_true_homolougs)
        tf5.checkSize()
        tf5.printFirstLinesInDict(2)

        FileMaker().fromDictToFile(gene_ids_WB_and_true_homolougs, "true-homologs0-30")

        human_genes_and_c_elegans_WB_id = FileReader(FileReader.research_path + r"\Executors",
                                                     r"\true-homologs0-30",
                                                     FileType.TSV).fromFileToDictWithPluralValues(1, 0)
        print("There are " + str(len(human_genes_and_c_elegans_WB_id)) + " human genes with homologs")

    # uses the true homologs files created by the former function to filter the list of genes
    @staticmethod
    def filterGenesAccordingToReversedBlast(results_file_names: list, data_file_name: str, new_data_file_name: str):
        # first we make a dictionary of human genes as keys and their true C.elegans orthologs as list of values
        human_genes_and_c_elegans_orthologs = {}
        for file_name in results_file_names:
            f = open(file_name)
            for line in f:
                genes = line.rstrip("\n").split("\t")
                c_elegans_gene = genes[0]
                human_gene = genes[1]
                DataExtracter.add_to_dictionary(human_genes_and_c_elegans_orthologs, human_gene, c_elegans_gene)
            f.close()
        FileMaker().fromPluralValuedDictToFile(human_genes_and_c_elegans_orthologs, "human-and-c-elegans-genes")

        # now we filter the data file
        data_file = open(data_file_name)
        new_data_file = open(new_data_file_name, FileMode.WRITE.value)
        for row in data_file:
            orthologs_to_be_written = []
            line = row.rstrip("\n").split("\t")
            human_gene_id = line[0]
            human_gene_name = line[1]
            c_elegans_orthologs = FileReader.fromStringToTuplesList(line[2])
            phenotype = line[3]
            if human_gene_name in human_genes_and_c_elegans_orthologs:
                true_orthologs = human_genes_and_c_elegans_orthologs[human_gene_name]
                for gene_tuple in c_elegans_orthologs:
                    maybe_ortholog = gene_tuple[0]  # the c_elegans_name
                    if maybe_ortholog in true_orthologs:
                        orthologs_to_be_written.append(gene_tuple)
            else:
                print("human gene name " + human_gene_name + " wasn't found in the dictionary of true orthologs :(")
                print("meaning no one of his claimed orthologs was found to be real... check it")
                continue
            new_data_file.write(human_gene_id + "\t" +
                                human_gene_name + "\t" +
                                str(orthologs_to_be_written).strip('[]') + "\t" +
                                phenotype + "\n")
        data_file.close()
        new_data_file.close()

    @staticmethod
    def get_human_genes_and_true_orthologs():
        human_genes_and_c_elegans_WB_id = FileReader(FileReader.research_path + r"\Executors",
                                                     r"\true-homologs0-18",
                                                     FileType.TSV).fromFileToDictWithPluralValues(1, 0)
        print("There are " + str(len(human_genes_and_c_elegans_WB_id)) + " human genes with homologs")
        FileMaker().fromPluralValuedDictToFile(human_genes_and_c_elegans_WB_id, "human-genes-and-true-orthologs")

    @staticmethod
    def checkExtraGenes():
        # obtaining the new C.elegans genes id (number) list
        c_elegans_current_genes_id_number = FileReader(FileReader.research_path + r"\Executors",
                                                       r"\fixed-relevant-c-elegans-genes-and-their-id-110619",
                                                       FileType.TSV).get_genes_list(1)

        # obtaining the former C.elegans genes id (number) list
        c_elegans_former_genes_id_number = FileReader(FileReader.research_path + r"\Test\Files",
                                                      r"\relevant-c-elegans-genes-and-their-id",
                                                      FileType.TSV).get_genes_list(1)

        print("length of now-genes is: " + str(len(c_elegans_current_genes_id_number)) + " and length of former genes "
                                                                                         "is: " + str(
            len(c_elegans_former_genes_id_number)))
        extra_genes = []
        for gene in c_elegans_current_genes_id_number:
            if gene not in c_elegans_former_genes_id_number:
                print(gene + " does not exist in the former list")
                extra_genes.append(gene.rstrip("\n"))
        print("There are " + str(len(extra_genes)) + " extra genes ids_number that needs a blast")

        fm = FileMaker()
        fm.fromListToFile(extra_genes, "extra-c-elegans-gened-id-number-120619")

    @staticmethod
    def fromAccessionsToBLAST():
        # extra_genes_id_number_and_accessions_dict = FileReader(FileReader.research_path + r"\Data",
        #                                                        r"\extra-gene-ids-number-accession-numbers.txt",
        #                                                        FileType.TSV).fromFileToDictWithPluralValues(0,1)
        # print("Size of genes and accessions is " + str(len(extra_genes_id_number_and_accessions_dict)))
        # extra_genes_id_number_and_chosen_accessions = DataExtracter.from_multiple_accessions_to_one(
        #     extra_genes_id_number_and_accessions_dict)
        # print("Size of genes and chosen accessions is " + str(len(extra_genes_id_number_and_chosen_accessions)))
        # FileMaker().fromDictToFile(extra_genes_id_number_and_chosen_accessions, "extra_genes_ids_number_and_
        # chosen_accessions-120619")

        # accessions = FileReader(FileReader.research_path + r"\Executors",
        #                         r"\extra_genes_ids_number_and_chosen_accessions-120619",
        #                         FileType.TSV).get_genes_list(1)
        # print("Got " + str(len(accessions)) + " accessions! let the work begin...")
        #
        # bp = BioPython()
        # accessionsAndSeqs = bp.make_accession_number_and_seq_dict(accessions, "ORIGIN", "translation=")
        # FileMaker().fromDictToFile(accessionsAndSeqs, "extra_chosen_accessions_and_sequences")

        accessionsAndSeqs = FileReader(FileReader.research_path + r"\Test\Files\accessions-and-sequences",
                                       r"\accessions-and-sequences-part-helper",
                                       FileType.TSV).fromFileToDict(0, 1)
        BioPython().blastp_by_accessions("blastp", "nr", accessionsAndSeqs)

    # one time function
    @staticmethod
    def fixConservedDomainRatioFile(data_file_path, data_file_name, delete_first_line):
        # first we need the c-elegans genes-number of domains dictionary
        fr1 = FileReader(FileReader.research_path + r"\Executors",
                         r"\c-elegans-genes-and-conserved-domains-110619",
                         FileType.TSV)
        c_elegans_genes_domains_dic = fr1.fromFileToDict(0, 1)
        tf1 = TestFunctions("fromFileToDict to read c.elegans genes and domains",
                            dictionary=c_elegans_genes_domains_dic)
        tf1.checkSize()
        tf1.printFirstLinesInDict(2)

        # second we need the human genes-number of domains dictionary
        fr2 = FileReader(FileReader.research_path + r"\Executors",
                         r"\human-genes-and-conserved-domains-110619",
                         FileType.TSV)
        human_genes_domains_dic = fr2.fromFileToDict(0, 1)
        tf2 = TestFunctions("fromFileToDict to read human genes and domains",
                            dictionary=human_genes_domains_dic)
        tf2.checkSize()
        tf2.printFirstLinesInDict(2)

        # now we change the data
        DataExtracter.fix_conserved_domain_info(data_file_path, data_file_name, c_elegans_genes_domains_dic,
                                                human_genes_domains_dic, 4, "data-190619-short-fixed-domains",
                                                delete_first_line)

    # receives (1) data file path, (2) data file name, (3) url address, (4) new file name, (5) boolean value to indicate
    # whether or not we need to disregard the first line of the data, it retrieves the worm genes names from the data,
    # pass it to function that returns a dictionary of said genes and phenotypes, and copies the data with phenotype
    # to each worm gene to the new file
    @staticmethod
    def add_c_elegans_phenotypes(data_file_path, data_file_name, url, new_file, delete_first_line: bool = False):
        # first we achieve a list of all relevant C.elegans genes
        cElegansGenes = FileReader(data_file_path, data_file_name, FileType.TSV).get_genes_list(2)
        print("List of genes with", str(len(cElegansGenes)), "have been obtained:", cElegansGenes[0], ",",
              cElegansGenes[1])

        genesAndPhenotypes = DataExtracter().add_c_elegans_phenotypes(cElegansGenes, url)

        in_file = open(data_file_path + data_file_name)
        out_file = open(new_file, FileMode.WRITE.value)
        if delete_first_line:
            in_file.readline()
        for line in in_file:
            info = line.rstrip("\n").split("\t")
            cElegansGene = info[2]
            if cElegansGene in genesAndPhenotypes:
                phenotypes: list = genesAndPhenotypes[cElegansGene]
            else:
                phenotypes = ["Not Found"]
            out_file.write(line.rstrip("\n") + "\t" + ", ".join(phenotypes) + "\n")
        in_file.close()
        out_file.close()

    # receives (1) a dictionary of human-worm pairs, if by file or if by a human-worm-id dictionary and for each pair,
    # and the species whose genes we wish to blast, and returns whether the reversed blast confirms this pair as
    # orthologous.
    @staticmethod
    def pair_pipeline(list_of_pairs_file_path="", list_of_pair_file_name="", dic_of_optional_orthologs: dict = None,
                      key_species="C.elegans"):
        if not dic_of_optional_orthologs:
            pairs_in_names = FileReader(list_of_pairs_file_path, list_of_pair_file_name,
                                        FileType.TSV).fromFileToDictWithPluralValues(0, 1, True)
        else:
            orthologs_key_pairs_id = DataExtracter().from_plural_valued_dict_to_plural_valued_reversed_dict(dic_of_optional_orthologs)
            print("reversed dic:", orthologs_key_pairs_id)
            pairs_in_names = DataExtracter().convert_dic(orthologs_key_pairs_id, True)
            print("pairs in names:", pairs_in_names)

        print("List of all pairs have been obtained! we have " + str(len(pairs_in_names)) + " pairs to check")

        # accession_numbers_and_hit_ids, hit_ids_and_hsps = DataExtracter.get_hit_ids_for_accession_numbers(
        #     31,
        #     FileReader.research_path + r"\Test\Files",
        #     r"\blast-results-complete-200619")
        #
        # FileMaker().fromPluralValuedDictToFile(accession_numbers_and_hit_ids, "accession numbers and hit ids")
        # FileMaker().fromPluralValuedDictToFile(hit_ids_and_hsps, "hit ids and hsp")
        # tf = TestFunctions("get hit ids for accession numbers", dictionary=accession_numbers_and_hit_ids)
        # tf.checkSize()
        # tf.printRandomLinesInDict(5)

        true_matches, false_matches = {}, {}
        de = DataExtracter()
        for key_gene_name in pairs_in_names:
            ortholog_genes_names: list = pairs_in_names[key_gene_name]
            for ortholog_gene_name in ortholog_genes_names:
                c_elegans_gene_name = key_gene_name if key_species == "C.elegans" else ortholog_gene_name
                human_gene_name = ortholog_gene_name if key_species == "C.elegans" else key_gene_name
                print("Now working on " + key_gene_name + " and " + ortholog_gene_name)
                key_gene_id = Ensembl.get_gene_id_by_gene_name(key_gene_name, key_species)
                result = de.check_reversed_blast_hit_ids(key_gene_name, key_gene_id, ortholog_gene_name, key_species)
                genes_tuple = (ortholog_gene_name, key_gene_name)
                status_tuple = (de.get_sources(c_elegans_gene_name, human_gene_name),
                                de.get_conserved_domains_ratio_of_pair(c_elegans_gene_name, human_gene_name),
                                de.get_pair_cd_length(c_elegans_gene_name, human_gene_name),
                                c_elegans_gene_name + " " + de.get_c_elegans_description_for_gene_id(
                                    Ensembl.get_gene_id_by_gene_name(c_elegans_gene_name, "C.elegans")))
                if result:
                    true_matches[genes_tuple] = status_tuple
                else:
                    false_matches[genes_tuple] = status_tuple
                print("The domains ratio, number of sources, human gene length and C.elegans gene length: ",
                      str(status_tuple))

        print("true matches:", str(len(true_matches)), "false matches:", str(len(false_matches)))
        return true_matches

    @staticmethod
    def filterGenesTest():
        DataExtracter.filter_genes_by_name(FileReader.research_path + r"\Executors",
                                           r"\data-250619-fixedDomains-short-lethalFiltered",
                                           FileReader.research_path + r"\Data",
                                           r"\data-010719-sorted",
                                           "data-010719-filteredGenesByUnwantedGenes")

    # receives (1) file type to know whether to print to a file or to console, read all variants and extracts for each
    # variant the sequence of its human gene and its C.elegans ortholog, runs pairwise alignment and returns data
    # regarding the alignments of the two sequences. also filter out all genes for which the amino acid mutated is not
    # conserved
    # not used method, outdated 150919
    @staticmethod
    def get_variants_dict(file_type=FileType.CONSOLE):
        human_genes_and_sequences = {}
        human_genes_and_variants = FileReader(FileReader.research_path + r"\Data\variants",
                                              r"\all-missense-nonsense").fromHGMDtoDict()
        data = FileReader(FileReader.research_path + r"\Data",
                          r"\data-250619-fixed-domains").readData(1, False)
        genes_names = human_genes_and_variants.keys()

        orthologs_dic = FileReader(FileReader.research_path + r"\Data",
                                   r"\data-250619-fixed-domains").fromFileToDictWithPluralValues(1, 2, False)
        c_elegans_id_and_accessions = FileReader(
            FileReader.research_path + r"\Test\Files",
            r"\c-elegans-genes-and-longest-accession_number_complete").fromFileToDict(0, 1, False)
        # new file: FileReader.research_path + r"\Test\Files\c-elegans-gene-ids-and-accession-numbers"

        accessions_and_sequences = FileReader(FileReader.research_path + r"\Test\Files",
                                              r"\all-accession-numbers-and-sequences").fromFileToDict(0, 1)
        # new file: FileReader.research_path + r"\Extraction\c-elegans-accession-numbers-and-sequences""

        mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                           r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(9, [11, 18, 12], delete_first=True)

        f = open("data-with-variants-230919", FileMode.WRITE.value) if file_type != FileType.CONSOLE else ""
        print(human_genes_and_variants)
        print("Number of genes: " + str(len(human_genes_and_variants)))
        for human_gene_name in genes_names:
            human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
            human_seq = HttpRequester().get_human_protein_sequence_from_uniprot(human_gene_id)
            human_genes_and_sequences[human_gene_name] = human_seq
            mmp_data = mmp_data_by_gene_name[human_gene_name] if human_gene_name in mmp_data_by_gene_name \
                else "Not mention in MMP"
            print("Human gene name: " + human_gene_name + ", seq: " + human_seq)
            orthologs_id_WB = orthologs_dic[human_gene_name]
            print("Ortholog gene id WB: ", str(orthologs_id_WB))
            for ortholog_id_WB in orthologs_id_WB:
                ortholog_id_number = Ensembl.get_ncbi_id_by_gene_id(ortholog_id_WB)
                if not ortholog_id_number:
                    ortholog_id_number = BioPython.get_gene_id(ortholog_id_WB)
                    if not ortholog_id_number:
                        print("Couldn't find C.elegans' id (number), moving on to the next gene")
                        continue
                print("id: " + ortholog_id_number)
                try:
                    c_elegans_accession_number = c_elegans_id_and_accessions[ortholog_id_number]
                except:
                    print("Couldn't find C.elegans' accession number, moving on to the next gene")
                    # right now no way to get accession number other that DAVID
                    continue
                print("accession: " + c_elegans_accession_number)
                try:
                    c_elegans_seq = accessions_and_sequences[c_elegans_accession_number]
                except:
                    c_elegans_seq = BioPython.get_aa_seq_from_entrez(c_elegans_accession_number)
                print("C.elegans seq: " + c_elegans_seq)

                variants = human_genes_and_variants[human_gene_name]
                print("variants: " + ", ".join(variants))
                for variant in variants:
                    print("For variant: " + variant)
                    former_aa, place, current_aa = Strings.fromVariantStringToTuple(variant)
                    result, count, c_elegans_location, alignment_conservation_score = \
                        BioPython.pairwise_alignment_inspector(human_seq,
                                                               c_elegans_seq,
                                                               Strings.fromNameToSymbol(former_aa),
                                                               place)
                    if result.startswith("not conserved") or result.startswith("not in the human"):
                        # filter if amino acid is not conserved
                        continue
                    print("For gene: " + human_gene_name + " with variant " + variant + ", the amino acid is " + result)
                    line = data[human_gene_name] + "\t" + variant + "\t" + result + "\t" + str(
                        c_elegans_location) + "\t" + str(alignment_conservation_score) + "\t" + str(count) + "\t" + \
                           mmp_data + "\n"
                    FileMaker().write_to(file_type, line, f)
        if file_type != FileType.CONSOLE:
            f.close()

    # receives (1) file type to know whether to print to a file or to console, (2) human genes and  variants dictionary,
    # read all variants and extracts for each variant the sequence of its human gene and its C.elegans ortholog,
    # runs pairwise alignment and returns data regarding the alignments of the two sequences. also filter out all genes
    # for which the amino acid mutated is not conserved
    @staticmethod
    def get_variants_data(file_type,
                          human_genes_and_variants):
        conserved_variants = 0

        genes_names = human_genes_and_variants.keys()

        mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                           r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(9, [11, 18, 12], delete_first=True)

        f = open("data-with-variants", FileMode.WRITE.value) if file_type != FileType.CONSOLE else ""
        print(human_genes_and_variants)
        print("Number of human genes: " + str(len(human_genes_and_variants)))
        for human_gene_name in genes_names:
            human_seq = BioPython().get_aa_seq_by_human_gene_name(human_gene_name)
            if not human_seq:
                continue
            mmp_data = mmp_data_by_gene_name[human_gene_name] if human_gene_name in mmp_data_by_gene_name \
                else "No mention in MMP"
            try:
                true_matches_pairs = list(executor.find_me_orthologs([human_gene_name]).keys())
                orthologs_names = []
                for i in range(len(true_matches_pairs)):
                    pair = true_matches_pairs[i]
                    if pair[1] == human_gene_name:
                        orthologs_names.append(pair[0])
            except:
                orthologs_names_string: str = input("What are the C.elegans orthologs for the human gene {}? Please "
                                                    "write the genes in a list format".format(human_gene_name))
                if orthologs_names_string.startswith("[") and orthologs_names_string.endswith("]"):
                    orthologs_names = ast.literal_eval(orthologs_names_string)
                else:
                    print("Orthologs are not valid. They need to be given as a list, try again")
                    continue
            orthologs_id_WB = Ensembl.get_c_elegans_genes_ids_by_genes_names(orthologs_names)
            print("Human gene name:", human_gene_name, ", seq:", human_seq, ", orthologs gene ids WB:", orthologs_id_WB)

            for ortholog_id_WB in orthologs_id_WB:
                if not ortholog_id_WB:
                    print("Ortholog id WB couldn't be reached")
                    continue
                print("for ortholog:", ortholog_id_WB)
                c_elegans_seq = BioPython().get_c_elegans_aa_seq(ortholog_id_WB)
                if not c_elegans_seq:
                    continue
                print("C.elegans seq:", c_elegans_seq)

                variants = human_genes_and_variants[human_gene_name]
                print("variants: " + ", ".join(variants))
                for variant in variants:
                    print("For variant: " + variant)
                    former_aa, place, current_aa = Strings.fromVariantStringToTuple(variant)
                    result, count, c_elegans_location, alignment_conservation_score = \
                        BioPython.pairwise_alignment_inspector(human_seq,
                                                               c_elegans_seq,
                                                               Strings.fromNameToSymbol(former_aa),
                                                               place)
                    print("For gene: " + human_gene_name + " with variant " + variant + ", the amino acid is " + result)
                    if not result.startswith("conserved"):
                        # filter if amino acid is not conserved
                        continue
                    conserved_variants += 1
                    line = variant + "\t" + result + "\t" + str(c_elegans_location) + "\t" + \
                        str(alignment_conservation_score) + "\t" + str(count) + "\t" + mmp_data + "\n"
                    FileMaker().write_to(file_type, line, f)
        if not conserved_variants:
            print("No variants are conserved")
        if file_type != FileType.CONSOLE:
            f.close()

    # receives (1) file type to know whether to print to a file or to console, human genes and variants dictionary, read
    # all variants and extracts for each variant the sequence of its human gene and its C.elegans ortholog, runs
    # pairwise alignment and returns data regarding the alignments of the two sequences.
    @staticmethod
    def get_variants_data_for_server(human_genes_and_variants):
        conserved_variants = 0
        output = ''
        genes_names = human_genes_and_variants.keys()

        mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                           r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(9, [11, 18, 12], delete_first=True)

        print("Human genes: " + str(list(genes_names)))
        for human_gene_name in genes_names:
            human_seq = BioPython().get_aa_seq_by_human_gene_name(human_gene_name)
            if not human_seq:
                continue
            mmp_data = mmp_data_by_gene_name[human_gene_name] if human_gene_name in mmp_data_by_gene_name \
                else "No mention in MMP"

            true_matches_pairs = list(executor.find_me_orthologs([human_gene_name]).keys())
            orthologs_names = []
            for i in range(len(true_matches_pairs)):
                pair = true_matches_pairs[i]
                orthologs_names.append(pair[0])

            orthologs_id_WB = Ensembl.get_c_elegans_genes_ids_by_genes_names(orthologs_names)
            print("Human gene name:", human_gene_name, ", seq:", human_seq, ", orthologs gene ids WB:",
                   orthologs_id_WB)

            for ortholog_id_WB in orthologs_id_WB:
                if not ortholog_id_WB:
                    print("Ortholog id WB couldn't be reached")
                    continue
                print("for ortholog:", ortholog_id_WB)
                c_elegans_seq = BioPython().get_c_elegans_aa_seq(ortholog_id_WB)
                if not c_elegans_seq:
                    continue
                print("C.elegans seq:", c_elegans_seq)

                variants = human_genes_and_variants[human_gene_name]
                print("variants: " + ", ".join(variants))
                for variant in variants:
                    print("For variant: " + variant)
                    former_aa, place, current_aa = Strings.fromVariantStringToTuple(variant)
                    result, count, c_elegans_location, alignment_conservation_score = \
                        BioPython.pairwise_alignment_inspector(human_seq,
                                                               c_elegans_seq,
                                                               Strings.fromNameToSymbol(former_aa),
                                                               place)
                    print("For gene: " + human_gene_name + " with variant " + variant + ", the amino acid is " + result)
                    conserved_variants += 1
                    line = variant + "\t" + result + "\t" + str(c_elegans_location) + "\t" + \
                        str(alignment_conservation_score) + "\t" + str(count) + "\t" + mmp_data + "\n"
                    output += line
            if not conserved_variants:
                print("No variants are conserved")
        return output

    # keys = human gene names
    @staticmethod
    def add_mmp_record_to_data(data_file_path, data_file_name, key_index, value_indexes, new_file_name):
        data = FileReader(data_file_path, data_file_name).readData(1, True)
        mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                           r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(key_index, value_indexes, delete_first=True)
        f = open(new_file_name, FileMode.WRITE.value)
        gene_names = data.keys()
        for gene_name in gene_names:
            if gene_name in mmp_data_by_gene_name:
                f.write(data[gene_name] + "\t" + "".join(mmp_data_by_gene_name[gene_name]) + "\n")
            else:
                f.write(data[gene_name] + "\t" + "Not mentioned in MMP" + "\n")
        f.close()

    @staticmethod
    def find_me_orthologs(list_of_human_genes,
                          genes_in_names: bool = True,
                          sources_bar: int = 3,
                          length_bar: int = 10,
                          domains_range: tuple = (0.5, 2)):
        if genes_in_names:
            list_of_human_genes_names = list_of_human_genes
            list_of_human_genes = DataExtracter().convert_list_for_human(list_of_human_genes_names)
            print("Genes Ids:" + ", ".join(list_of_human_genes))

        # from human gene id to C.elegans gene id dictionary
        orthologs_dic = FileReader(FileReader.research_path + r"\Data",
                                   r"\ortholist_master",
                                   FileType.TSV).fromFileToDictWithPluralValues(4, 0, True)
        relevant_orthologs_dic = DataExtracter.get_specific_dic_of_orthologs(list_of_human_genes, orthologs_dic)
        DataExtracter.check_if_dict_not_empty(relevant_orthologs_dic)
        print("Orthologs: " + str(relevant_orthologs_dic))

        # now we have the ortholist-orthologs to our human genes.
        # next step - filtration by sources
        # sources dic is a dictionary with a tuple-keys of (human gene id, c.elegans gene id) and values of number of
        # sources supporting the pair is homologous
        sources_dic = FileReader(FileReader.research_path + r"\Data",
                                 r"\ortholist_master",
                                 FileType.TSV).fromFileToDictWithTupleKey(4, 0, 6, True)
        filtered_by_sources_orthologs = DataExtracter.filter_dic_by_sources(relevant_orthologs_dic,
                                                                            sources_dic,
                                                                            sources_bar)
        DataExtracter.check_if_dict_not_empty(filtered_by_sources_orthologs, "after sources filtration")
        print(str(len(filtered_by_sources_orthologs)) + " genes are left " + str(filtered_by_sources_orthologs),
              "after filtration by sources")

        # now orthologs are filtered by sources
        # next step - filtration by length ratio
        filtered_by_length_orthologs = DataExtracter.filter_genes_by_length_differences(filtered_by_sources_orthologs,
                                                                                        length_bar)
        DataExtracter.check_if_dict_not_empty(filtered_by_length_orthologs, "after length filtration")
        print(str(len(filtered_by_length_orthologs)) + " genes are left " + str(filtered_by_length_orthologs),
              "after filtration by length")

        # now we have genes filtered by size and sources.
        # next step: by domains ratio
        filtered_by_conserved_domains = executor.filter_by_conserved_domains_ratio(filtered_by_length_orthologs,
                                                                                   domains_range)
        DataExtracter.check_if_dict_not_empty(filtered_by_conserved_domains, "after domains filtration")
        print(str(len(filtered_by_conserved_domains)) + " genes are left " + str(filtered_by_conserved_domains),
              "after filtration by domains")

        # now we have genes filtered by size, sources and domains ratio.
        # next step: by opposite blast
        true_matches = executor.pair_pipeline(dic_of_optional_orthologs=filtered_by_conserved_domains)
        if not true_matches:
            exit()
        return true_matches

    def get_shinjini_data(self, human_genes_names, c_elegans_genes_names):
        de = DataExtracter()
        for i in range(len(human_genes_names)):
            human_gene_name = human_genes_names[i]
            c_elegans_gene_name = c_elegans_genes_names[i]
            conserved_domains_value = de.get_conserved_domains_ratio_of_pair(c_elegans_gene_name,
                                                                             human_gene_name)
            human_seq = BioPython().get_aa_seq_by_human_gene_name(human_gene_name)
            try:
                c_elegans_id_wb = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
            except:
                print("Couldn't find id for", c_elegans_gene_name)
                continue
            c_elegans_seq = BioPython().\
                get_c_elegans_aa_seq(c_elegans_id_wb)
            conservation_score = BioPython.get_conservation_score(human_seq, c_elegans_seq)
            print(human_gene_name, c_elegans_gene_name, conservation_score, conserved_domains_value)

    # this function is built to answer Ronen's list of genes:
    @staticmethod
    def check_if_gene_has_ortholog(file_path, file_name):
        # first read the genes id
        fd = FileReader(file_path, file_name)
        genes = fd.get_genes_list(0)[1:]  # without the headline
        genes_and_ortholog_data = fd.fromFileToDict(0, 4)
        result_list = executor.find_me_orthologs_for_worm(genes, False, 1, 5, (0.1, 10))
        count = 0
        for worm_gene in genes_and_ortholog_data:
            res = 1 if worm_gene in result_list else 0
            count += res
            print(worm_gene, "input:", genes_and_ortholog_data[worm_gene], "output:", res)
        print(count, "had ortholog out of", len(genes_and_ortholog_data))

    # pipeline that provides you with humans genes that are orthologous for your worm ones
    @staticmethod
    def find_me_orthologs_for_worm(list_of_worm_genes,
                                   genes_in_names: bool = True,
                                   sources_bar: int = 3,
                                   length_bar: int = 10,
                                   domains_range: tuple = (0.5, 2)):
        if genes_in_names:
            list_of_worm_genes_names = list_of_worm_genes
            list_of_human_genes = DataExtracter().convert_list_for_c_elegans(list_of_worm_genes_names)
            print("Genes Ids:" + ", ".join(list_of_human_genes))

        # from C.elegans gene id to human gene id dictionary
        orthologs_dic = FileReader(FileReader.research_path + r"\Data",
                                   r"\ortholist_master",
                                   FileType.TSV).fromFileToDictWithPluralValues(0, 4, True)
        relevant_orthologs_dic = DataExtracter.get_specific_dic_of_orthologs(list_of_worm_genes, orthologs_dic)
        DataExtracter.check_if_dict_not_empty(relevant_orthologs_dic)
        print("Orthologs: " + str(relevant_orthologs_dic))

        # now we have the ortholist-orthologs to our human genes.
        # next step - filtration by sources
        # sources dic is a dictionary with a tuple-keys of (human gene id, c.elegans gene id) and values of number of
        # sources supporting the pair is homologous
        sources_dic = FileReader(FileReader.research_path + r"\Data",
                                 r"\ortholist_master",
                                 FileType.TSV).fromFileToDictWithTupleKey(0, 4, 6, True)
        filtered_by_sources_orthologs = DataExtracter.filter_dic_by_sources(relevant_orthologs_dic,
                                                                            sources_dic,
                                                                            sources_bar)
        DataExtracter.check_if_dict_not_empty(filtered_by_sources_orthologs, "after sources filtration")
        print(str(len(filtered_by_sources_orthologs)) + " genes are left " + str(filtered_by_sources_orthologs),
              "after filtration by sources")

        # now orthologs are filtered by sources
        # next step - filtration by length ratio
        filtered_by_length_orthologs = DataExtracter.filter_genes_by_length_differences(filtered_by_sources_orthologs,
                                                                                        length_bar, "c_elegans")
        DataExtracter.check_if_dict_not_empty(filtered_by_length_orthologs, "after length filtration")
        print(str(len(filtered_by_length_orthologs)) + " genes are left " + str(filtered_by_length_orthologs),
              "after filtration by length")

        # now we have genes filtered by size and sources.
        # next step: by domains ratio
        filtered_by_conserved_domains = executor.filter_by_conserved_domains_ratio(filtered_by_length_orthologs,
                                                                                   domains_range)
        DataExtracter.check_if_dict_not_empty(filtered_by_conserved_domains, "after domains filtration")
        print(str(len(filtered_by_conserved_domains)) + " genes are left " + str(filtered_by_conserved_domains),
              "after filtration by domains")

        # now we have genes filtered by size, sources and domains ratio.
        # next step: by opposite blast
        true_matches = executor.pair_pipeline(dic_of_optional_orthologs=filtered_by_conserved_domains, key_species="human")
        if not true_matches:
            exit()
        return true_matches

exec = executor()

# executor.c_elegans_conserved_domains_file()
# executor.human_conserved_domains_file()
# executor.addDomainsScoreInfoToFile()

# executor.getOnlyGenesWithHumanOrtholog()

# fixing the filtration by size:
# executor.filterGenesBySize(FileReader.research_path + r"\Data",
#                            r"\human_cd_length.txt",
#                            FileReader.research_path + r"\Data",
#                            r"\c_elegans_genes_cds_length.txt"
#                            "genes-orthologs-phenotypes-filteredBySize-100619")
# executor.human_conserved_domains_file()
# executor.c_elegans_conserved_domains_file()
# executor.addDomainsScoreInfoToFile()
# executor.checkExtraGenes()
# executor.fromAccessionsToBLAST()

# executor.getOnlyGenesWithHumanOrtholog()
# executor.filterGenesAccordingToReversedBlast([FileReader.research_path + r"\Executors\true-homologs0-30"],
#                                               FileReader.research_path + r"\Executors\genes_id-genes_names-orthologs-conditions-110619",
#                                               "filterized-data-190619")

# executor.fixConservedDomainRatioFile(FileReader.research_path + r"\Data",
#                                      r"\data-190619-short",
#                                      True)

# executor.add_c_elegans_phenotypes(FileReader.research_path + r"\Executors",
#                                r"\data-190619-fixed-domains",
#                                "http://rest.wormbase.org/rest/widget/gene/",
#                                "data-230619-fixed-ratio-with-C-elegans-phenotypes",
#                                False)

# executor.pair_pipeline(FileReader.research_path + r"\Data",
#                        r"\positive-control-orthologs-pairs")

# executor.filterGenesTest()

# filter_out_duplicated_genes(FileReader.research_path + r"\Data\genes-to-hagar",
#                             r"\genes-sent-to-Hagar",
#                             FileReader.research_path + r"\Data\genes-to-hagar",
#                             r"\rest-genes",
#                             FileReader.research_path + r"\Data\genes-to-hagar\genes-to-hagar-5-duplication-filtered")
#
#

# executor.add_mmp_record_to_data(FileReader.research_path + r"\Data",
#                                 r"\data-220719",
#                                 9,
#                                 [11, 18, 12],
#                                 "data-with-mmp-220719")

# lst = ['KIF5B', 'TNNI1', 'SF3A2', 'EMC8', 'SLCO4C1']

# executor.get_variants_dict(FileType.FILE)

# exec.get_shinjini_data(human_genes_names=['CAPZA1', 'CAPZA2', 'CAPZB'],
#                        c_elegans_genes_names=['cap-1', 'cap-1', 'cap-2'])

print(exec.find_me_orthologs_for_worm(['WBGene00002240'], False, sources_bar=1))
