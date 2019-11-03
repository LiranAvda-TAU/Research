from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.BioPython import BioPython
from Code.Utils.Ensembl import Ensembl


class DataExtracter:

    # receives two dict: one converting human gene to c.elegans gene and one converting human gene to its variants
    # not sure
    @staticmethod
    def find_genes_variants_and_homologous(homologousDic: dict, variantasDic: dict):
        unitedDic = {}
        count = 0
        keys = variantasDic.keys()
        for key in keys:
            humanGene = key[0]
            variantName = key[1]
            if humanGene in homologousDic:
                unitedDic[humanGene, variantName] = [homologousDic[humanGene]] + variantasDic[key]
            else:
                if count < 100:
                    print("key: " + str(key) + str(variantasDic[key]))
                    count += 1
        return unitedDic

    # receivrd a dict of human genes and their conditions as values, and returns a set of all possible conditions
    @staticmethod
    def get_conditions_set(conditionsFilteredDict: dict, condition_column):
        s = set()
        values = conditionsFilteredDict.values()
        for value in values:
            s.add(value[condition_column])

        return s

    # receives two dicts: one converting HGNC gene to its conditions and one converting human gene to c.elegans gene
    # and returns a dict with HGNC name for key and orthologs and conditions as values
    @staticmethod
    def find_genes_with_valid_condition_and_homolog(HGNCGenesWConditions: dict, HomologousGenes: dict):
        genes = {}
        for HCGNGene in HGNCGenesWConditions.keys():  # gene has valid condition
            if HCGNGene in HomologousGenes:  # gene has homoloug gene
                # new dict will have HCGN gene name for key and homoloug genes + condition as value
                genes[HCGNGene] = HomologousGenes[HCGNGene] + [HGNCGenesWConditions[HCGNGene][0]]
        return genes

    # receives (1) a dictionary containing gene name as keys, and orthologs and phenotypes as values, (2) a
    # dictionary for human genes and their lengths, (3) dictionary for c.elegans genes and
    # their lengths and deletes genes for which the length differences between the human gene and its ortholog
    # are to high
    @staticmethod
    def filter_genes_with_size_differences(d: dict, human_genes: dict, c_elegans_genes: dict, p: int):
        genesAndOrthologs = {}
        genesAndPhenotypes = {}
        d_keys = list(d.keys())
        for key in d_keys:
            try:
                humanGeneLength = float(human_genes[key])
            except:
                print("gene " + key + " wasn't found...")
                continue
            d_values = d[key]
            for value in d_values:
                if isinstance(value, tuple):  # item is an ortholog
                    c_elegans_gene = value[0]
                    try:
                        cElegansGeneLength = float(c_elegans_genes[c_elegans_gene])
                    except:
                        print("gene " + c_elegans_gene + " wasn't found...")
                        continue
                    if (humanGeneLength * p) / 100 <= cElegansGeneLength:
                        if key in genesAndOrthologs:
                            genesAndOrthologs[key].append(value)
                        else:
                            genesAndOrthologs[key] = [value]
                    else:  # genes have too different sizes
                        print("C.elegans gene's length is " + str((cElegansGeneLength * 100) / humanGeneLength) +
                              " percent of the human gene's length")
                elif isinstance(value, str):
                    genesAndPhenotypes[key] = value
                else:
                    print("Problem has occured in function: filter_genes_with_size_differences")
                    exit()
        return genesAndOrthologs, genesAndPhenotypes

    # receives (1) a dictionary of c.elegans genes id (number) as keys and accession numbers as values, enact the
    # function that chooses the longest isoform of all accession number, and returns a dictionary with c.elegans gene
    # ids and the chosen longest accession number
    @staticmethod
    def from_multiple_accessions_to_one(genesAndAccessions: dict):
        genesAndChosenAccessions = {}
        keys = genesAndAccessions.keys()
        count = 0
        for gene_id in keys:
            count += 1
            accessions_set = genesAndAccessions[gene_id]
            chosen_isoform = BioPython.from_multiple_accessions_to_one_single_id(accessions_set)
            genesAndChosenAccessions[gene_id] = chosen_isoform
            if count % 100 == 0:
                print("got to", str(count), "with", gene_id)
        return genesAndChosenAccessions

    # receives (1) list of human/c.elegans genes id (WB\ENSG) and (2) a term to search the needed element, make an
    # http request to NCBI to extract the number of domains for the gene, and returns a dictionary of the gene id
    #  and its number of domains
    @staticmethod
    def get_conserved_domains(genes: list):
        genesAndConservedDomain = {}
        parsed_number = DataExtracter().get_conserved_domain_per_gene(genes)
        for gene in genes:
            if parsed_number is not None:
                genesAndConservedDomain[gene] = parsed_number
        return genesAndConservedDomain

    # receives (1) a gene id(WB/ENSG) and (2),(3) terms to extract only the needed number out of the line, makes an
    # http request to extract the domains of each gene, parses it to an integer and returns it
    @staticmethod
    def get_conserved_domain_per_gene(gene, start_term: str = "Conserved Domains (", end_term: str = ")"):
        hr = HttpRequester("https://www.ncbi.nlm.nih.gov/gene/?term=" + gene)
        data = hr.makeRequest()
        if isinstance(data, str):  # request was successful
            index = data.find(start_term)
            shorter_data = data[index + len(start_term):]
            end_index = shorter_data.find(end_term)
            number = shorter_data[:end_index]
            try:
                parsed_number = int(number)
            except ValueError:
                print("for gene: " + gene + " letter extracted is " + number)
                parsed_number = None
        return parsed_number

    # receives a worm gene name and a human gene name, get conserved domains for both genes, and returns their ratio
    @staticmethod
    def get_conserved_domains_ratio_of_pair(c_elegans_gene_name, human_gene_name):
        print("c elegance gene: " + c_elegans_gene_name + " and human gene: " + human_gene_name)
        try:
            c_elegans_id = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
        except:
            print("no id for gene:", human_gene_name)
            return None
        try:
            human_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        except:
            print("no id for gene:", human_gene_name)
            return None
        if c_elegans_id is None or human_id is None:
            print("id not found")
            return None

        de = DataExtracter()
        c_elegans_domains = de.get_conserved_domain_per_gene(c_elegans_id)
        human_domains = de.get_conserved_domain_per_gene(human_id)
        return float(c_elegans_domains) / float(human_domains)

    # receives worm gene name and human gene name, reads the ortholist file to a dictionary with
    # (worm gene id, human gene id) keys and number of sources as values, and returns the number of sources if the tuple
    # is in the dictionary, otherwise returns -1
    @staticmethod
    def get_sources(c_elegans_gene_name, human_gene_name):
        c_elegans_gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
        human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)

        orthologs = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                               r"\ortholist_master",
                               FileType.TSV).fromFileToDictWithPluralValues(0, 4, True)
        sources = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                             r"\ortholist_master",
                             FileType.TSV).fromFileToDictWithTupleKey(0, 4, 6, True)

        if c_elegans_gene_id in orthologs:
            if human_gene_id in orthologs[c_elegans_gene_id]:
                return int(sources[(c_elegans_gene_id, human_gene_id)])
            else:
                print(orthologs[c_elegans_gene_id])
                return -1
        return -1

    # receives a (1) dictionary with gene name as key and the gene id as value, and (2) a list of gene names, and
    # returns a list of the genes ids
    @staticmethod
    def from_human_genes_names_to_human_genes_id(genes_names: list):
        genes_ids = []
        for gene_name in genes_names:
            gene_id = Ensembl.get_human_gene_id_by_gene_name(gene_name)
            if gene_id:
                genes_ids.append(gene_id)
        return genes_ids

    # receives (1) a dictionary of c.elegans genes and number of domains, (2) a dictionary of human genes ids and
    # number of domains, (3) a dictionary with human genes ids as keys and c.elegans orthologs as values, and
    # returns a new dic with the human genes and their orthologs if they are inside the range of ratio
    @staticmethod
    def filter_by_conserved_domains(c_elegans_domains_dic: dict,
                                    human_domains_dic: dict,
                                    orthologs_dic: dict,
                                    domains_range):
        new_orthologs_dic = {}
        human_genes = orthologs_dic.keys()
        for human_gene_id in human_genes:
            try:
                c_elegans_orthologs = orthologs_dic[human_gene_id]
            except KeyError:
                print("human gene id", human_gene_id, "doesn't exist in the file of orthologs")
                continue
            if human_gene_id in human_domains_dic:
                human_domains = human_domains_dic[human_gene_id]
            else:
                human_domains = DataExtracter().get_conserved_domain_per_gene(human_gene_id)

            for ortholog in c_elegans_orthologs:
                try:
                    c_elegans_domains = c_elegans_domains_dic[ortholog]
                except:
                    c_elegans_domains = DataExtracter().get_conserved_domain_per_gene(ortholog)
                ratio = float(c_elegans_domains) / float(human_domains)
                if domains_range[0] <= ratio <= domains_range[1]:
                    print("Ratio between", ortholog, "and",  human_gene_id, "is", str(ratio), ", thus passes the"
                          " domains ratio bar!")
                    DataExtracter.add_to_dictionary(new_orthologs_dic, human_gene_id, ortholog)
                else:
                    print("Ratio between", human_gene_id, "and", ortholog + " is", str(ratio), ", thus not good enough")
        return new_orthologs_dic

    # receives (1) a dictionary of c.elegans genes and number of domains, (2) a dictionary of human genes ids and
    # number of domains, (3) a dictionary with human genes names as keys and c.elegans orthologs as values, and (4)
    # a dictionary of human genes names as keys and human genes ids as values, and adds to the orthologs the ratio
    # of domains between C.elegans and humans
    @staticmethod
    def add_conserved_domains_info(c_elegans_domains_dic: dict, human_domains_dic: dict, orthologs_dic: dict,
                                   human_name_id_dic: dict):
        new_orthologs_dic = {}
        human_genes = orthologs_dic.keys()
        for human_gene_name in human_genes:
            try:
                human_gene_id = human_name_id_dic[human_gene_name]
            except KeyError:
                print(human_gene_name + " has no equivalent gene id in the dictionary")
                continue
            try:
                c_elegans_orthologs = FileReader.fromStringToTuplesList(orthologs_dic[human_gene_name])
            except KeyError:
                print(human_gene_name + "doesn't exist in the  file of orthologs")
                continue
            if human_gene_id in human_domains_dic:
                human_domains = human_domains_dic[human_gene_id]
                orthologs_with_domains = DataExtracter.add_domains_similarity_score(int(human_domains),
                                                                                    c_elegans_domains_dic,
                                                                                    c_elegans_orthologs)
                new_orthologs_dic[human_gene_id] = orthologs_with_domains
            else:
                print("human gene " + human_gene_id + " doesn't exist in the domains dictionary")
                pass
        return new_orthologs_dic

    # receives (1) number of human gene domains, (2) a dictionary of c.elegans genes ids (WB?) and number of domains as
    # values and (3) list of C.elegans orthologs, and adds the ratio of number of domains between each C.elegans
    # ortholog gene and the human gene domains to the list of orthologs
    @staticmethod
    def add_domains_similarity_score(human_domains: int, c_elegans_domains_dic: dict, orthologs: list):
        orthologs_with_domains = []
        for ortholog, sources in orthologs:
            if ortholog in c_elegans_domains_dic:
                ortholog_domains = int(c_elegans_domains_dic[ortholog])
                domainsScore = min(ortholog_domains, human_domains) / max(ortholog_domains, human_domains)
                orthologs_with_domains.append((ortholog, sources, domainsScore))
            else:
                print("c elegans gene " + ortholog + " doesn't exist in dictionary of c.elegans domains")
        return orthologs_with_domains

    # receives (1) the number of results files, (2) the file path, (3) the file name, and (4) an optional start file
    # index, goes through the files and extract for each accession number the given hit ids and their stats, and returns
    # adictionary of each accession number and their hit ids, and each hit id and it stats
    @staticmethod
    def get_hit_ids_for_accession_numbers(end_file: int, file_path, file_name, start_file: int = 0):
        accession_numbers_and_hit_ids = {}
        hit_ids_and_hsp = {}
        fd = FileReader(file_path, file_name, FileType.TSV)
        accession_number = "default"
        hit_id = "default"
        for result_index in range(start_file, end_file):
            f = fd.readResultsFile()
            print("created a file to extract the results!")
            line = f.readline()  # accession-number...
            while line:  # line is not an empty string, thus we haven't reached the end of the file yet
                if line.startswith("Accession Number"):
                    accession_number = line.rstrip()[len("Accession Number: "):]
                    print("now working on accession number: " + accession_number)
                    line = f.readline()
                elif line.startswith("****Alignment****"):
                    next_line = f.readline().rstrip()
                    hit_id = next_line[next_line.find("|") + 1: next_line.find("|", next_line.find("|") + 1)]
                    extra_letter = next_line[next_line.find(hit_id) + len(hit_id) + 1:
                    next_line.find(hit_id) + len(hit_id) + 2]
                    if extra_letter:
                        hit_id = hit_id + "_" + extra_letter
                    DataExtracter.add_to_dictionary(accession_numbers_and_hit_ids, accession_number, hit_id)
                    line = f.readline()
                elif line.startswith("****HSP****"):
                    e_value = f.readline().rstrip()[len("e value: "):]
                    score = f.readline().rstrip()[len("hsp score: "):]
                    DataExtracter.add_to_dictionary(hit_ids_and_hsp, hit_id, e_value + "+" + score)
                    line = f.readline()
                else:
                    line = f.readline()
            f.close()
        return accession_numbers_and_hit_ids, hit_ids_and_hsp

    # receives a (1) dictionary, (2) a key, and (3) a value, and added the value to a list of values for the
    # given key
    @staticmethod
    def add_to_dictionary(d: dict, key: str, value):
        if key in d:
            d[key].append(value)
        else:
            d[key] = [value]

    # receives (1) a dictionary with list of values for each key, and return a reversed dictionary
    @staticmethod
    def from_plural_valued_dict_to_plural_valued_reversed_dict(d: dict):
        reversed_dict = {}
        keys = d.keys()
        for key in keys:
            value_list = d[key]
            for value in value_list:
                DataExtracter().add_to_dictionary(reversed_dict, value, key)
        return reversed_dict

    # receives (1) the data dictionary of orthologs and (2) a domain multiplying factor, creates a score for each
    # C.elegans gene and returns a dictionary of C.elegans genes and their orthology score
    @staticmethod
    def scoring_function(orthologs_dictionary: dict, domain_factor):
        c_elegans_genes_and_scores = {}
        human_genes = list(orthologs_dictionary.keys())
        for human_gene in human_genes:
            c_elegans_orthologs = FileReader.fromStringToTuplesList(orthologs_dictionary[human_gene])
            for homolog in c_elegans_orthologs:
                gene_name = homolog[0]
                sources = homolog[1]
                domain_equality = homolog[2]
                score = int(sources) + domain_factor * float(domain_equality)
                c_elegans_genes_and_scores[gene_name] = score

    @staticmethod
    def add_description_from_worm_base(data_file_path, new_data_file_path, column, empty_description, ortholog_column,
                                       start_term: str = "\"text\":\"", end_term: str = "\",\"evidence\""):
        in_file = open(data_file_path)
        out_file = open(new_data_file_path, FileMode.WRITE.value)
        for line in in_file:
            lst = line.strip("\n").split("\t")
            description = lst[column]
            if description == empty_description:
                record = HttpRequester("http://rest.wormbase.org/rest/widget/gene/" +
                                       lst[ortholog_column] + "/overview").makeRequest()
                info_index = record.find("concise_description")
                shorter_info = record[info_index + len("concise_description"):]
                start_index = shorter_info.find(start_term)
                end_index = shorter_info.find(end_term)
                description = shorter_info[start_index + len(start_term): end_index]

                for i in range(column):
                    out_file.write(lst[i] + "\t")
                out_file.write(description + "\n")
            else:
                out_file.write(line)
        in_file.close()
        out_file.close()

    # receives (1) list of genes, (2) url address to the site where we can find phenotypes, and (3),(4) terms to extract
    # only the relevant phenotypes from the line given, and returns a dictionary of genes and phenotypes
    @staticmethod
    def add_c_elegans_phenotypes(list_of_genes: list, url, search_term: str = "\"class\":\"phenotype\"",
                                 label_term: str = "\"label\":"):
        genes_and_phenotypes = {}
        for gene in list_of_genes:
            print("gene: " + gene)
            info = HttpRequester(url + gene + "/" + "phenotype").makeRequest()
            search_index = info.find(search_term)
            while -1 < search_index < info.find("\"phenotype_not_observed\":") and \
                            search_index < info.find("\"phenotype_by_interaction\""):
                info = info[search_index:]
                label_index = info.find(label_term)
                info = info[label_index:]
                phenotype = info[len(label_term):info.find(",")].strip("\"")
                print(phenotype)
                if gene in genes_and_phenotypes:
                    genes_and_phenotypes[gene].append(phenotype)
                else:
                    genes_and_phenotypes[gene] = [phenotype]
                search_index = info.find(search_term)

        return genes_and_phenotypes

    @staticmethod
    def get_hit_ids(c_elegans_gene_id, c_elegans_gene_name, c_elegans_id_multiple_accessions, accession_number_to_hit_ids_dic):
        try:
            c_elegans_id_number = Ensembl.get_ncbi_id_by_gene_id(c_elegans_gene_id)
            print("gene id is: " + c_elegans_id_number)
        except:
            print("Couldn't find gene id number for gene: " + c_elegans_gene_name)
            exit()
        try:
            c_elegans_accession_numbers = set(c_elegans_id_multiple_accessions[c_elegans_id_number])
            c_elegans_accession_number = BioPython.from_multiple_accessions_to_one_single_id(
                c_elegans_accession_numbers)
        except:
            c_elegans_accession_number = input("What is the accession number for id: " + c_elegans_id_number + "\n")
        print("Gene's accession number is: " + c_elegans_accession_number)
        try:
            hit_ids = accession_number_to_hit_ids_dic[c_elegans_accession_number]
        except:
            print("got here!")
            seq = BioPython.get_aa_seq(c_elegans_gene_id)
            if not seq:
                seq = BioPython().get_aa_seq_of_longest_isoform(c_elegans_accession_number)
                if not seq:
                    print("Couldn't find gene sequence for gene: " + c_elegans_gene_name)
                    exit()
            hit_ids = BioPython().pipeline_blast_with_seq("blastp", "nr", seq)
        if not hit_ids:
            print("Couldn't find hit ids for", c_elegans_gene_name)
        print("Length of hit ids list is: " + str(len(hit_ids)))
        return hit_ids

    # receives (1) worm gene name, (2) worm gene id(WB), (3) human gene name, (4) dictionary of worm genes id ans their
    #  accessions numbers, (5) dictionary of accession number as key and hit ids as value, (6) dictionary of hit id as
    # keys and human gene name as value and (7) dictionary of id(WB) to id(number), and runs a reverse blast on the worm
    # genes to see if we get a match in the form of the human gene
    # scheme: worm gene id(WB) -> worm gene id(number) -> accession numbers -> accession number -> hit ids -> human
    # gene name -> check if our human gene is in there
    @staticmethod
    def check_reversed_blast_hit_ids(c_elegans_gene_name, c_elegans_gene, human_gene_name,
                                     c_elegans_id_multiple_accessions, accession_number_to_hit_ids_dic,
                                     hit_ids_to_human_genes_names_dic):

        hit_ids = DataExtracter.get_hit_ids(c_elegans_gene,
                                            c_elegans_gene_name,
                                            c_elegans_id_multiple_accessions,
                                            accession_number_to_hit_ids_dic)
        for hit_id in hit_ids:
            print("current hit id:", hit_id)
            try:
                hit_human_gene_name = hit_ids_to_human_genes_names_dic[hit_id]
                if hit_human_gene_name == human_gene_name:
                    print(c_elegans_gene + " and " + human_gene_name + " are orthologs indeed!")
                    return True
                else:
                    print("hit human gene name is", hit_human_gene_name, "thus there is no match")
            except:
                # need to extract gene's name
                try:
                    hit_human_gene_name = BioPython().get_gene_name_from_protein_accession([hit_id])[hit_id]
                    print("The gene's name: " + hit_human_gene_name)
                    if hit_human_gene_name == human_gene_name:
                        print(c_elegans_gene_name + " and " + human_gene_name + " are orthologs indeed!")
                        return True
                except:
                    print("Couldn't find the human gene name to hit id: " + hit_id)
        print("No match between " + c_elegans_gene_name + " and " + human_gene_name)
        return False

    # receives (1) list of genes, (2) path of data file, (3) name of data file, (4) name for filtered file,
    # (5) phenotype column and (6) filtering word, and for each gene, if its name is in the genes list, and the
    # filtering word is in its phenotype description, it copies the info line to the new filtered file
    @staticmethod
    def filter_genes(short_list_genes, data_path, data_name, filtered_name, phenotype_column, filtering_word):
        print("So the function works...")
        f = open(data_path + data_name, FileMode.READ.value)
        filtered_file = open(filtered_name, FileMode.WRITE.value)
        f.readline()
        for line in f:
            info = line.strip("\n").split("\t")
            gene_name = info[0]
            if gene_name in short_list_genes:
                print(gene_name + " is in short list!")
                print(info[phenotype_column])
                if filtering_word in info[phenotype_column]:
                    filtered_file.write(line)
                else:
                    continue
            else:
                continue
        f.close()
        filtered_file.close()

    # receives (1) file path of unwanted genes, (2) said file name, (3) path of data, (4) name of data, (5) new data
    # file name, and copies only the data lines for human genes that are not included in the unwanted genes list, and
    # only if the matching C.elegans gene is not included in the unwanted genes for the human genes, to the new file
    @staticmethod
    def filter_genes_by_name(unwanted_genes_file_path, unwanted_genes_file_name, data_path, data_name, new_data_name):
        unwanted_genes_dic = FileReader(unwanted_genes_file_path,
                                        unwanted_genes_file_name).fromFileToDictWithPluralValues(0, 2, False)
        data = open(data_path + data_name)
        new_data = open(new_data_name, FileMode.WRITE.value)
        data.readline()
        count = 0
        read = 0
        stayed = 0
        for line in data:
            read += 1
            lst = line.strip("\n").split("\t")
            if lst[0] not in unwanted_genes_dic:
                new_data.write("\t".join(lst) + "\n")
                stayed += 1
            else:
                if lst[2] not in unwanted_genes_dic[lst[0]]:
                    new_data.write("\t".join(lst) + "\n")
                    stayed += 1
                else:
                    print(line.strip("\n"))
                    count += 1
        data.close()
        new_data.close()
        print(count)
        print(read)
        print(stayed)

    # receives (1) list of keys, and (2) dictionary of keys and values and returns a new dictionary containing only the
    # keys and values for the keys given in input (1)
    @staticmethod
    def get_specific_dic_of_orthologs(list_of_human_genes, human_to_c_elegans_dic):
        dic = {}
        for gene in list_of_human_genes:
            try:
                dic[gene] = human_to_c_elegans_dic[gene]
            except:
                print("Couldn't find orthologs for", gene)
        return dic

    # receives (1) dictionary of human-c.elegans orthologs, (2) dictionary containing tuples of human gene id and
    # c.elegans gene id as keys and number of sources supporting this is a homologous pair as values, and (3) the
    # efficient number of sources bar, and only the pairs that have sufficiently high number of sources are copied into
    # a new filtered dictionary
    @staticmethod
    def filter_dic_by_sources(orthologs_dic, sources_dic, bar):
        filtered_dic = {}
        keys = orthologs_dic.keys()
        for key in keys:
            for ortholog in orthologs_dic[key]:
                if int(sources_dic[(key, ortholog)]) >= bar:
                    print(key + " with its ortholog " + ortholog + " have passes the sources bar!")
                    DataExtracter.add_to_dictionary(filtered_dic, key, ortholog)

        return filtered_dic

    # receives (1) list of genes names of human or c.elegans and (2) the string "human" or "c.elegans", and returns a
    # new list with the relevant genes ids.
    @staticmethod
    def convert_list(genes_list, subject):
        new_list = []
        for gene in genes_list:
            if subject == "human":
                gene_id = Ensembl.get_human_gene_id_by_gene_name(gene)
            elif subject == "c.elegans":
                gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(gene)
            if gene_id:
                new_list.append(gene_id)
            else:
                print("The name", gene, "is not listed as a", subject, "gene in our sources")
                exit()

        return new_list

    def convert_list_for_human(self, genes_list):
        return self.convert_list(genes_list, "human")

    def convert_list_for_c_elegans(self, genes_list):
        return self.convert_list(genes_list, "C.elegans")

    @staticmethod
    def get_pair_cd_length(human_gene, c_elegans_gene):
        humans_id_cd_length = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                         r"\human_cd_length.txt").get_genes_cd_length(1, 2, True)

        c_elegans_cd_length = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                         r"\c_elegans_genes_cd_length.txt").get_genes_cd_length(1, 2, True)

        try:
            human_gene_length = humans_id_cd_length[human_gene]
        except:
            print("gene " + human_gene + "'s length wasn't found...")
            human_gene_length = None
        try:
            c_elegans_gene_length = c_elegans_cd_length[c_elegans_gene]
        except:
            print("gene " + c_elegans_gene + "'s length wasn't found...")
            c_elegans_gene_length = None
        return human_gene_length, c_elegans_gene_length

    # receives (1) a dictionary containing human genes name as keys, and orthologs as values, (2) a ratio bar, and
    #  deletes genes for which the length differences between the human gene and its ortholog
    # are too high
    @staticmethod
    def filter_genes_by_length_differences(d: dict, p: int):
        humans_id_cd_length = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                         r"\human_cd_length.txt").get_genes_cd_length(0, 2, True)

        c_elegans_cd_length = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                         r"\c_elegans_genes_cd_length.txt").get_genes_cd_length(0, 2, True)

        genes_and_orthologs = {}
        d_keys = list(d.keys())
        for key in d_keys:
            try:
                human_gene_length = float(humans_id_cd_length[key])
            except:
                print("gene " + key + " wasn't found...")
                continue
            d_values = d[key]
            for c_elegans_gene in d_values:
                try:
                    c_elegans_gene_length = float(c_elegans_cd_length[c_elegans_gene])
                except:
                    print("gene " + c_elegans_gene + " wasn't found...")
                    continue
                if (human_gene_length * p) / 100 <= c_elegans_gene_length and \
                                        (c_elegans_gene_length * p) / 100 <= human_gene_length:
                    print("Gene:", key, "has length of", str(int(human_gene_length)), "bp, and its C.elegans",
                          "ortholog", c_elegans_gene, "has length of", str(int(c_elegans_gene_length)),
                          "bp, so they pass the length bar")
                    DataExtracter.add_to_dictionary(genes_and_orthologs, key, c_elegans_gene)
                else:  # genes have too different sizes
                    print("C.elegans gene's length is", c_elegans_gene_length, "and human gene length is",
                          human_gene_length, "and the C.elegans gene is", str((c_elegans_gene_length * 100) / human_gene_length),
                            "percent of the human gene's length")
        return genes_and_orthologs

    # receives (1) a possible-plural-valued-dictionary of worm-human orthologs and (2) a boolean value indicating if we
    # are interested in converting id to name or name to id, and returns a new worm-human dictionary with the converted keys and values
    def convert_dic(self, dic, from_id_to_name):
        converted_dic = {}
        keys = dic.keys()
        for key in keys:
            values = dic[key]
            for value in values:
                if from_id_to_name:
                    DataExtracter().add_to_dictionary(converted_dic,
                                                      Ensembl.get_gene_name_by_gene_id(key),
                                                      Ensembl.get_gene_name_by_gene_id(value))
                else:
                    DataExtracter().add_to_dictionary(converted_dic,
                                                      Ensembl.get_c_elegans_gene_id_by_gene_name(key),
                                                      Ensembl.get_human_gene_id_by_gene_name(value))
        return converted_dic

    # receives (1) a dict and checks if it is empty. if so, it exits the function
    @staticmethod
    def check_if_dict_not_empty(dic: dict, step=""):
        if not dic:
            print("dictionary has no genes left", step)
            exit()

    ########### irrelevant functions ###########

    @staticmethod
    def fix_conserved_domain_info(data_path, data_name, c_elegans_domains_dic, human_domains_dic, domain_column,
                                  fixed_file, delete_first_line: bool = False):
        in_file = open(data_path + data_name)
        out_file = open(fixed_file, FileMode.WRITE.value)
        if delete_first_line:
            in_file.readline()
        for line in in_file:
            row = line.rstrip("\n").split("\t")
            c_elegans_gene = row[2]
            human_gene = row[0]
            row[domain_column] = str(
                float(c_elegans_domains_dic[c_elegans_gene]) / float(human_domains_dic[human_gene]))
            out_file.write("\t".join(row) + "\n")
        in_file.close()
        out_file.close()

    @staticmethod
    def fix_conserved_domains_file():
        human_id_ENSG = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                   r"\human-genes-and-conserved-domains-230619",
                                   FileType.TSV).fromFileToDict(0, 1)
        c_elegans_id_WB = FileReader(r"C:\Users\Liran\PycharmProjects\Research\Data",
                                     r"\c-elegans-genes-and-conserved-domains-230619",
                                     FileType.TSV).fromFileToDict(0, 1)
        f = open(r"C:\Users\Liran\PycharmProjects\Research\Executors\data-230619-fixed-ratio-with-C-elegans-phenotypes")
        new_file = open("data-250619-fixed-domains", FileMode.WRITE.value)
        f.readline()
        for line in f:
            info = line.strip("\n").split("\t")
            human_id = info[0]
            c_elegans_id = info[2]
            info[4] = str(float(c_elegans_id_WB[c_elegans_id]) / float(human_id_ENSG[human_id]))
            new_file.write("\t".join(info) + "\n")
        f.close()
        new_file.close()

