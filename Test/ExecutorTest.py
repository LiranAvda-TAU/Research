from Code.Enum.FileType import FileType
from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileReader import FileReader
from Code.Utils import Extra
from Executors.executor import executor

exec = executor()

# executor.c_elegans_conserved_domains_file()
# executor.human_conserved_domains_file()
# executor.add_domains_score_info_to_file()

# executor.get_only_genes_with_human_ortholog()

fix_filtration = False
if fix_filtration:
    # fixing the filtration by size:
    executor.filter_genes_by_size(FileReader.research_path + r"\Data",
                               r"\human_cd_length.txt",
                               FileReader.research_path + r"\Data",
                               r"\c_elegans_genes_cds_length.txt"
                               "genes-orthologs-phenotypes-filteredBySize-100619")
    executor.human_conserved_domains_file()
    executor.c_elegans_conserved_domains_file()
    executor.add_domains_score_info_to_file()
    executor.check_extra_genes()
    executor.from_accessions_to_blast()


# executor.get_only_genes_with_human_ortholog()
# executor.filterGenesAccordingToReversedBlast([FileReader.research_path + r"\Executors\true-homologs0-30"],
#                                               FileReader.research_path + r"\Executors\genes_id-genes_names-orthologs-conditions-110619",
#                                               "filterized-data-190619")

fix_conserved_domains_ratio = False
if fix_conserved_domains_ratio:
    executor.fix_conserved_domain_ratio_file(FileReader.research_path + r"\Data",
                                         r"\data-190619-short",
                                         True)

add_phenotypes = False
if add_phenotypes:
    executor.add_c_elegans_phenotypes(FileReader.research_path + r"\Executors",
                                   r"\data-190619-fixed-domains",
                                   "http://rest.wormbase.org/rest/widget/gene/",
                                   "data-230619-fixed-ratio-with-C-elegans-phenotypes",
                                   False)

check_orthologs_pairs = False
if check_orthologs_pairs:
    print("#CHECKING ORTHOLOGS PAIRS#")
    true_matches = executor.pair_pipeline(FileReader.research_path + r"\Data", r"\positive-control-orthologs-pairs")
    # true_matches = executor.pair_pipeline(FileReader.research_path + r"\Data", r"\positive-control-orthologs-pairs-test")
    print("results:")
    for pair in true_matches:
        human_length, worm_length = true_matches[pair][2]
        print(pair, human_length, worm_length, worm_length/human_length*100, sep="\t")

# executor.filterGenesTest()

filter_out_duplicates = False
if filter_out_duplicates:
    Extra.filter_out_duplicated_genes(FileReader.research_path + r"\Data\genes-to-hagar",
                                r"\genes-sent-to-Hagar",
                                FileReader.research_path + r"\Data\genes-to-hagar",
                                r"\rest-genes",
                                FileReader.research_path + r"\Data\genes-to-hagar\genes-to-hagar-5-duplication-filtered")



add_mmp_records = False
if add_mmp_records:
    executor.add_mmp_record_to_data(FileReader.research_path + r"\Data",
                                    r"\data-220719",
                                    9,
                                    [11, 18, 12],
                                    "data-with-mmp-220719")

# lst = ['KIF5B', 'TNNI1', 'SF3A2', 'EMC8', 'SLCO4C1']

get_variants_dict = False
if get_variants_dict:
    executor.get_variants_dict(FileType.FILE)

# help_shinjini = False
# if help_shinjini:
#     exec.get_shinjini_data(human_genes_names=['CAPZA1', 'CAPZA2', 'CAPZB'],
#                            c_elegans_genes_names=['cap-1', 'cap-1', 'cap-2'])


# print(exec.find_me_orthologs_for_worm(['WBGene00013355'], False, sources_bar=1))

priti_request = False
if priti_request:
    exec.check_if_gene_has_ortholog(file_path=FileReader.research_path + r"\Data",
                                    file_name=r"\C.elegans-kinase-phosphatase-genes.xlsx",
                                    sheet_name="kinase")

improved_check_control_orthologs_pairs = False
if improved_check_control_orthologs_pairs:
#CHECKING POSITIVE CONTROL ORTHOLOG PAIRS
    fd = FileReader(FileReader.research_path + r"\Data", r"\positive-control-orthologs-pairs")
    c_elegans_human_orthologs = fd.from_file_to_dict(delete_first_line=True)
    for worm in c_elegans_human_orthologs:
        worm_length, human_length = DataExtracter().get_pair_cd_length(c_elegans_human_orthologs[worm], worm)
        print(c_elegans_human_orthologs[worm], worm, human_length, worm_length, worm_length/human_length*100, sep="\t")
