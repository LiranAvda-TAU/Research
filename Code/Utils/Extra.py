from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileReader import FileReader
from Code.Files.FileMaker import FileMaker
from Code.Utils.Ensembl import Ensembl


def splitDict(d):
    split_index = len(d) // 2
    items = list(d.items())
    d1 = items[:split_index]
    d2 = items[split_index:]
    return d1, d2


def filter_out_duplicated_genes(sent_file_path, sent_path_name, in_file_path, in_file_name, out_file_name):
    filtered_lst = []
    sent_genes = FileReader(sent_file_path, sent_path_name).from_file_to_list(False)
    lst = FileReader(in_file_path, in_file_name).from_file_to_list(False)
    # to filter out self-duplicates
    lst_without_duplicates = list(set(lst))
    for gene in lst_without_duplicates:
        if gene not in sent_genes:
            filtered_lst.append(gene)
    FileMaker().fromListToFile(filtered_lst, out_file_name)


def priti_request(file_path1, file_name1, file_path2, file_name2):
    data1 = {}
    data2 = {}
    worm_ids1 = FileReader(file_path1, file_name1).get_genes_list(0)
    worm_ids2 = FileReader(file_path2, file_name2).get_genes_list(0)
    for worm_id in worm_ids1:
        gene_name = Ensembl.get_gene_name_by_gene_id(worm_id)
        gene_description = DataExtracter.get_c_elegans_description_for_gene_id(worm_id)
        dtc = True if "dtc" in gene_description or "distal tip cell" in gene_description else False
        data1[worm_id] = (gene_name, gene_description, dtc)
    for worm_id in worm_ids2:
        gene_name = Ensembl.get_gene_name_by_gene_id(worm_id)
        gene_description = DataExtracter.get_c_elegans_description_for_gene_id(worm_id)
        dtc = True if "dtc" in gene_description or "distal tip cell" in gene_description else False
        data2[worm_id] = (gene_name, gene_description, dtc)

    common = set()
    specific = set()
    for gene in data1:
        if gene in data2:
            common.add(gene)
        else:
            specific.add(gene)
    for gene in data2:
        if gene not in data1:
            specific.add(gene)

# priti_request(r"C:\Users\Liran\Downloads", r"\priti.txt")
