from Code.Files.FileReader import FileReader
from Code.Files.FileMaker import FileMaker


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