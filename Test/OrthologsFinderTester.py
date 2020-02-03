from Executors import executor

list_of_human_genes_names = 'A2M', 'ABL1', 'ADCY5', 'AGPAT2'
list_of_human_genes_ids = 'ENSG00000164330', 'ENSG00000183878', 'ENSG00000099246', 'ENSG00000188487'
list_of_c_elegans_genes_names = 'linc-144', 'T27B1.9', 'F46H6.10', 'C42D8.15', 'samt-1', 'srz-89', 'fbxc-18', 'cec-7', 'glb-14'
list_of_c_elegans_genes_ids = 'WBGene00000287', 'WBGene00017010', 'WBGene00047390', 'WBGene00047453', 'WBGene00047647', 'WBGene00165664', 'WBGene00167828', 'WBGene00195616', 'WBGene00219694'


def test_find_me_orthologs_for_human_by_name():
    executor.executor.find_me_orthologs_for_human(list_of_human_genes_names)


def test_find_me_orthologs_for_human_by_id():
    executor.executor.find_me_orthologs_for_human(list_of_human_genes_ids, False)


def test_find_me_orthologs_for_worm_by_name():
    executor.executor.find_me_orthologs_for_worm(list_of_c_elegans_genes_names)


def test_find_me_orthologs_for_worm_by_id():
    executor.executor.find_me_orthologs_for_worm(list_of_c_elegans_genes_ids, False)


test_find_me_orthologs_for_human_by_name()
test_find_me_orthologs_for_human_by_id()
test_find_me_orthologs_for_worm_by_name()
test_find_me_orthologs_for_worm_by_id()
