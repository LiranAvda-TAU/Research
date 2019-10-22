from Code.Enum.FileType import FileType
from Executors.executor import executor
import sys


def main(program):

    if program == 'Homology':
        ###### HOMOLOGOUS PIPELINE SECTION ########
        # input = ['TCP1']
        true_matches = executor.find_me_orthologs(['TCP1F'],
                                                  genes_in_names=True,
                                                  sources_bar=3,
                                                  length_bar=10,
                                                  domains_range=(0.3, 2))
        print("We are left with " + str(len(true_matches)) + " matches")
        print(str(true_matches))

    elif program == 'Variants':
        ### VARIANTS ###
        executor.get_variants_data(FileType.CONSOLE,
                                   {'TCP1': ['Asn284Ser', 'Ala453Glu']})
    else:
        print("No program called", sys.argv[1])
        exit()

main(sys.argv[1])
