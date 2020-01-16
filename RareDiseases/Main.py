from Code.Files.FileMaker import FileMaker
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester


def getNumberOfAffectedPopulations(file_path, file_name, url, start_term, end_term):
    diseases_and_numbers = {}
    diseases = FileReader(file_path, file_name).from_file_to_list(True)
    for disease in diseases:
        info: str = HttpRequester(url + putHyphenInName(disease)).make_request()
        start_index = info.find(start_term)
        end_index = info.find(end_term)
        data = info[start_index + len(start_term):end_index]
        number = input("for disease: " + disease + ":" + data + "\n")
        diseases_and_numbers[disease] = number
    FileMaker().from_dict_to_file(diseases_and_numbers, "diseases-and-affected-populations")


def putHyphenInName(name: str):
    hyphened_name = ""
    hyphened_name += name[0]
    for letter in name[1:]:
        if letter.isupper():
            hyphened_name += "-"
            hyphened_name += letter
        elif letter == " ":
            continue
        else:
            hyphened_name += letter
    return hyphened_name


getNumberOfAffectedPopulations(FileReader.research_path + r"\RareDiseases",
                               r"\diseases",
                               "https://rarediseases.org/rare-diseases/",
                               "wpcf-field-rd_affected_populations",
                               "id=\"related-disorders\"")
