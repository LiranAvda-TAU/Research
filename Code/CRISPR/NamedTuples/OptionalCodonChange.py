from typing import NamedTuple


class OptionalCodonChange(NamedTuple):
    codon: str
    number_of_mutations: int
    dict_of_mutations: dict
    usage: float



