from typing import NamedTuple


class RestrictionMutation(NamedTuple):
    mutated_strand: str
    restriction_sites: dict
    number_of_mutations: int
    mutated_sites: list
