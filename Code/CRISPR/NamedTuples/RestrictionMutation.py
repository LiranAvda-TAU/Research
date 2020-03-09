from typing import NamedTuple
from Code.CRISPR.NamedTuples.RestrictionSite import RestrictionSite


class RestrictionMutation(NamedTuple):
    mutated_strand: str
    restriction_site: RestrictionSite
    number_of_mutations: int
    mutated_sites: list
