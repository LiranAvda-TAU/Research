from typing import NamedTuple
from Code.CRISPR.NamedTuples.RestrictionSite import RestrictionSite


class RestrictionMutation(NamedTuple):
    restriction_site: RestrictionSite
    number_of_mutations: int
    mutated_sites: list
