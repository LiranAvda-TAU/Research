from typing import NamedTuple
from Code.CRISPR.NamedTuples.SequenceSites import SequenceSites

class RestrictionSite(NamedTuple):
    index: SequenceSites
    site: str
