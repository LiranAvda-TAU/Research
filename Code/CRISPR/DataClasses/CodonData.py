from typing import NamedTuple

from Code.CRISPR.DataClasses.SequenceSites import SequenceSites


class CodonData(NamedTuple):
    codon: str
    codon_sites: SequenceSites




