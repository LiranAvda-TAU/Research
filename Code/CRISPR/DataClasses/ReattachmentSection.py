from Code.CRISPR.DataClasses.SequenceSites import SequenceSites
from Code.CRISPR.Enum.DNASection import DNASection
from typing import NamedTuple


class ReattachmentSection(NamedTuple):
    section_type: DNASection
    section_sites: SequenceSites
    number_of_mutations: int



