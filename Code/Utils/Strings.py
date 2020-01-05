class Strings:

    # receives a string and return a tuple of (former amino acid, place, current amino acid)
    @staticmethod
    def from_variant_string_to_tuple(s):
        former_amino_acid = ''
        current_amino_acid = ''
        place = ''
        next_amino_acid = False
        for char in s:
            if 'A' <= char <= 'z' and not next_amino_acid:
                former_amino_acid += char
            elif 'A' <= char <= 'z' and next_amino_acid:
                current_amino_acid += char
            else:
                place += char
                next_amino_acid = True
        place = int(place)
        return former_amino_acid, place, current_amino_acid

    @staticmethod
    def fromNameToSymbol(amino_acid_name):
        symbols = {'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                   'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S',
                   'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y', 'Term': '-'}
        return symbols[amino_acid_name]

    @staticmethod
    def fromFastaSeqToSeq(fastaSeq):
        seq = ''
        fastaSeq = fastaSeq[fastaSeq.find("\n") + len("\n"):]
        for ch in fastaSeq:
            if ch != "\n":
                seq += ch
        return seq

    @staticmethod
    def areAminoAcidsSimilar(first, second):
        symbols = {'A': 0, 'C': 2, 'D': 3, 'E': 3, 'F': 1, 'G': 5, 'H': 4, 'I': 0,
                   'K': 4, 'L': 0, 'M': 0, 'N': 2, 'P': 5, 'Q': 2, 'R': 4, 'S': 2,
                   'T': 2, 'V': 0, 'W': 1, 'Y': 1, '-': 6}
        return True if symbols[first]-symbols[second] == 0 else False

    @staticmethod
    def getAminoAcidInLocationInAlignment(location: int, sequence_alignment):
        if location > len(sequence_alignment.replace("-", "")):
            return len(sequence_alignment)
        index = 0
        while location > 0:
            if 'A' <= sequence_alignment[index] <= 'Z':
                location -= 1
            index += 1
        return index

