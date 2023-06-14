from pathlib import Path
import gzip
import typing
from textwrap import wrap

ExtendedIUPACProtein = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
IUPACProtein = 'ACDEFGHIKLMNPQRSTVWY'
IUPACAmbiguousDNA = 'GATCRYWSMKHBVDN'
IUPACUnambiguousDNA = 'GATC'
ExtendedIUPACDNA = 'GATCBDSW'
IUPACAmbiguousRNA = 'GAUCRYWSMKHBVDN'
IUPACUnambiguousRNA = 'GAUC'
TRANSLATION_TABLE = {
    1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    3: 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    5: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
    6: 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    9: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    10: 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    11: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    12: 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    13: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
    14: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    15: 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    16: 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    21: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    22: 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    23: 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    24: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
    25: 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    26: 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    27: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    28: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    29: 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    30: 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    31: 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    32: 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    33: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'
}


class Fasta:
    def __init__(self, path: Path):
        """
        :param fasta_string: ">this_is_a_header this is a description\nthis\nis\na\sequence\n"
        :return: object of type Record
        """
        if not path.is_file():
            raise f"{path} is not a valid file"
        self.path = path.absolute()
        self.gzipped = True if '.gz' in path.suffixes else False
        self.records = []

    def parse(self):
        if self.gzipped:
            with gzip.open(self.path, 'rt') as fh:
                self.records += [Record(i) for i in fh.read().split('>') if i]
        else:
            with open(self.path, 'rt') as fh:
                self.records += [Record(i) for i in fh.read().split('>') if i]


class Record:
    def __init__(self, fasta_string: str):
        """
        :param fasta_string: A fasta string that has been split by '>'
        :return: object of type Record
        """
        fasta_string = fasta_string.split('\n')
        self.header = fasta_string[0].replace('>', '').strip()
        self.accession, self.description = self.header.split(' ', 1)
        self.seq = Seq(''.join(fasta_string[1:]))

    def __len__(self):
        return self.seq.length

    def __repr__(self):
        return self.accession

    def __getitem__(self, key: slice):
        return Record(f"{self.header}\n{self.seq.seq[key]}")

    def as_fasta(self, description=False, flatten=False, wrap_width=80):
        header = f"{self.accession} {self.description}" if description else self.accession
        return f">{header}\n{self.seq if flatten else self.seq.seq_wrap(wrap_width)}"


class Seq:
    def __init__(self, fasta_lines: str, gencode: int = 11,
                 mol: typing.Literal['dna', 'rna', 'protein', 'unknown'] = 'unknown'):
        """
        :param fasta_lines: "this\nis\na\sequence\n"
        :return: object of type Seq
        """
        self.seq = fasta_lines.strip()
        self.mol = mol
        self.gencode = gencode
        self.length = len(self.seq)

    def __len__(self):
        return self.length

    def __repr__(self):
        return self.seq

    def __getitem__(self, key: slice):
        return Seq(fasta_lines=self.seq[key], gencode=self.gencode, mol=self.mol)

    def determine_mol(self):
        """
        Cheeky fun to determine molecule with Jaccard similarity on a set of the sequence
        :return: str of either 'dna', 'rna' or 'protein'
        """
        seq = set(self.seq.upper())

        return [
            k for k, v in sorted(
                {
                    'dna': jaccard_similarity(seq, IUPACUnambiguousDNA),
                    'rna': jaccard_similarity(seq, IUPACUnambiguousRNA),
                    'protein': jaccard_similarity(seq, IUPACProtein)
                }.items(),
                key=lambda item: item[1], reverse=True
            )
        ][0]

    def translate(self, to_stop=False):
        """
        Cheeky fun to determine molecule with Jaccard similarity on a set of the sequence
        :return: Object of type Seq
        """
        if self.mol == 'unknown':
            self.mol = self.determine_mol()
        if self.mol == 'dna':
            table = dict(zip(
                [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],
                TRANSLATION_TABLE[self.gencode])
            )
            protein = ''
            for i in range(0, self.length, 3):
                aa = table.get(self.seq[i: i + 3], '*')
                protein += aa
                if to_stop and aa == '*':
                    return Seq(protein, mol='protein', gencode=self.gencode)
            return Seq(protein, mol='protein', gencode=self.gencode)

    def seq_wrap(self, width: int = 80) -> str:
        """
        Wraps sequence to a specified width
        """
        return '\n'.join(wrap(self.seq, width))

    def seq_flatten(self) -> str:
        """
        Removes newlines from sequence
        """
        return self.seq.replace('\n', '')


def jaccard_similarity(a, b) -> float:
    """
    Cheeky fun to determine molecule with Jaccard similarity on a set of the sequence
    :return: str of either 'dna', 'rna' or 'protein'
    """
    return len(a.intersection(b)) / len(a.union(b))


def jaccard_distance(a: set, b: set) -> float:
    """
    Cheeky fun to determine molecule with Jaccard similarity on a set of the sequence
    :return: str of either 'dna', 'rna' or 'protein'
    """
    return len(a.symmetric_difference(b)) / len(a.union(b))
