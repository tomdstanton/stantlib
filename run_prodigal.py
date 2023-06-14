#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_prodigal
-----------------
Function to run and parse the Prodigal gene prediction command
line tool (Hyatt and Houser).
Included the CDS class to store predicted genes and objects.
-----------------
Requires:
- Python >=3.7
-----------------
Tom Stanton, 2021
"""

# ------------ Functions ------------#
def run_prodigal(fasta_string):
    cmd = ['prodigal', '-a', '/dev/stdout', '-f', 'gff']
    child = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    child.stdin.write(fasta_string.encode())
    out, err = child.communicate()
    if err:
        print(err.decode("utf-8"))
    if out:
        return out.decode("utf-8").strip()


def parse_prodigal(seq_dict):
    fasta_string = '\n'.join([f'>{k} {v.description}\n{str(v.seq)}' for k, v in record_dict.items()])
    gff, proteins = run_prodigal(fasta_string).split('>', 1)
    protein_lines = proteins.split('>')
    gff_lines = gff.split('\n')[3:-1]
    return [CDS(seq_dict, gff_lines[i], protein_lines[i]) for i in range(len(protein_lines))]


# ------------ Classes ------------#
class CDS(object):
    """Creates CDS objects from Prodigal fasta.faa and gff"""

    def __init__(self, seq_dict, gff_line, protein_line):
        self.seqname, self.source, self.feature, self.start, \
        self.end, self.score, self.strand, self.frame, \
        attr = gff_line.split('\t')
        self.ID, self.partial, self.start_type, self.rbs_motif, \
        self.rbs_spacer, self.gc_cont, self.conf, self.score, \
        self.cscore, self.sscore, self.rscore, self.uscore, \
        self.tscore = re.findall('(?<==).*?(?=;)', attr)

        self.protein_seq = Seq(''.join(protein_line.split('\n')[1:]))
        if self.strand == '+':
            self.feature = SeqFeature(FeatureLocation(int(self.start) - 1, int(self.end)), type=self.feature, strand=+1)
            self.nucleotide_seq = seq_dict[self.seqname].seq[self.feature.location.start:self.feature.location.end]
        else:
            self.feature = SeqFeature(FeatureLocation(int(self.start) - 1, int(self.end)), type=self.feature, strand=-1)
            self.nucleotide_seq = seq_dict[self.seqname].seq[
                                  self.feature.location.start:self.feature.location.end].reverse_complement()

        aa_codon_pairs = [str(self.protein_seq)[i:i + 2]
                          for i in range(len(self.protein_seq))
                          if len(str(self.protein_seq)[i:i + 2]) > 1]

        dna_codon_pairs = [str(self.nucleotide_seq)[i:i + 6]
                           for i in range(0, len(self.nucleotide_seq), 3)
                           if len(str(self.nucleotide_seq)[i:i + 6]) > 5]

        self.codon_pairs = {i: (dna_codon_pairs[i], aa_codon_pairs[i]) for i in range(len(dna_codon_pairs))}


