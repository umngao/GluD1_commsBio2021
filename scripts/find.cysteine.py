#!/usr/bin/env python

from Bio import SeqIO, Seq
import numpy as np


def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

fa='/homes/lianggao/HMW_glutenin/data/210605.all.glud1x.glud1y.aa.cmb.fa'
fasta_sequences=SeqIO.parse(open (fa), 'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    vec = findOccurrences(sequence, 'C')
    vec = np.asarray(vec)
    string_ints = [str(int) for int in vec]
    print  (name, '\t'.join(string_ints))


