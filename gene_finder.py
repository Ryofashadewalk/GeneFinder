# -*- coding: utf-8 -*-
"""
The genefinder function

@author: Alex Core

"""

from load import load_seq
dna = load_seq("./data/X73525.fa")
import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    '''
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    '''
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'

import doctest
doctest.run_docstring_examples(get_complement,globals(),verbose=True)



def get_reverse_complement(dna):
    """
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    a=''
    n=len(dna)-1
    while n>=0:
        a=a+get_complement(dna[n])
        n=n-1
    return a


import doctest
doctest.run_docstring_examples(get_reverse_complement,globals(),verbose=True)

def rest_of_ORF(dna):
    """
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    n=0
    a=''
    #c=''
    while n<=(len(dna)-3):
        x=dna[n]
        y=dna[n+1]
        z=dna[n+2]
        if x+y+z == 'TAG':
            n=n+3
            return a
        if x+y+z == 'TAA':
            n=n+3
            return a
        if x+y+z == 'TGA':
            n=n+3
            return a
        n=n+3
        a=a+x+y+z
    if n<= (len(dna)-1):
        a=a+dna[n]
        if n<=(len(dna)-2):
            a=a+dna[n+1]
            n=n+1
    return a
import doctest
doctest.run_docstring_examples(rest_of_ORF,globals(),verbose=True)



def find_all_ORFs_oneframe(dna):
    """
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGXYZ")
    ['ATGCATGAATGTAGA', 'ATGTGXYZ']
    >>> find_all_ORFs_oneframe("ATGGTAATGATGATGTGA")
    ['ATGGTAATGATGATG']
    """
    k=0
    c=[]
    while k<=(len(dna)-3):
        if (dna[k]+dna[k+1]+dna[k+2])=='ATG':
            a=rest_of_ORF(dna[k:])
            c.append(a)
            k=k+len(a)
        if k<=(len(dna)-3):
                k=k+3
    return c
    # TODO: implement this

doctest.run_docstring_examples(find_all_ORFs_oneframe,globals(),verbose=True)


def find_all_ORFs(dna):
    """

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    k=0
    d=[]

    while(k<3):
        #if (dna[k]+dna[k+1]+dna[k+2]) =='ATG':
        t=find_all_ORFs_oneframe(dna[k:])
        for orf in t:
            d.append(orf)
        k=k+1
    return d
    # TODO: implement this
doctest.run_docstring_examples(find_all_ORFs,globals(),verbose=True)



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    y=[]
    I=find_all_ORFs(dna)
    for orf in I:
        y.append(orf)
    rdna= get_reverse_complement(dna)
    I=find_all_ORFs(rdna)
    for orf in I:
        y.append(orf)
    return y

    # TODO: implement this

doctest.run_docstring_examples(find_all_ORFs_both_strands,globals(),verbose=True)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    n=0
    a=""
    y=[]
    I=find_all_ORFs_both_strands(dna)
    # print(I)
    while(n<=len(I)-1):
        b=I[n]
        if len(a)<len(b):
            a=b
        n=n+1
    return a
    # TODO: implement this
    pass
doctest.run_docstring_examples(longest_ORF,globals(),verbose=True)



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    n=0
    maxlength=0
    while(n<num_trials):
        newdna=shuffle_string(dna)
        b=longest_ORF(newdna)
        if maxlength<len(b):
            maxlength=len(b)
        n=n+1
    return maxlength
    # TODO: implement thiss



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    n=0
    am=''
    while n<=(len(dna)-3):
        x=dna[n]
        y=dna[n+1]
        z=dna[n+2]
        codon=x+y+z
        suck = aa_table[codon]
        am=am+suck
        n=n+3
    return am
doctest.run_docstring_examples(coding_strand_to_AA,globals(),verbose=True)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold= longest_ORF_noncoding(dna,1500)
    print(threshold)
    z=find_all_ORFs_both_strands(dna)
    final=[]
    for ORF in z:
        if(len(ORF)>threshold):
            final.append(coding_strand_to_AA(ORF))
    return final

    
print(gene_finder(dna))


if __name__ == "__main__":
    import doctest
    #print(longest_ORF("ATGCGAATGTAGCATCAAA"))
    #print(find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA"))
    # doctest.testmod()
