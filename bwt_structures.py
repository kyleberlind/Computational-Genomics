from collections import Counter
from compsci260lib import *


# Note that the '$' character will be used to designate the end of a given string.

def forward_bwt(seq):
    """
    forward_bwt(seq) takes as input a string containing the EOF
    character to which the BWT must be applied. The method should
    then return the result of the BWT on the input string.
    
    For example:
    print forward_bwt('GATTACA$') --> 'ACTGA$TA'
    """
    
    string_array = build_bwt_array(seq)
    bwt_string = ""
    for string in string_array:
        bwt_string += string [-1]
    #print "forward: ", bwt_string
    return bwt_string
    
    
def build_bwt_array(seq):
    """helper function which will build the BWT matix out of arrays 
        by rotating the fits letter of the sting to the back of the
        string and appending the permutation to a list"""
    
    string_array = []
    for i in range(len(seq)):
        seq = seq[1:] + seq[0]
        string_array.append(seq)

    return sorted(string_array)
    
    
def reverse_bwt(seq):
    """
    reverse_bwt(seq) takes as input a string containing the EOF
    character to which the reverse of the BWT must be applied.
    The method should then return the result of the reversal on
    the input string.
    
    For example:
    print reverse_bwt('ACTGA$TA') --> 'GATTACA$'
    
    this thing runs in like 2n^2 time
    """    
    
    orig_seq = list(seq)
    orig_seq_sorted = sorted(orig_seq)
    BWT_matrix = []
    i = 0
    while i < len(orig_seq)-1:
        for j in range(len(orig_seq)):
            if BWT_matrix == []:
                [BWT_matrix.append([x]) for x in orig_seq_sorted]   
            BWT_matrix[j] = [orig_seq[j]] + BWT_matrix[j]
        BWT_matrix = sorted(BWT_matrix)   
        i+=1
    for element in  BWT_matrix:
        if element[-1]== "$":
            reverse_BWT = "".join(element)   
    return reverse_BWT
            
            
    


# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.

def make_suffix_array(seq):
    """
    Returns the suffix array of the input string.
    
    For example:
    print make_suffix_array('GATTACA$') --> [7, 6, 4, 1, 5, 0, 3, 2]
    """
    suffixes = {}
    for x in range(len(seq)):
        suffixes[seq[x:]] = x
    suffix_array = [suffixes[suffix] for suffix in sorted(suffixes.keys())]
    return suffix_array

def rank(bwt_seq):
    """
    Takes as input a string transformed by the BWT. Returns a
    dictionary with characters as keys and lists as values.
    Each list contains the total number of occurrences for the 
    corresponding character up until each position in the
    BWT-transformed string (i.e., its rank).
    
    For example:
    print rank('ACTGA$TA')['$'] --> [0, 0, 0, 0, 0, 1, 1, 1]
    print rank('ACTGA$TA')['A'] --> [1, 1, 1, 1, 2, 2, 2, 3]
    print rank('ACTGA$TA')['C'] --> [0, 1, 1, 1, 1, 1, 1, 1]
    print rank('ACTGA$TA')['G'] --> [0, 0, 0, 1, 1, 1, 1, 1]
    print rank('ACTGA$TA')['T'] --> [0, 0, 1, 1, 1, 1, 2, 2]
    """
    rank = {}
    characters = set(bwt_seq)
    for character in characters:
        rank[character] = [0]
    rank[bwt_seq[0]] = [1]
    for letter in bwt_seq[1:]:
        for k, v in rank.items():
            v.append(v[-1] + (k == letter))
    return rank

def count_smaller_chars(seq):
    """
    Takes as input a string. Returns a dictionary with characters
    as keys and integers as values. The integers track the
    number of characters in the input string which are
    lexicographically smaller than the corresponding character
    key.
    
    For example, using an input DNA sequence like 'GATTACA':
    print count_smaller_chars('GATTACA')['A'] --> 0  (A, being lexicographically first in a DNA sequence, should always return 0)
    print count_smaller_chars('GATTACA')['C'] --> 3  (C, being second, should return the number of A's, which here is 3)
    print count_smaller_chars('GATTACA')['G'] --> 4  (G, being third, should return the number of A's or C's, which here is 4)
    print count_smaller_chars('GATTACA')['T'] --> 5  (T, being fourth, should return the number of A's or C's or G's, which here is 5)
    """
    characters = set(seq)
    cntr = Counter(seq)
    total = 0
    counts = {}
    for character in sorted(characters):
        counts[character] = total
        total += cntr[character]
    return counts

def make_all(reference):
    """
    Takes as input a reference string. Returns the data
    structures necessary to perform efficient exact string
    matching searches.
    """
    counts = count_smaller_chars(reference)
    reference = reference + '$'
    suffix_array = make_suffix_array(reference)
    bwt = forward_bwt(reference)
    ranks = rank(bwt)
    return bwt, suffix_array, ranks, counts

if __name__ == '__main__':
    # example sequences for forward_bwt() and reverse_bwt()
    # you may change them to different sequences to test your code.
    seq1 = 'GATTACA$'
    seq2 = 'ACTGA$TA'
    
    forward_bwt(seq1)
    reverse_bwt(seq2)
