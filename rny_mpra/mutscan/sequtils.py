import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import itertools

def rc(seq):
    return str(Seq(seq).reverse_complement())


def count_mutations(seq, consensus_seq):
    """
    Determines the number of differences between a sequence seq and consensus sequence con_seq
    *** must be the same length and/or aligned
    """
    try:
        return sum(np.array(list(seq)) != np.array(list(consensus_seq)))
    except(ValueError):
        return np.nan


def get_mutations(variant_sequence, consensus_seq):
    called_mutations = ''

    if len(variant_sequence) != len(consensus_seq):
        return np.nan

    for ii, (a, b) in enumerate(zip(variant_sequence, consensus_seq)):
        if a != b:
            called_mutations += '%s-%s-%s_' % (b, ii, a)

    if called_mutations != '':
        return called_mutations[:-1]
    else:
        return called_mutations


def make_all_mutations(seq, start_idx = None, end_idx = None):
    '''
    Make a list of all possible point mutations of a sequence
    '''
    mut_list = []
    for i, c in enumerate(seq.upper()):
        for n in ['A', 'T', 'C', 'G']:
            if c != n:
                if start_idx is not None:
                    if i < start_idx:
                        continue
                if end_idx is not None:
                    if i > end_idx:
                        continue
                mut_list.append(c + '-' + str(i) + '-' + n)
    return mut_list

def make_all_double_mutations(seq, start_idx = None, end_idx = None):

    double_mut_list = []

    for i, c1 in enumerate(seq.upper()):
        for n1 in ['A', 'T', 'C', 'G']:
            if c1 != n1:
                for j, c2 in enumerate(seq.upper()):
                    if i!=j:
                        for n2 in ['A', 'T', 'C', 'G']:
                            if c2 != n2:
                                if start_idx is not None:
                                    if i<start_idx or j<start_idx:
                                        continue
                                if end_idx is not None:
                                    if i>end_idx or j>end_idx:
                                        continue
                                double_mut_list.append(c1 + '-' + str(i) + '-' + n1 + '_' + c2 + '-' + str(j) + '-' + n2)

    return double_mut_list

def find_mutations(variant_mutations, mutations_to_find, mutations_to_avoid = None):
    # Find all variants which contain the listed mutations
    tofind_set = set(mutations_to_find.split('_'))
    variant_mutation_set = set(variant_mutations.split('_'))

    if mutations_to_avoid is not None:
        try:  # Allow for input as mutation string format or as iterable that can be converted to set
            avoid_set = set(mutations_to_avoid.split('_'))
        except(AttributeError):
            avoid_set = set(mutations_to_avoid)

        avoid_set = avoid_set - tofind_set  # Don't want to avoid the things we're trying to find!

        for m in variant_mutation_set:
            if m in avoid_set:
                return False

    for m in tofind_set:
        if m not in variant_mutation_set:
            return False

    return True


def make_mutated_seq(mutation_list, reference_sequence):
    mutated_sequence = reference_sequence
    for mut in mutation_list.split('_'):
        from_nt, pos, to_nt = mut.split('-')
        pos = int(pos)

        if from_nt != mutated_sequence[pos]:
            raise ValueError('Encountered unexpected nucleotide in reference sequence.')
        mutated_sequence = replace_letter(mutated_sequence, pos, to_nt)

    return mutated_sequence

def replace_letter(string, position, newletter):
    return string[:position]+newletter+string[position+1:]


def combine_mutation_strings(ms1, ms2):
    start_nt = []
    pos = []
    end_nt = []

    for ms in (ms1, ms2):
        for m in ms.split('_'):
            s,p,e = m.split('-')
            start_nt.append(s)
            assert int(p) not in pos
            pos.append(int(p))
            end_nt.append(e)

    pos,start_nt,end_nt = zip(*sorted(zip(pos,start_nt,end_nt)))

    out = ''
    for i in range(len(pos)):
        out+='%s-%s-%s_' % (start_nt[i],str(pos[i]),end_nt[i])
    return out[:-1]

def check_mutation_string(mut_string, mutations_of_interest, acceptable_background):
    split = set(mut_string.split('_'))
    moi_split = mutations_of_interest.split('_')

    # Check to make sure all mutations of interest are in the mutation string
    for m in moi_split:
        if m not in split:
            return False
        else: # Added so if mutation of interest is not acceptable, you don't flag it later
            split.remove(m)
    
    for m in split: # Check to make sure all other mutations are in acceptable list
        if m not in acceptable_background:
            return False
        
    # If all mutations are in acceptable_background or are mutiations of interest, and we capture all ...
    # mutations of interest, return True.
    return True
        
def calculate_GC(seq):
    return (seq.count('G')+seq.count('C'))/float(len(seq))

def calculate_CtoG(seq):
    if seq.count('G') == 0:
        return len(seq)+1
    return float(seq.count('C'))/seq.count('G')

def offset_mutations(mut_string, offset):
    if len(mut_string)==0:
        return mut_string
    out_string = ''
    muts = mut_string.split('_')
    for m in muts:
        old, p, new = m.split('-')
        new_p = int(p)+offset
        out_string+='%s-%s-%s_' % (old, str(new_p), new)
    return out_string[:-1] # Trim off the last "_"