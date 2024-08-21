import numpy as np

comp = { 'A':'T','T':'A',
         'G':'C','C':'G',
         'N':'N'}

def rev_comp(seq):
    rc = [ comp[N] for N in reversed(seq) ]
    return "".join(rc)

def per_position_AT(seqs):
    seq_length = len(seqs[0])
    AT_count = np.zeros(seq_length)
    for seq in seqs:
        if len(seq)!=seq_length:
            raise RuntimeError('Sequences not the same length')
        for i,c in enumerate(seq):
            if c in ('A','T','a','t','u','U'):
                AT_count[i]+=1
    return AT_count / len(seqs)

def per_position_GC(seqs):
    return 1-per_position_AT(seqs)

def per_position_G(seqs):
    seq_length = len(seqs[0])
    G_count = np.zeros(seq_length)
    for seq in seqs:
        assert len(seq)==seq_length
        for i,c in enumerate(seq):
            if c in ('G','g'):
                G_count[i]+=1
    return G_count / len(seqs)

def per_position_purine(seqs):
    seq_length = len(seqs[0])
    purine_count = np.zeros(seq_length)
    for seq in seqs:
        assert len(seq)==seq_length
        for i,c in enumerate(seq):
            if c in ('A','G','a','g'):
                purine_count[i]+=1
    return purine_count / len(seqs)

def per_position_pyrimidine(seqs):
    return 1-per_position_purine(seqs)

def count_flanking_motif(seqs, motif='GA'):
    assert len(seqs[0])==40
    counter = 0
    for seq in seqs:
        if seq[19:21] == motif:
            counter +=1
    return counter/len(seqs)