import subprocess
import numpy as np
import time
import pdb
from rend_utilities.wig_utils import check_nearby_peak
from endmap.calc_endseq_columns import get_sequences_from_position
from Bio import SeqIO

BS168_genome = SeqIO.to_dict(SeqIO.parse("reference_files/B_subtilis_168_NC_000964.3.fa", "fasta"))['NC_000964.3'].seq

def make_kmers(seq,k,write_interval = 1,output_file = None):
    '''
    Given a string-formatted sequence, return all overlapping kmers
    :param seq: str, Input sequence to be broken into kmer
    :param k: int, length of kmers to produce
    :param output_file: str, file handle for kmer output
    :return: list, all kmers in order
    '''
    kmers = []
    for ii in range(len(seq)-k):
        kmers.append(seq[ii:ii+k])

    counter = 0
    if output_file is not None:
        with open(output_file,'w') as f:
            for kmer in kmers:
                if counter % write_interval == 0:
                    f.write(str(kmer)+'\n')
                counter +=1

    return kmers


def fold_all(sequences, constraints=None):
    '''
    :param sequences: File containing sequences, each separated by line break
    :param constraints: File containing constraints, each separated by line break
    :return:
    '''

    proc = subprocess.Popen(['RNAfold','-i',sequences],stdout=subprocess.PIPE)

    deltaG_array = []

    RNAFold_output = proc.stdout.readlines()
    ii = 0
    for line in RNAFold_output:
        if isinstance(line, bytes):
            line = line.decode()
        if line[0] in ['A','U','C','G']:
            seq = line.rstrip()
        elif line[0]=='.' or line[0]=='(':
            fields = line.split()
            if fields[1]=='(':
                deltaG = fields[2][:-1]
            else:
                deltaG = fields[1][1:-1]

            deltaG_array.append(deltaG)
            ii += 1

    return deltaG_array


def metastructure(forward_positions, reverse_positions, dG_array_f, dG_array_r, halfrange, ref_peaks_for, ref_peaks_rev, return_seqs = False):
    '''

    If site is within distance windowsize of the end of their transcript, we discard the peak.

    ga = GenomeArray containing average MFE per position within a sliding window
    halfrange = half-length of the range to be considered around peak (NOT the averaging window for MFE calculation)
    ref_peaks = positions of reference peaks which define transcript boundaries
    '''
    skip_count = 0
    cleavage_sites = np.zeros((len(forward_positions) + len(reverse_positions), halfrange * 2))
    ii = 0
    sequences = []
    for site in forward_positions:
        if check_nearby_peak(site, ref_peaks_for,
                             halfrange):  # Check if there is a peak nearby the cleavage site being considered
            skip_count += 1
            # pdb.set_trace()
            cleavage_sites[ii, :] = np.nan
            ii += 1
            
            continue  # If there is a peak in the WT sample near this position, skip it. Don't want to fold things next to TSS or intrinsic terminators.

        cleavage_sites[ii] = dG_array_f[site - halfrange: site + halfrange]
        sequences.append(get_sequences_from_position(site, '+', 75, BS168_genome))
        ii += 1

    for site in reverse_positions:
        if check_nearby_peak(site, ref_peaks_rev, halfrange):
            skip_count += 1
            cleavage_sites[ii, :] = np.nan
            ii += 1

            continue
        cleavage_sites[ii] = dG_array_r[site - halfrange: site + halfrange][::-1]
        sequences.append(get_sequences_from_position(site, '-', 75, BS168_genome))
        ii += 1

    print('Skipped %s out of %s sites.' % (skip_count, ii))
    if return_seqs:
        return np.nanmean(cleavage_sites, axis=0), sequences
    else:
        return np.nanmean(cleavage_sites, axis=0)
    


def get_all_mfe_structures(sequences):   
    '''
    For use in extracting structural properties of RNase III substrates
    '''

    proc = subprocess.Popen(['RNAfold','-i',sequences],stdout=subprocess.PIPE)


    deltaG_array = []
    dotbracket_list = []
    
    RNAFold_output = proc.stdout.readlines()
    ii = 0
    for line in RNAFold_output:
        if isinstance(line, bytes):
            line = line.decode()

        if line[0] in ['A','U','C','G']:
            seq = line.rstrip()
        elif line[0]=='.' or line[0]=='(':
            fields = line.split()
            dotbracket_list.append(fields[0])
            if fields[1]=='(':
                deltaG = fields[2][:-1]
            else:
                deltaG = fields[1][1:-1]

            deltaG_array.append(deltaG)
            ii += 1

    return dotbracket_list, deltaG_array


def get_structural_distance(p1, p2, dotbracket):
    traversed_structure = dotbracket[p1:p2]
    d = 0 # d is the structural distance counting only base pairs
    overhang = False
    unpaired_5prime = 0
    oh_l = 0

    for i,c in enumerate(traversed_structure):

        if c=='(':
            d+=1
        elif c==')':
            d-=1
        elif c=='.' and d==0 and i<6:  # Assuming that it won't form a full stem-loop in first 6 nt 
            print('Unpaired 5prime!')
            unpaired_5prime+=1                # and also won't have more than 6 unpaired
        elif c=='.' and d==0 and i==6:
            raise RuntimeError('6 nt isnt enough for unpaired threshold!')

        if overhang:
            oh_l+=1

        if d-unpaired_5prime==0 and i>10:  # i 10 means you won't call "overhang" in first 10 iterations, so very short stems at start will be ignored
            overhang = True
        if overhang and d>0:
            print('Warning - sites appear to be in different stems!')
            return np.nan, np.nan

    if not overhang:
        print('Warning - no 3\' overhang!')
        return np.nan, np.nan
    else:
        return oh_l, d
    

def get_first_stem_end(dotbracket):
    d = 0 # d is the structural distance counting only base pairs
    in_stem = False
    for i,c in enumerate(dotbracket):
        if c=='(':
            d+=1
            in_stem=True
        elif c==')':
            d-=1
        elif c=='.' and d==0 and in_stem:  # Assuming that it won't form a full stem-loop in first 6 nt 
            return i
    return len(dotbracket)


def get_positionwise_paired_fraction(dotbrackets, seqs, positions):
    paired_array = np.zeros(20)
    seq_list = []

    paired_nt_arrays = {'A':np.zeros(20),'T': np.zeros(20),'C':np.zeros(20),'G':np.zeros(20)}

    for j,db in enumerate(dotbrackets):
        p = positions[j]-10 # starting position
        in_duplex=True
        for i in range(20):
            if not in_duplex:
                paired_array[i] += 0.0
            elif db[p+i] == '(':
                paired_array[i] += 1.0
                paired_nt_arrays[seqs[j][p+i]][i] += 1.0
            elif db[p+i] == '.':
                paired_array[i] += 0.0
            elif db[p+i] == ')':
                in_duplex = False


    for key,arr in paired_nt_arrays.items():
        paired_nt_arrays[key] = arr/paired_array  # Scale values by the total number of pairs at each position


    return paired_array/len(dotbrackets), paired_nt_arrays
            
