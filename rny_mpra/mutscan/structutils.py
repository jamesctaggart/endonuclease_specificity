import subprocess
import pdb

def pair_type(n1,n2):
    if (n1 == 'G' and n2 == 'C') or (n1 == 'C' and n2 == 'G'):
        return 'GC'
    if (n1 == 'G' and n2 == 'T') or (n1 == 'T' and n2 == 'G'):
        return 'w'
    if (n1 == 'A' and n2 == 'T') or (n1 == 'T' and n2 == 'A'):
        return 'AT'
    if (n1 == 'A' and n2 == 'G') or (n1 == 'G' and n2 == 'A') or \
            (n1 == 'A' and n2 == 'C') or (n1 == 'C' and n2 == 'A') or \
            (n1 == 'T' and n2 == 'C') or (n1 == 'C' and n2 == 'T'):
        return 'u'
    

def fold_all(sequences):
    '''
    :param sequences: File containing sequences, each separated by line break
    '''

    proc = subprocess.Popen(['RNAfold','-i',sequences],stdout=subprocess.PIPE)

    deltaG_array = []
    fold_array = []
    RNAFold_output = proc.stdout.readlines()
    ii = 0
    for line in RNAFold_output:
        if isinstance(line, bytes):
            line = line.decode()
        if line[0] in ['A','U','C','G']:
            seq = line.rstrip()
        elif line[0]=='.' or line[0]=='(':
            fields = line.split()
            fold = fields[0]
            if fields[1]=='(':
                deltaG = fields[2][:-1]
            else:
                deltaG = fields[1][1:-1]

            deltaG_array.append(float(deltaG))
            fold_array.append(fold)

            ii += 1

    return deltaG_array, fold_array