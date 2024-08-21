import pandas as pd
import pdb
import pickle
from collections import Counter

def write_mochi_location_set(site_df, out_file, sequence_set = 'NC_000964.3'):
    f = open(out_file,'w')
    f.write('SEQ_NAME\tSTART\tEND\tSTRAND\n')
    for row in site_df.iterrows():
        fields = row[1]
        if fields['strand'] == '+':
            f.write('%s\t%s\t%s\t%s\n' % (sequence_set, str(fields['position']-1),str(fields['position']),fields['strand']))
        elif fields['strand'] == '-':
            f.write('%s\t%s\t%s\t%s\n' % (sequence_set, str(fields['position']),str(fields['position']+1),fields['strand']))
        else:
            raise ValueError('Strand must be + or -')
    f.close()

def subsample_wig(wigF, wigR):
    read_dict = {}

    with open(wigF,'r') as f:
        for line in f:
            if line.find('=')!=-1:
                continue
            fields = line[:-1].split('\t')
            p = fields[0]
            rc = fields[1]
            read_dict[int(p)] = float(rc)

    with open(wigR,'r') as f:
        for line in f:
            if line.find('=')!=-1:
                continue
            fields = line[:-1].split('\t')
            p = fields[0]
            rc = fields[1]
            read_dict[int('-'+p)] = float(rc)

    read_counter = Counter(read_dict)

    return read_counter

def wigs_from_counter(counter, wigF_out, wigR_out, header=None):
    items = sorted(counter.items())

    outF = open(wigF_out, 'w')
    outR = open(wigR_out, 'w')

    if header is not None:
        outF.write(header)
        outR.write(header)

    fstring = ''
    rstring = ''

    for item in items:
        if item[0]<0:
            rstring = '%s\t%s\n' % (str(abs(item[0])),str(item[1]))+rstring
        elif item[0]>0:
            fstring+='%s\t%s\n' % (str(abs(item[0])),str(item[1]))
    
    outF.write(fstring)
    outR.write(rstring)

    outF.close()
    outR.close()