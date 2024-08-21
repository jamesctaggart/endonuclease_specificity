import timeit
from collections import defaultdict
import os
from Bio import SeqIO
import glob

def split_fastq(fastq_file,reads_per_file=1000000):
    line_counter = 0
    file_counter = 0
    with open(fastq_file,'r') as f:
        for line in f:
            if line_counter%(reads_per_file*4)==0: # Split into 1M read batches by default
                file_counter+=1
                outF = open('%s_split%s' % (fastq_file,str(file_counter)),'w')
            outF.write(line)
            line_counter+=1


def error_file_from_fastq(reads1, reads2, umi_start=0, umi_end=16, vbc_start=39, vbc_end=54, sampled_UMIs=None,
                              has_linker=True, constant_region_start=17, constant_region_end=39, read2_length=75,
                              collapse_UMI=True, constant_region_seq ='TTGATGGTGCCTACAGGCTAGC', error_file_dir = '', append=True, report=False,
                              reference_seq='GCTTACTTTAAAAAGCCACGCAACACGGTTCTCGTCACAGACGAAGGAGCCGCAAAGAAGTTATTAAGGGATGAATAATCCCTCAATATAAATATCTCTCACTTATTTAAAGGAGGAAACAATCATGGCAGTATAA'):
    start = timeit.default_timer()

    if sampled_UMIs == None:
        sampled_UMIs = set() # dict to check against for umi collapse
    UMI_duplicate_count = 0

    total_read1_count = 0
    total_read2_count = 0

    incorrect_linker_count = 0

    if not append:
        error_file = open(error_file_dir+'temp_error_file','w') # to reduce memory footprint
    else:
        error_file = open(error_file_dir+'temp_error_file','a')
    reference_seq = reference_seq[:read2_length] # Trim down to read length
    
    for rec1, rec2 in zip(reads1, reads2):

        total_read1_count += 1
        total_read2_count += 1

        UMI = str(rec1.seq[umi_start:umi_end])
        VBC = str(rec1.seq[vbc_start:vbc_end])

        if has_linker:
            linker_region = str(rec1.seq[constant_region_start:constant_region_end])
            if linker_region != constant_region_seq:  # Check linker1 sequence in amplicon
                incorrect_linker_count += 1
                continue  # Remove reads where linker isn't present since they aren't expected to prime in RNAseq reverse transcription
        
        if collapse_UMI:
            if UMI in sampled_UMIs:  # Only want to consider each variant once, skipping repeated UMIs
                UMI_duplicate_count += 1  # Strictly speaking, this shouldn't matter since UMIs should have same variant barcode.
                continue
            else:
                sampled_UMIs.add(UMI)

        variant_seq = str(rec2.seq[:read2_length])
        error_file.write(VBC+'\t'+variant_seq+'\n')

    assert total_read1_count == total_read2_count

    if report:
        if has_linker:
            print('%s reads removed due to linker mismatch of %s total processed' % (
            incorrect_linker_count, total_read1_count))
        if collapse_UMI:
            print('Collapsed %s reads out of %s remaining reads processed' % (
            UMI_duplicate_count, total_read1_count - incorrect_linker_count))
        else:
            print('UMI collapse not enabled.')

        stop = timeit.default_timer()
        print('Total time to generate error file: ', round(stop - start, 3), 'seconds.')

    error_file.close()
    return sampled_UMIs


def get_consensus_variant(variant_dict, min_count, multimax_threshold):

    if sum(list(variant_dict.values()))<min_count:
        return None
    
    # If all variant sequences the same, can just return the one sequence
    if len(variant_dict)==1:
        return [key for key in variant_dict][0]

    # Otherwise we take top variant, ignoring cases where the second-highest is over max value times max_threshold
    max_seq = max(variant_dict, key=lambda key: variant_dict[key])
    no_max_dict = {k:variant_dict[k] for k in variant_dict if k!=max_seq}

    if multimax_threshold*max(variant_dict.values())<max(no_max_dict.values()):  # E.g. if multimax_threshold = 0.5, second-highest variant must be less than half of max
        return None
    
    return max_seq


def read_error_file(bc_seq_file):
    error_dict = dict()
    # Read in barcode and variant seq from file
    with open(bc_seq_file,'r') as f:
        for line in f:
            fields = line[:-1].split('\t')
            VBC = fields[0]
            variant_seq = fields[1]
            if VBC in error_dict:
                error_dict[VBC][variant_seq]+=1
            else:
                error_dict[VBC] = defaultdict(int)
                error_dict[VBC][variant_seq]+=1

    return error_dict


def extract_variant_sequences(fastq1, fastq2, multimax_threshold, min_count, umi_start=0, umi_end=16, vbc_start=39, vbc_end=54,
                              has_linker=True, constant_region_start=17, constant_region_end=39, read2_length=75,
                              collapse_UMI=True, constant_region_seq ='TTGATGGTGCCTACAGGCTAGC', error_file_dir = '/Volumes/T7/endonuclease_site_mutational_scanning_data/', report=True,
                              reference_seq='GCTTACTTTAAAAAGCCACGCAACACGGTTCTCGTCACAGACGAAGGAGCCGCAAAGAAGTTATTAAGGGATGAATAATCCCTCAATATAAATATCTCTCACTTATTTAAAGGAGGAAACAATCATGGCAGTATAA'):
    '''
    Extract the variant barcode (VBC) and gDNA sequence from raw fastq data
    If collapse_UMI is True, will only take the first instance of each UMI
    If has_linker=True, require presence of linker sequence

    reads: generator producing reads, corresponding to read1 in paired end
    returns: OrderedDict of barcodes
    '''

    # Split input fastq
    split_fastq(fastq1,reads_per_file=1000000)
    split_fastq(fastq2,reads_per_file=1000000)

    split_files1 = sorted(glob.glob(fastq1+'_split*'))
    split_files2 = sorted(glob.glob(fastq2+'_split*'))
    # # pdb.set_trace()

    sampled_UMIs = set()

    os.remove(error_file_dir+'temp_error_file') # Clear temporary error file, if it exists.

    for f1,f2 in zip(split_files1,split_files2):
        # Read in split fastqs  
        reads1 = SeqIO.parse(f1, "fastq")
        reads2 = SeqIO.parse(f2, "fastq")
    
        # Generate file which contains VBC and measured sequence for each read (tab delimited)
        sampled_UMIs = error_file_from_fastq(reads1, reads2, umi_start=umi_start, umi_end=umi_end, vbc_start=vbc_start, vbc_end=vbc_end,
                                has_linker=has_linker, constant_region_start=constant_region_start, constant_region_end=constant_region_end, read2_length=read2_length,
                                collapse_UMI=collapse_UMI, constant_region_seq =constant_region_seq, error_file_dir = error_file_dir,
                                reference_seq=reference_seq,append=True,report=report,sampled_UMIs=sampled_UMIs)

    # Remove split fastq files
    for f in split_files1:
        os.remove(f)
    for f in split_files2:
        os.remove(f)

    # Collapse barcodes down to the top variant
    reconcile_errors = read_error_file(error_file_dir+'temp_error_file')

    variant_seqs = dict() # initialize final dict to return
    multi_max_barcodes = []
    for VBC in reconcile_errors:
        max_seq = get_consensus_variant(reconcile_errors[VBC],multimax_threshold=multimax_threshold, min_count=min_count)
        if max_seq is None:
            multi_max_barcodes.append(VBC)
        else:
            variant_seqs[VBC] = max_seq

    print('Dropped %s of %s barcodes due to multimax' % (str(len(multi_max_barcodes)),str(len(variant_seqs)+len(multi_max_barcodes))))

    # os.remove(error_file_dir+'temp_error_file')
    return variant_seqs

def quantify_barcodes(reads, umi_start=0, umi_end=16, vbc_start=39, vbc_end=54, collapse_UMI=True, has_linker=True, constant_region_start=17, constant_region_end=39, constant_region_seq ='TTGATGGTGCCTACAGGCTAGC'):
    '''

    Extract the variant barcode (VBC) and quantify
    If collapse_UMI is True, will only take the first instance of each UMI

    reads: generator producing reads, corresponding to read in RNAseq data
    VBC2seq: dictionary which takes variant barcode and returns the gDNA sequence of the variant
    returns: Dictionary with counts per barcode
    '''
    start = timeit.default_timer()

    sampled_UMIs = set()

    UMI_duplicate_count = 0
    incorrect_linker_count = 0
    total_read_count = 0
    variant_counts = defaultdict(int)

    for read in reads:
        total_read_count += 1
        UMI = str(read.seq[umi_start:umi_end])
        VBC = str(read.seq[vbc_start:vbc_end])

        if has_linker:

            linker_region = str(read.seq[constant_region_start:constant_region_end])
            if linker_region != constant_region_seq:  # Check linker1 sequence in amplicon
                incorrect_linker_count += 1
                continue  # Remove reads where linker isn't present since they aren't expected to prime in RNAseq reverse transcription

        if collapse_UMI:
            if UMI in sampled_UMIs:  # Only want to consider each variant once, skipping repeated UMIs
                UMI_duplicate_count += 1  # Strictly speaking, this shouldn't matter since UMIs should have same variant barcode.
                continue
            else:
                sampled_UMIs.add(UMI)

        variant_counts[VBC] += 1

    if has_linker:
        print('%s reads removed due to linker mismatch of %s total processed' % (incorrect_linker_count, total_read_count))
    if collapse_UMI:
        print('Collapsed %s reads out of %s remaining reads processed' % (
        UMI_duplicate_count, total_read_count - incorrect_linker_count))
    stop = timeit.default_timer()
    print('Total time running: ', round(stop - start, 3), 'seconds.')

    return variant_counts


def map_umi_to_vbc(reads, umi_start=0, umi_end=16, vbc_start=39, vbc_end=54,
                   has_linker=True, linker_start=17, linker_end=39, linker_seq='TTGATGGTGCCTACAGGCTAGC'):
    '''
    Extract umi and variant barcodes from fastq-formated reads and generate
    dictionary mapping umi to variant barcode


    UMI = Unique molecular index
    VBS = Variant barcode sequence
    Used to check how many instances of each variant barcode sequence there are
    Creates a dictionary with keys: UMI, values: VBS

    :param reads1: first read of paired end sequencing, read from fastq with SeqIO.parse

    '''
    start = timeit.default_timer()

    total_read1_count = 0
    incorrect_linker_count = 0
    UMI2VBC = dict()

    for rec in reads:

        total_read1_count += 1

        UMI = str(rec.seq[umi_start:umi_end])  # Unique Molecular Index sequence
        VBC = str(rec.seq[vbc_start:vbc_end])  # Variant BarCode sequence

        if has_linker:  # Require perfect match to linker if specified True
            linker_region = str(rec.seq[linker_start:linker_end])
            if linker_region != linker_seq:
                incorrect_linker_count += 1
                continue

        UMI2VBC[UMI] = VBC

    print('%s reads removed due to linker mismatch of %s total processed' % (incorrect_linker_count, total_read1_count))

    stop = timeit.default_timer()
    print('Total time running: ', round(stop - start, 3), 'seconds.')

    return UMI2VBC