from mutscan.importutils import extract_variant_sequences, quantify_barcodes
from Bio import SeqIO
import argparse
import pickle

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_read1", type=str)
    parser.add_argument("fastq_read2", type=str)
    parser.add_argument("outfile", type=str)
    parser.add_argument("--umi_start",type=int, nargs='?', default=0)
    parser.add_argument("--umi_end", type=int, nargs='?', default=16)
    parser.add_argument("--vbc_start", type=int, nargs='?', default=39)
    parser.add_argument("--vbc_end", type=int, nargs='?', default=54)
    parser.add_argument("--read2_length", type=int, nargs='?', default=75)
    parser.add_argument("--constant_region_start", type=int, nargs='?', default=17)
    parser.add_argument("--constant_region_end", type=int, nargs='?', default=39)
    parser.add_argument("--constant_region_seq", type=str, nargs='?', default='TTGATGGTGCCTACAGGCTAGC')

    args = parser.parse_args()

    # reads1 = SeqIO.parse(args.fastq_read1, "fastq")
    # reads2 = SeqIO.parse(args.fastq_read2, "fastq")

    outfile = open(args.outfile, 'wb')

    bc_to_seq = extract_variant_sequences(args.fastq_read1, args.fastq_read2, multimax_threshold=0.25, min_count=3, \
                                          read2_length = args.read2_length, error_file_dir='/Users/james/Dropbox (HMS)/T7/endonuclease_site_mutational_scanning_data/pickles/') # Add other args here if changing

    pickle.dump(bc_to_seq,outfile)  # Return counts per variant barcode


if __name__ == '__main__':
    main()