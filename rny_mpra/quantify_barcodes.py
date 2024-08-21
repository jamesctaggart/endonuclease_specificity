from mutscan.importutils import extract_variant_sequences, quantify_barcodes
from Bio import SeqIO
import pickle
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", type=str)
    parser.add_argument("outfile", type=str)
    parser.add_argument("--umi_start",type=int, nargs='?', default=0)
    parser.add_argument("--umi_end", type=int, nargs='?', default=16)
    parser.add_argument("--vbc_start", type=int, nargs='?', default=39)
    parser.add_argument("--vbc_end", type=int, nargs='?', default=54)
    parser.add_argument("--constant_region_start", type=int, nargs='?', default=17)
    parser.add_argument("--constant_region_end", type=int, nargs='?', default=39)
    parser.add_argument("--constant_region_seq", type=str, nargs='?', default='TTGATGGTGCCTACAGGCTAGC')

    parser.add_argument('--linker', action='store_true')
    parser.add_argument('--no-linker', dest='linker', action='store_false')
    parser.set_defaults(linker=True)

    args=parser.parse_args()

    outfile = open(args.outfile,'w')

    reads_RNA = SeqIO.parse(args.fastq, "fastq")
    counts_per_VBC = quantify_barcodes(reads_RNA, umi_start=args.umi_start, umi_end=args.umi_end,
                vbc_start=args.vbc_start, vbc_end=args.vbc_end, collapse_UMI=True, has_linker=args.linker,
                constant_region_start=args.constant_region_start, constant_region_end=args.constant_region_end,
                                       constant_region_seq =args.constant_region_seq)

    pickle.dump(counts_per_VBC,outfile)  # Return counts per variant barcode


if __name__ == '__main__':
    main()