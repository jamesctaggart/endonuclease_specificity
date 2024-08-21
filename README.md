# Analysis of stabilized end sequencing and endoribonuclease cleavage MRPA

Scripts used for identification of endoribonuclease cleavage positions from stabilized end sequencing data and analysis of endonuclease cleavage MPRA data.

## Analysis of stabilized end sequencing data
Included are scripts used to go from Rend-seq, 5'-end, or 3'-end sequencing wig files to called positions of cleavage. These scripts are contained within endmap, with additional helper functions contained within rend_utilities. Examples for how these scripts are used are contained within the three example python scripts:
1. example_initialize_endmaps.py demonstrates how to convert .wig files into EndMap objects, which are described in endmap/datastructures.py and used in downstream steps.
2. example_process_endmaps.py shows how to integrate EndMap objects generated for each genetic background and sequencing data type of interest into a single DataFrame which can be used to call positions of cleavage.
3. example_call_sites.py shows how to take this DataFrame, calculate sensitivity values for each position, and extract positions of cleavage.

Additionally provided is an annotation file specific to Bacillus subtilis that is used when reading in .wig files from this species, CDS_168.mochiview_v13.txt.

## Analysis of MPRA (rny_mpra)
Included are scripts used to process and visualize data from our MPRA to study the effects of mutation on the processing of an RNase Y substrate.

Each experiment involves three sequencing libraries, two for gDNA/RNA barcode quantification and one for mapping of barcodes to variant sequences. Raw fastq files from these libraries are processed using extract_variant_sequences.py (for barcode to variant sequence mapping) or quantify_barcodes.py (for gDNA/RNA barcode quantification). These scripts generate pickles which can be loaded and processed using functions in mutscan to generate a MutScanPool object, which summarizes the results a particular experiment.

mutscan contains functions and classes used in analysis of data derived from our mpra, including the MutScanPool object used to represent each experiment. Included are functions used to convert the previously-described pickled data into MutScanPool objects, summarized in example_loading_msp.py. mutscan additionally contains all functions used in the downstream analysis/visualization of these data, examples of which can be seen in example_plotting.py.

