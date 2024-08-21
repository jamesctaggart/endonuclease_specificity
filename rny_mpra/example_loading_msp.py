def process_msp_data(msp, threshold, neutral_threshold, roi_start=None, roi_end=None, fold_start=None, fold_end=None):
    print('Thresholding and normalizing data...')
    msp.threshold_data(threshold, 'gDNA_count')
    msp.normalize_data()
    print('Calling mutations...')
    msp.call_mutations()
    print('Flagging non-neutral mutations...')
    msp.flag_nonneutral(column='norm_RNA:gDNA', neutral_threshold=neutral_threshold)
    print(msp.identifier, len(msp.data))
    print('Calculating double mutant dataframe...')
    msp.calculate_double_mutant_dataframe(start_idx=roi_start, end_idx=roi_end,column='norm_RNA:gDNA')
    msp.calculate_RNA_secondary_structure(start_idx=fold_start, end_idx=fold_end)


msp = MutScanPool(unmutated_cdn_seq, 0) # The second argument, stored as .zero_position, is not currently used.
msp.RNAquant_pickle_path = '/path/to/pickled/data/'
msp.DNAquant_pickle_path = '/path/to/pickled/data/'
msp.bcmap_pickle_path = '/path/to/pickled/data/'
msp.raw_path = '/path/to/raw/fastq/'
msp.identifier = 'aprE_cdn_cut_xyl_nov2021'

msp.read_from_pickle('bc_mapping_pickle','RNA_counts_pickle','gDNA_counts_pickle')
process_msp_data(msp, threshold = count_threshold, neutral_threshold=neutral_thresh, roi_start=30, roi_end=59, fold_start=20, fold_end=59). # start/end positions depend on experiment