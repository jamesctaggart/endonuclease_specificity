

def check_nearby_peaks(peak_position, reference_set, threshold):
    '''
    returns True if position "peak_position" is within distance "threshold" of any peak in set of peaks ("reference_set")
    :param peak_position: int, position of interest
    :param reference_set: set, contains int-formatted positions peak_position is to be checked against
    :param threshold: int, distance threshold by which less than or equal should return True
    :return: bool, True if nearby peak is present in reference_set, else False
    '''
    for peak in reference_set:
        if abs(float(peak_position)-float(peak))<=threshold:
            return True
    return False


def parse_endenrich_file(end_enrichment_fh):  # Currently not used?
    '''
    Parse end enrichment files generated by Rend_seq_end_enrichment_File_header_20180206.m (JBL)
    :param end_enrichment_fh:  file name of end_enrichment file
    :return: dictionary containing peak positions, with keys ['3f', '3r', '5f', '5r']
    '''

    peaks = {'3f':[],'3r':[],'5f':[],'5r':[]}

    with open(end_enrichment_fh,'r') as f:
        for line in f:
            fields = line.split()
            if fields[0] == '1':
                peaks['3f'].append(int(fields[1]))
            elif fields[0] == '2':
                peaks['3r'].append(int(fields[1]))
            elif fields[0] == '3':
                peaks['5f'].append(int(fields[1]))
            elif fields[0] == '4':
                peaks['5r'].append(int(fields[1]))
            else:
                raise ValueError('Unknown peak type in end_enrichment file')

    return peaks


def parse_zscores(zscore_peak_file, mode):

    '''
    Given a table generated by call_peak_z_score.py, parse either enriched or depleted values by z-score.
    :param zscore_peak_file:
    :param mode: tells function to either detect peaks with positive or negative Z-scores
    :return: Two sets, corresponding to the enriched (or depleted) peaks on the plus and minus strands.
    '''

    if mode in ('negZ', 'depleted', 'negative'):  # Accept multiple values inputs for mode
        mode = 'negZ'
    elif mode in ('posZ', 'enriched', 'positive'):
        mode = 'posZ'

    plus_strand_peaks = set()
    minus_strand_peaks = set()

    with open(zscore_peak_file, 'rU') as f:
        for line in f:
            fields = line.split('\t')
            if fields[1] != mode:
                continue
            if fields[0] == '+':
                plus_strand_peaks.add(int(fields[2][:-1]))
            elif fields[0] == '-':
                minus_strand_peaks.add(int(fields[2][:-1]))

    return plus_strand_peaks, minus_strand_peaks


def is_in_coding(position, strand, plus_genes, minus_genes):
    '''
    Returns True if a peak is within a coding region, based on reference dictionary
    containing key: gene_name, vals: (start, end)
    Note: the "start" in this interval is always the smaller number, not biological start (e.g. start site of transcription)
    This means So 5'-3' orientation not consistent in these intervals betwen strand!

    :param position: position of peak in genome
    :param strand: strand of position
    :param plus_genes: reference dictionary containing key: gene_name, vals: (start, end)
    :param minus_genes: reference dictionary containing key: gene_name, vals: (start, end)
    '''
    if strand == '+':
        for gene in plus_genes.values():
            if int(position) > int(gene[0]) and int(position) < int(gene[1]):
                return True
        return False
    elif strand == '-':
        for gene in minus_genes.values():
            if int(position) > int(gene[0]) and int(position) < int(gene[1]):
                return True
        return False



def is_near_coding(position, strand, threshold, plus_genes, minus_genes):

    if strand == '+':
        for gene in plus_genes.values():
            if int(position) > int(gene[0])-threshold and int(position) < int(gene[1])+threshold:
                return True
        return False
    elif strand == '-':
        for gene in minus_genes.values():
            if int(position) > int(gene[0])-threshold and int(position) < int(gene[1])+threshold:
                return True
        return False

