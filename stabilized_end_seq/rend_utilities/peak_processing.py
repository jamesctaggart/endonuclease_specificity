from plastid import GenomeArray, GenomicSegment
import numpy.ma as ma


def peak_height(peak_position, strand, genome_array, mode='position', chrom = 'NC_000964.3'):
    if mode == 'position':
        peak_gs = GenomicSegment(chrom, peak_position-1, peak_position, strand)
        return genome_array[peak_gs]
    elif mode == 'max':
        peak_gs = GenomicSegment(chrom, peak_position-3, peak_position+3, strand)
        return max(genome_array[peak_gs])


def integrated_peak_height(peak_position, strand, threshold, genome_array, mask=None, return_mask=False, chrom ='NC_000964.3'):
    """
    Given a particular peak position, calculate the sum of all positions with a value
    of at least threshold * peak_height (where threshold in [0,1]).
    :param peak_position: Postion of peak determined by peak calling (used as center of window)
    :param genome_array: GenomeArray object (data imported from wiggle file)
    :return: Sum of peak
    """

    local_region = GenomicSegment(chrom, peak_position-3, peak_position+3, strand)
    peak = genome_array[local_region]
    if mask is None:
        mask = [peak <= (threshold*max(peak))]

    peak = ma.array(peak, mask=mask)

    if return_mask:
        return ma.sum(peak), mask
    else:
        return ma.sum(peak)

def peak_height_ratio(peak_position, strand, normalization, ds1, ds2, normalize = True, mode='integrated'):

    if mode == 'single':
        ratio = peak_height(peak_position, strand, ds1, mode='max') / peak_height(peak_position, strand, ds2, mode='max')
    elif mode == 'integrated':
        numerator, mask = integrated_peak_height(peak_position, strand, 0.5, ds1, return_mask=True)
        ratio = numerator / integrated_peak_height(peak_position, strand, 0.5, ds2, mask=mask)
    else:
        raise ValueError('Please specify peak quantification mode as "single" or "integrated"')

    if normalize == True:
        return ratio / normalization
    else:
        return ratio