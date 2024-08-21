import numpy as np

def pitted_mean(array, x, L, g):
    '''

    :param array: data in array
    :param x: position of interest
    :param L: Averaging length
    :param g: gap
    :return:
    '''
    g = g + 1

    N = (L)-(2.0*g)
    half_window_length = (L/2)-g
    return (sum(array[x-half_window_length-g+1:x-g+1])+sum(array[x+g:x+half_window_length+g])) / N


def pitted_stdev(array, x, L, g):
    g = g + 1
    half_window_length = (L/2)-g
    stdev_window = np.concatenate((array[x-half_window_length-g+1:x-g+1],array[x+g:x+half_window_length+g]),axis=None)
    return np.std(stdev_window)


def calc_z(array, x, L, g, min_depth=True, reference_array=None, depth_threshold=0.25):
    '''
    If calculating z score from something where mean is expected to be near zero (e.g. subtracted two normalized datasets)
    use reference_array, which would just be one of the two datasets, as a threshold on depth instead.
    '''
    pm = pitted_mean(array, x, L, g)
    if reference_array is not None:
        ref_pm = pitted_mean(reference_array, x, L, g)
    else:
        ref_pm = pm

    if min_depth is True:
        if ref_pm < depth_threshold:
            return 0  # Force it to return a z-score of 0 if the mean is below a threshold to avoid spurious peaks from counting noise

    return (array[x] - pm) / pitted_stdev(array, x, L, g)


def collapse_adjacent_peaks(peak_array, data_array, dist, choose_smaller = False):
    '''
    :param dist: distance below which peaks will be combined
    '''
    collapsed_array = [peak_array[0]]
    for ii in range(1, len(peak_array)):
        if (peak_array[ii]-peak_array[ii-1])<dist:
            if choose_smaller is False:
                if data_array[peak_array[ii]] > data_array[peak_array[ii-1]]:
                    collapsed_array[-1] = peak_array[ii]
            else:
                if data_array[peak_array[ii]] < data_array[peak_array[ii-1]]:
                    collapsed_array[-1] = peak_array[ii]
        else:
            collapsed_array.append(peak_array[ii])
    return np.array(collapsed_array)

