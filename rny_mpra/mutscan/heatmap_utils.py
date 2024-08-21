import numpy as np
from collections import defaultdict

def nest_dm_dict(dm_dict, log2_average=True):
    '''
    Split dictionary of double mutants m1_m2 into nested dictionary
    that can be called with dict[m1][m2]
    :param dm_dict: Dictionary of form "m1_m2":[list - value for each barcoded variant]
    :return: Nested double mutant dictionary
    '''

    nested_dict = defaultdict(dict)
    for dm in dm_dict:
        m1 = dm.split('_')[0]
        m2 = dm.split('_')[1]
        if log2_average:
            nested_dict[m1][m2] = np.log2(np.median(dm_dict[dm]))
        else:
            nested_dict[m1][m2] = dm_dict[dm]

    return nested_dict