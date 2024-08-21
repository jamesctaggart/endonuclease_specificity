import math
import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import operator as op
import matplotlib.pyplot as plt
from mutscan.plotting import get_fig_axes

def p_of_n_errors(rate,n,L):
    probability = math.comb(L,n)*math.pow(rate,n)*math.pow(1-rate,L-n)
    return probability

def p_of_specific_sequence(rate,n,L):
    probability = math.pow(rate/3.0,n)*math.pow(1-rate,L-n)
    return probability

def number_of_instances(rate,n,L,library_size):
    return p_of_specific_sequence(rate,n,L)*library_size


def get_123nmut_count(msp):
    df = msp.data
    return len(df.loc[df['n_mutations']==0]), len(df.loc[df['n_mutations']==1]), \
        len(df.loc[df['n_mutations']==2]), len(df.loc[df['n_mutations']==3])

def get_nmut_prob_dist(msp,axes=None,**kwargs):
    fig, ax = get_fig_axes(axes)
    data = msp.data.n_mutations.to_list()
    ax.hist(data, density=True, histtype='step',bins=range(min(data), max(data) + 1, 1))
    # ax.step(df['n_mutations'],**kwargs)