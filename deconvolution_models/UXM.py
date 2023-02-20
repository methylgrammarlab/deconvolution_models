'''
portions of this code are copied from
Loyfer, N., Magenheim, J., Peretz, A. et al.
A DNA methylation atlas of normal human cell types. Nature 613, 355–364 (2023). https://doi.org/10.1038/s41586-022-05580-6

https://github.com/nloyfer/UXM_deconv

Portions of the code are licensed under:

                    GNU AFFERO GENERAL PUBLIC LICENSE
                       Version 3, 19 November 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>

'''
import numpy as np
import pandas as pd
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *
from scipy import optimize

def uxm(atlas, samp, weights=None):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :param atlas: the atlas nested list
    :param weights: optional weight per region
    :return: the mixture coefficients
    """
    #TODO: remove missing sites from both sample and atlas

    if weights is not None: #TODO:fix
        ready_samp = [s*weights[i] for (i, s) in enumerate(samp)]
        ready_atlas = [s*weights[i] for (i, s) in enumerate(atlas)]
    else:
        ready_samp = samp
        ready_atlas = atlas
    mixture, residual = optimize.nnls(ready_atlas, ready_samp) #TODO add maxiter?
    mixture /= np.sum(mixture)
    return mixture

