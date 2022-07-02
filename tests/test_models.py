import pytest
from epiread_tools.naming_conventions import *
import numpy as np
from deconvolution_models.celfie import em as celfie
from deconvolution_models.celfie_plus import READMeth as celfie_plus
from deconvolution_models.epistate import READMeth as epistate
from deconvolution_models.epistate_plus import READMeth as epistate_plus
import math
def test_celfie_atlas():
    pass

def test_celfie_plus_atlas():
    pass

def test_epistate_atlas():
    pass

def test_simple_deconvolution():
    beta_tm = [np.zeros((2,2)), np.zeros((2,3))]
    [x.fill(UNMETHYLATED) for x in beta_tm]
    beta_tm[0][:,0], beta_tm[1][:,1] = METHYLATED, METHYLATED
    mixture = [np.zeros((100,2)),np.zeros((100,3))]
    [x.fill(UNMETHYLATED) for x in mixture]
    mixture[0][:20,:], mixture[1][:20,:] = METHYLATED, METHYLATED
    r = celfie_plus(mixture, beta_tm, num_iterations=1000,
                    convergence_criteria=0.001)
    alpha, _ = r.two_step()
    assert (alpha-np.array([0.2, 0.8])).sum() < 0.03
    x = np.array([np.sum(m==METHYLATED) for m in mixture]) #TODO:fix
    x_depths = np.array([m.shape[0] for m in mixture])
    y = np.array([100,0,0,100]).reshape((2,2))
    y_depths = [np.array([100,100])]
    alpha, _ = celfie(
            x, x_depths, y, y_depths, 1000, 0.001
        )
    assert (alpha-np.array([0.2, 0.8])).sum() < 0.03

