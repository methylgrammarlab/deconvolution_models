'''
MIT License

Copyright (c) 2022 irene unterman and ben berman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import pytest
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *
import numpy as np
from deconvolution_models.celfie import em as celfie
from deconvolution_models.celfie_ish import CelfieISH as celfie_plus
# from deconvolution_models.epistate import READMeth as epistate
# from deconvolution_models.epistate_plus import Epistate as epistate_plus

@pytest.fixture
def plus_deterministic_input():
    beta_tm = [np.zeros((2,2)), np.ones((2,3))]
    beta_tm[0][0,:], beta_tm[1][0,:] = 1,0
    mixture = [np.zeros((100,2)),np.zeros((100,3))]
    [x.fill(UNMETHYLATED) for x in mixture]
    mixture[0][:20,:], mixture[1][:80,:] = METHYLATED, METHYLATED
    return beta_tm, mixture

@pytest.fixture
def deterministic_input():
    y_depths = np.array([100,100,100,100]).reshape((2,2))
    y = np.identity(2)*100
    x_depths = np.array([100,100])
    x = np.array([20,80])
    return x, x_depths, y, y_depths

def test_celfie_atlas():
    pass

def test_celfie_plus_atlas():
    pass

def test_epistate_atlas():
    pass

def test_simple_celfie_deconvolution(deterministic_input):
    x, x_depths, y, y_depths = deterministic_input
    alpha, *_ = celfie(
        x.reshape((1, -1)), x_depths.reshape((1, -1)),y.T, y_depths.T,
        1000, 0.001
    )
    assert np.sum(alpha-np.array([0.2, 0.8]))<0.1

def test_simple_celfie_plus_deconvolution(plus_deterministic_input):
    beta_tm, mixture = plus_deterministic_input
    r = celfie_plus(mixture, beta_tm, num_iterations=1000,
                    convergence_criteria=0.001)
    alpha, _ = r.two_step()
    assert np.sum(alpha-np.array([0.2, 0.8]))<0.1