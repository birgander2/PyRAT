from __future__ import print_function

import numpy as np
from scipy.stats import chi2
from . import PolsarWorker


class ChangeDet(PolsarWorker):
    """
    Incoherent change detection by comparing covariance matrices over 
    a temporal base line.
    
    Implements "A Test Statistic in the Complex Wishart Distribution
    and Its Application to Change Detection in Polarimetric SAR Data"
    by Conradsen et al.

    Returns two new layers:
    1. Change probabiity (between 0 and 1; values close to 1 indicate change)
    2. Change detection mask, obtained by thresholding above probability
    
    :author: Marc Jaeger
    :param c1: Covariance matrices in the first data set (4D np.ndarray)
    :type array: complex64
    :param n1: The effective number of looks in the first data set (per-pixel 2D np.ndarray or scalar)
    :type array: float32
    :param c2: Covariance matrices in the second data set (4D np.ndarray)
    :type array: complex64
    :param n2: The effective number of looks in the second data set (per-pixel 2D np.ndarray or scalar)
    :type array: float32
    :param p_thresh: The threshold for the binary change detection mask. Changes are indicated by values close to 1.0; default to 0.99.
    :type array: float32
    """

    def __init__(self, *args, **kwargs):
        super(ChangeDet, self).__init__(*args, **kwargs)
        self.name = "POLSAR CHANGE DETECTOR"
        self.blockprocess = True
        self.blockoverlap = 0
    
        if 'n1' not in self.__dict__:
            self.n1 = 16
        if 'n2' not in self.__dict__:
            self.n2 = 16
        if 'p_thresh' not in self.__dict__:
            self.p_thresh = 0.99


    def filter(self, data, *args, **kwargs):
    
        if (len(data) == 4):
            c1 = data[0]
            n1 = data[1]
            c2 = data[2]
            n2 = data[3]
        else:
            c1 = data[0]
            n1 = self.n1
            c2 = data[1]
            n2 = self.n2

        p = c1.shape[0]

        lnq = p*((n1+n2)*np.log(n1+n2)-n1*np.log(n1)-n2*np.log(n2)) + \
              n1*np.log(self.block_det(c1)) + \
              n2*np.log(self.block_det(c2)) - \
              (n1+n2)*np.log(self.block_det(c1+c2))

        rho = 1 - (2*p**2-1)*(1/n1 + 1/n2 - 1/(n1+n2))/(6*p)

        o_2 = p**2 * (p**2-1) * \
              (1/n1**2 + 1/n2**2 - 1/(n1+n2)**2)/(24*rho**2)  - \
              0.25 * p**2 * (1-1/rho)**2

        lnq *= -2*rho

        pfa = (1-o_2)*chi2.cdf(lnq,p**2) + o_2*chi2.cdf(lnq,p**2+4)

        return (pfa, pfa>self.p_thresh)

