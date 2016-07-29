import pyrat
import numpy as np
from scipy.stats import chi2
from scipy.ndimage.filters import median_filter
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

        lnq = p * ((n1 + n2) * np.log(n1 + n2) - n1 * np.log(n1) - n2 * np.log(n2)) + \
              n1 * np.log(self.block_det(c1)) + \
              n2 * np.log(self.block_det(c2)) - \
              (n1 + n2) * np.log(self.block_det(c1 + c2))

        rho = 1 - (2 * p ** 2 - 1) * (1 / n1 + 1 / n2 - 1 / (n1 + n2)) / (6 * p)

        o_2 = p ** 2 * (p ** 2 - 1) * \
              (1 / n1 ** 2 + 1 / n2 ** 2 - 1 / (n1 + n2) ** 2) / (24 * rho ** 2) - \
              0.25 * p ** 2 * (1 - 1 / rho) ** 2

        lnq *= -2 * rho

        pfa = (1 - o_2) * chi2.cdf(lnq, p ** 2) + o_2 * chi2.cdf(lnq, p ** 2 + 4)

        return (pfa, median_filter(pfa > self.p_thresh, 3, mode='constant', cval=0))


@pyrat.docstringfrom(ChangeDet)
def changedet(*args, **kwargs):
    return ChangeDet(*args, **kwargs).run(*args, **kwargs)


class CoherChangeDet(PolsarWorker):
    """
    Simple coherent change detection from a PolInSAR covariance matrix.

    Does the following:
    - Masks invalid pixels (minimum covariance diagonal element smaller than noise level)
    - Thresholds maximum coherence (off-diagonal) for valid pixels
    - Low coherence indicates change...
    
    Returns one new layer with the change detection mask,
    
    :author: Marc Jaeger
    :param cov: Covariance matrices in the combined data set (4D np.ndarray)
    :type cov: complex64
    :param noise: The noise level in m^2. Either scalar or 2d, defaults to 1e-3 (-30 dB)
    :type noise: float32
    :param coh_thresh: Coherence threshold (defaults to 0.5)
    :type coh_thresh: float32
    """

    def __init__(self, *args, **kwargs):
        super(CoherChangeDet, self).__init__(*args, **kwargs)
        self.name = "COHERENT POLINSAR CHANGE DETECTOR"
        self.blockprocess = True
        self.blockoverlap = 0

        if 'noise' not in self.__dict__:
            self.noise = 1e-3
        if 'coh_thresh' not in self.__dict__:
            self.coh_thresh = 0.5

    def filter(self, cov, *args, **kwargs):
        cdim = cov.shape
        npol = cdim[0] // 2
        c_int = lambda n: cov[n, n, ...].real

        max_coh = np.zeros(cdim[2:], dtype='f8')
        inten = np.zeros(cdim[2:], dtype='f8')
        for n in range(npol):
            coh = np.abs(cov[n, n + npol, ...]) / np.sqrt(c_int(n) * c_int(n + npol))
            up_mask = (coh > max_coh)
            max_coh = (1 - up_mask) * max_coh + up_mask * coh
            inten = (1 - up_mask) * inten + 0.5 * up_mask * (c_int(n) + c_int(n + npol))

        mask = np.asarray((inten > self.noise) * (max_coh < self.coh_thresh), dtype='u1')

        return median_filter(mask, 3, mode='constant', cval=0)


@pyrat.docstringfrom(CoherChangeDet)
def coherchangedet(*args, **kwargs):
    return CoherChangeDet(*args, **kwargs).run(*args, **kwargs)
