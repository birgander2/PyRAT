from __future__ import print_function

import pyrat
import numpy as np
import numpy.linalg as la


class PolsarWorker(pyrat.FilterWorker):
    """
    A base class for manipulating PolSAR data. It provides a number of
    static utility functions for manipulating PolSAR covariance matrices.

    :author: Marc Jaeger
    """

    @staticmethod
    def block_det(cov):

        if (cov.ndim == 2):
            if (np.iscomplexobj(cov)):
                return cov.real
            else:
                return cov

        d = cov.shape[0]
        if d == 3:
            det = cov[0, 0, :, :] * cov[1, 1, :, :] * cov[2, 2, :, :] + \
                  cov[0, 1, :, :] * cov[1, 2, :, :] * cov[2, 0, :, :] + \
                  cov[0, 2, :, :] * cov[1, 0, :, :] * cov[2, 1, :, :] - \
                  cov[2, 0, :, :] * cov[1, 1, :, :] * cov[0, 2, :, :] - \
                  cov[2, 1, :, :] * cov[1, 2, :, :] * cov[0, 0, :, :] - \
                  cov[2, 2, :, :] * cov[1, 0, :, :] * cov[0, 1, :, :]
        else:
            det = la.det(cov.T).T

        return det.real
