import pyrat
import numpy as np
import numpy.linalg as la


def block_det(cov):
    """
    :author: Marc Jaeger
    """
    if cov.ndim == 2:
        if np.iscomplexobj(cov):
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


def C_to_T(array, pol):
    """
    Transforms a covariance matrices (3x3 or 4x4) into coherency matrices
    :author: Andreas Reigber

    :param array: The covariance matrices (3,3,n,m) or (4,4,n,m)
    :param pol: list with polarisation descriptor strings
    :return: coherency matrices (3,3,n,m) or (4,4,n,m)
    """
    oarray = np.empty_like(array)
    if array.shape[0] == 3:
        idx_hh = pol.index('HHHH*') % 3
        idx_vv = pol.index('VVVV*') % 3
        idx_xx = pol.index('XXXX*') % 3
        idx = np.array([idx_hh, idx_vv, idx_xx])
        D = np.array([[1.0, 1.0, 0.0],
                      [1.0, -1.0, 0.0],
                      [0.0, 0.0, np.sqrt(2)]], dtype='f4') / np.sqrt(2.0)
        oarray = array[idx, ...][:, idx, ...]
        oarray = np.rollaxis(np.rollaxis(oarray, 0, start=4), 0, start=4)
        oarray = np.matmul(D, np.matmul(oarray, np.conj(D.T)))
        oarray = np.rollaxis(np.rollaxis(oarray, 3), 3)
        newpol = ['P1P1*', 'P1P2*', 'P1P3*', 'P2P1*', 'P2P2*', 'P2P3*', 'P3P1*', 'P3P2*', 'P3P3*']
        for k, p in enumerate(newpol):
            pol[k] = p
    elif array.shape[0] == 4:
        idx_hh = pol.index('HHHH*') % 4
        idx_vv = pol.index('VVVV*') % 4
        idx_hv = pol.index('HVHV*') % 4
        idx_vh = pol.index('VHVH*') % 4
        idx = np.array([idx_hh, idx_vv, idx_hv, idx_vh])
        D = np.array([[1.0, 1.0, 0.0, 0.0],
                      [1.0, -1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0, 1.0],
                      [0.0, 0.0, -1j, 1j]], dtype='c8') / np.sqrt(2.0)
        oarray = array[idx, ...][:, idx, ...]
        oarray = np.rollaxis(np.rollaxis(oarray, 0, start=4), 0, start=4)
        oarray = np.matmul(D, np.matmul(oarray, np.conj(D.T)))
        oarray = np.rollaxis(np.rollaxis(oarray, 3), 3)
        newpol = ['P1P1*', 'P1P2*', 'P1P3*', 'P1P4*', 'P2P1*', 'P2P2*', 'P2P3*', 'P2P4*',
                  'P3P1*', 'P3P2*', 'P3P3*', 'P3P4*', 'P4P1*', 'P4P2*', 'P4P3*', 'P4P4*']
        for k, p in enumerate(newpol):
            pol[k] = p
    return oarray
