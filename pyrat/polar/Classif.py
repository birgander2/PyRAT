from __future__ import print_function
import pyrat
import numpy as np
from pyrat.tools import ProgressBar


class Wishart(pyrat.Worker):
    """
    Wishart k-means clustering into n classes using a random initialisation

    :author: Andreas Reigber
    """
    gui = {'menu': 'PolSAR|Classification', 'entry': 'Unsupervised Wishart'}
    para = [
        {'var': 'nclass', 'value': 8, 'type': 'int', 'range': [2, 99], 'text': '# of classes'},
        {'var': 'niter', 'value': 10, 'type': 'int', 'range': [2, 99], 'text': '# of iterations'}
    ]

    def run(self, *args, **kwargs):
        l_cov = pyrat.data.active
        outsize = pyrat.data.shape[-2:]

        # STEP0: Random initialisation
        l_init = self.layer_fromfunc(self.init_random, size=outsize, nclass=self.nclass)
        P = ProgressBar('  ' + self.name, self.niter)
        P.update(0)
        for iter in range(self.niter):
            # STEP1: Calculate cluster centres (and their frequency)
            pyrat.activate([l_cov, l_init], silent=True)
            cc, nc = self.layer_accumulate(self.calc_centers, nclass=self.nclass, combine=self.combine_mean)
            pyrat.delete(l_init, silent=True)

            # STEP2: Eliminate empty classes
            for k, n in enumerate(nc):
                if n == 0:
                    del cc[k]
                    del nc[k]

            # STEP3: Calculate class memberships
            pyrat.activate(l_cov, silent=True)
            l_init = self.layer_process(self.assign_classes, centers=cc)
            P.update(iter + 1)
        del P
        pyrat.activate(l_init)
        return l_init

    @staticmethod
    def assign_classes(data, centers=None, **kwargs):
        dist = np.empty((len(centers),) + data.shape[-2:])
        for k, center in enumerate(centers):
            dist[k, ...] = np.abs(np.trace(np.einsum('ij...,jk->ik...', data, np.linalg.inv(center)))) \
                           + np.log(np.abs(np.linalg.det(center)))
        out = np.argmin(dist, axis=0)
        return out

    @staticmethod
    def calc_centers(data, nclass=8, **kwargs):
        cov = data[0]
        cmshp = data[1]
        cc = []  # class centres
        nc = []  # number of members

        for k in range(nclass + 1):
            if k == 0:
                sel = cov[..., (cmshp < 1) | (cmshp > nclass + 1)]
            else:
                sel = cov[..., (cmshp == k)]
            nc.append(sel.shape[-1])  # number of members
            if sel.shape[-1] != 0:
                cc.append(np.mean(sel, axis=2))  # mean class centres
            else:
                cc.append(np.zeros(sel.shape[0:2], dtype=sel.dtype))
        return cc, nc

    @staticmethod
    def init_random(nclass=8, size=(1, 1), **kwargs):
        np.random.seed()
        return np.random.randint(1, nclass + 1, size=size)

    @staticmethod
    def combine_mean(values):
        npix = [0] * len(values[0][1])
        for val in values:
            npix = [x + y for x, y in zip(npix, val[1])]
        cc = [np.zeros_like(val[0][0])] * len(values[0][0])
        for val in values:
            cc = [c + x * y / n if n != 0 else c for c, n, x, y in zip(cc, npix, val[0], val[1])]
        return cc, npix


@pyrat.docstringfrom(Wishart)
def wishart(*args, **kwargs):
    return Wishart(*args, **kwargs).run(*args, **kwargs)
