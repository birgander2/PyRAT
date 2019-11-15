import pyrat
import numpy as np
import scipy as sp
import ast
import matplotlib.pyplot as plt

class EntalpClass(pyrat.FilterWorker):
    """
    Entropy / Alpha clustering into 9 physically motivated classes. For details see:
    "An Entropy Based Classification Scheme for Land Applications of Polarimetric SAR", Shane Robert Cloude et al,
    IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 35, NO. 1, JANUARY 1997

    This module expects as input the polarimetric entropy, alpha angle and anisotropy (PyRAT module entalpani)

    :author: Joel Amao
    """
    gui = {'menu': 'PolSAR|Classification', 'entry': 'Entropy/Alpha classification'}
    para = [
        {'var': 'alpha', 'value': True, 'type': 'bool', 'text': 'Use alpha-max (alpha-mean otherwise)'},
        {'var': 'ord', 'value': False, 'type':'bool', 'text': 'Re-arrange classes based on scattering type'}
    ]

    def __init__(self, *args, **kwargs):
        super(EntalpClass, self).__init__(*args, **kwargs)
        self.name = "H/a classification"
        self.blockprocess = True
        self.scaling_hint = 'labels'

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        tmp = np.zeros_like(array[0, ...], dtype='uint8')
        ent = array[0, ...]
        if self.alpha:
            alp = np.degrees(array[1, ...])
        else:
            alp = np.degrees(array[2, ...])

        if self.ord:
            # Wasser
            tmp[(ent < 0.5) & (alp < 42.5)] = 0                             # Low Entropy Surface Scattering
            tmp[(ent >= 0.5) & (ent < 0.9) & (alp < 40)] = 1                # Medium Entropy Surface Scattering
            tmp[(ent >= 0.9) & (alp < 40)] = 2                              # High Entropy Surface Scattering (Unfeas.)

            # Veg
            tmp[(ent < 0.5) & (alp >= 42.5) & (alp < 47.5)] = 3             # Low Entropy Dipole Scattering
            tmp[(ent >= 0.5) & (ent < 0.9) & (alp >= 40) & (alp < 50)] = 4  # Medium Entropy Vegetation Scattering
            tmp[(ent >= 0.9) & (alp >= 40) & (alp < 55)] = 5                # High Entropy Vegetation Scattering
            tmp[(ent >= 0.9) & (alp >= 55)] = 6                             # High Entropy Multiple Scattering

            # Geb
            tmp[(ent >= 0.5) & (ent < 0.9) & (alp >= 50)] = 7               # Medium Entropy Multiple Scattering
            tmp[(ent < 0.5) & (alp >= 47.5)] = 8                            # Low Entropy Multiple Scattering Events
        else:
            tmp[(ent >= 0.9) & (alp >= 55)] = 0                             # High Entropy Multiple Scattering
            tmp[(ent >= 0.9) & (alp >= 40) & (alp < 55)] = 1                # High Entropy Vegetation Scattering
            tmp[(ent >= 0.9) & (alp < 40)] = 2                              # High Entropy Surface Scattering (unfeas.)

            tmp[(ent >= 0.5) & (ent < 0.9) & (alp >= 50)] = 3               # Medium Entropy Multiple Scattering
            tmp[(ent >= 0.5) & (ent < 0.9) & (alp >= 0) & (alp < 50)] = 4   # Medium Entropy Vegetation Scattering
            tmp[(ent >= 0.5) & (ent < 0.9) & (alp < 40)] = 5                # Medium Entropy Surface Scattering

            tmp[(ent < 0.5) & (alp >= 55)] = 6                              # Low Entropy Multiple Scattering Events
            tmp[(ent < 0.5) & (alp >= 42.5) & (alp < 47.5)] = 7             # Low Entropy Dipole Scattering
            tmp[(ent < 0.5) & (alp < 42.5)] = 8                             # Low Entropy Surface Scattering

        meta['CH_name'] = ['Labels']
        if self.ord:
            meta['labels'] = {0: 'Low Entropy Surface Scattering', 1: 'Medium Entropy Surface Scattering',
                              2: 'High Entropy Surface Scattering', 3: 'Low Entropy Dipole Scattering',
                              4: 'Medium Entropy Vegetation Scattering', 5: 'High Entropy Vegetation Scattering',
                              6: 'High Entropy Multiple Scattering', 7: 'Medium Entropy Multiple Scattering',
                              8: 'Low Entropy Multiple Scattering Events'}
        else:
            meta['labels'] = {0: 'High Entropy Multiple Scattering', 1: 'High Entropy Vegetation Scattering',
                              2: 'High Entropy Surface Scattering', 3: 'Medium Entropy Multiple Scattering',
                              4: 'Medium Entropy Vegetation Scattering', 5: 'Medium Entropy Surface Scattering',
                              6: 'Low Entropy Multiple Scattering Events', 7: 'Low Entropy Dipole Scattering',
                              8: 'Low Entropy Surface Scattering'}
        meta['labels'] = str(meta['labels'])
        return tmp


@pyrat.docstringfrom(EntalpClass)
def entalpclass(*args, **kwargs):
    return EntalpClass(*args, **kwargs).run(*args, **kwargs)


class Wishart(pyrat.Worker):
    """
    Wishart k-means clustering into n classes using a random/given initialisation. If given 2 input layers, it will take
    the first one as the data layer and the second one as the initialisation layer, otherwise it will perform a random
    initialisation.

    :author: Andreas Reigber
    """
    gui = {'menu': 'PolSAR|Classification', 'entry': 'Unsupervised Wishart (K-means)'}
    para = [
        {'var': 'nclass', 'value': 8, 'type': 'int', 'range': [2, 99], 'text': '# of classes'},
        {'var': 'niter', 'value': 10, 'type': 'int', 'range': [2, 99], 'text': '# of iterations'}
    ]

    def run(self, *args, **kwargs):

        if isinstance(self.layer, list): # Means the init layer is also active
            l_cov = self.layer[0]
            l_init = self.layer[1]
        else:
            l_cov = self.layer
            l_init = []
            outsize = pyrat.data.shape[-2:]

        meta = pyrat.getmeta(layer=l_cov)

        # STEP0: Do random initialisation if l_init is empty
        if len(l_init) == 0:
            l_init = self.layer_fromfunc(self.init_random, size=outsize, nclass=self.nclass)

        P = pyrat.tools.ProgressBar('  ' + self.name, self.niter)
        P.update(0)
        for iter in range(self.niter):
            # STEP1: Calculate cluster centres (and their frequency)
            pyrat.activate([l_cov, l_init], silent=True)
            cc, nc = self.layer_accumulate(self.calc_centers, nclass=self.nclass, combine=self.combine_mean)
            pyrat.delete(l_init, silent=True)

            # STEP2: Eliminate empty classes
            cc = [cc[k] for k, n in enumerate(nc) if n != 0]
            nc = [nc[k] for k, n in enumerate(nc) if n != 0]

            # STEP3: Calculate class memberships
            pyrat.activate(l_cov, silent=True)
            l_init = self.layer_process(self.assign_classes, centers=cc)
            P.update(iter + 1)
        del P
        meta['Class Centers'] = cc
        pyrat.setmeta(meta, layer=l_init)
        pyrat.activate(l_init)
        return l_init

    @staticmethod
    def assign_classes(data, centers=None, **kwargs):
        dist = np.empty((len(centers),) + data.shape[-2:])
        for k, center in enumerate(centers):
            dist[k, ...] = np.abs(np.trace(np.einsum('ij...,jk->ik...', data, np.linalg.inv(center)))) \
                           + np.log(np.abs(np.linalg.det(center)))
        out = np.ndarray.astype(np.argmin(dist, axis=0), dtype=np.ubyte)
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


class ClassExtract(pyrat.FilterWorker):
    """
    -
    Input: Classification mask
    Output: Classification mask with just the given classes

    :author: Joel Amao
    """
    gui = {'menu': 'PolSAR|Classification', 'entry': 'Extract Class'}
    para = [
        {'var': 'classes', 'value': '0', 'type': 'str', 'text': 'Class to extract \n'
                                                                '(use commas for multiple classes, e.g. 3,4,...)'}
    ]

    def __init__(self, *args, **kwargs):
        super(ClassExtract, self).__init__(*args, **kwargs)
        self.name = "Extract Class"
        self.blockprocess = True
        self.scaling_hint = 'labels'

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        tmp = np.zeros_like(array, dtype='int16')
        tmp -= 9999
        classes = self.classes.split(',')
        classes = list(map(int, classes))

        for k in range(len(classes)):
            tmp[array == (classes[k])] = classes[k]

        meta['CH_name'] = ['Labels']
        # Remove labels from the metadata if possible
        if 'labels' in meta:
            new_dic = ast.literal_eval(meta['labels'])
            l = [{k:v for k, v in i.items() if k in set(classes)} for i in [new_dic]]
            meta['labels'] = str(l[0])
        return tmp


@pyrat.docstringfrom(ClassExtract)
def classextract(*args, **kwargs):
    return ClassExtract(*args, **kwargs).run(*args, **kwargs)

class ClassAssign(pyrat.FilterWorker):

    def __init__(self, *args, **kwargs):
        super(ClassAssign, self).__init__(*args, **kwargs)
        self.name = "Assign Classes"
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        centers = self.centers
        dist = np.empty((len(centers),) + array.shape[-2:])
        for k, center in enumerate(centers):
            dist[k, ...] = np.abs(np.trace(np.einsum('ij...,jk->ik...', array, np.linalg.inv(center)))) \
                           + np.log(np.abs(np.linalg.det(center)))
        out = np.ndarray.astype(np.argmin(dist, axis=0), dtype=np.ubyte)
        return out

@pyrat.docstringfrom(ClassAssign)
def classassign(*args, **kwargs):
    return ClassAssign(*args, **kwargs).run(*args, **kwargs)

