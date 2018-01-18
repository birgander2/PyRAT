import numpy as np
from itertools import product
from functools import reduce


class Blocxy:
    """
    Class for performing automatic blockwise processing of an n-dimensional numpy array.
    using overlapping patches, which are fused by a weigthing function
    (default: triangular weigths). User supplied weighting functions are supported.

    Example:
    array = numpy.random.rand(1000, 1000)
    blx = Blocxy(array, (32, 32), (16, 16))
    for block in blx.getiterblocks():
        block = myfunc(block)
        blx.setiterblocks(block)
    result = blx.getresult()

    :arg array: The image to filter (numpy.ndarray)
    :param blocksize: the blocksize for processing
    :type blocksize: tuple of integers
    :param stepsize: the offset between neigbouring blocks
    :type stepsize: tuple of integers
    :param margin: optional - the number of border pixels to ignore (default=0)
    :type margin: tuple of integers
    :param wfunc: a user-supplied 1D weighting function (optional)
    :type wfunc: function

    :author: Andreas Reigber
    """
    def __init__(self, array, blocksize, stepsize=None, margin=None, wfunc=None):
        if wfunc is None:
            wfunc = self.triangle
        if stepsize is None:
            stepsize = blocksize
        if isinstance(blocksize, int):
            blocksize = (blocksize,)
        if margin is None:
            margin = (0, ) * len(blocksize)
        if isinstance(margin, int):
            margin = (margin, )

        if isinstance(stepsize, int):
            stepsize = (stepsize,)
        blocksize = np.array(blocksize)
        stepsize = np.array(stepsize)
        shp = array.shape

        shptest = {len(shp), len(blocksize), len(stepsize), len(margin)}
        if len(shptest) != 1:
            raise ValueError("shapes do not match!")
        if np.any(blocksize > shp):
            raise ValueError("blocksize larger than array!")

        slices = [[] for i in range(len(blocksize))]
        for k in range(len(blocksize)):
            nblocks = ((shp[k] - blocksize[k]) // stepsize[k]) + 1
            for l in range(nblocks):
                start = l * stepsize[k]
                ende = start + blocksize[k]
                slices[k].append(slice(start, ende))
            slices[k].append(slice(shp[k] - blocksize[k], shp[k]))

        self.shape = tuple([len(s) for s in slices])
        self.blocksize = blocksize
        self.ndim = len(slices)
        slices = product(*slices)
        self.blocks = []
        self.slices = []
        for sl in slices:
            self.blocks.append(array[sl].copy())
            self.slices.append(sl)
        self.nblocks = len(self.blocks)

        wv = []
        for k in range(len(blocksize)):
            func = wfunc(blocksize[k] - 2*margin[k])
            func = np.append(np.append(np.zeros(margin[k]), func), np.zeros(margin[k]))
            wv.append(func)
        self.blockweight = reduce(np.multiply.outer, wv)
        self.iterindex = 0

        self.acc_test = np.zeros(self.nblocks, dtype=bool)
        self.acc_data = np.zeros_like(array)
        self.acc_blkw = np.zeros_like(array, dtype='f4')

    def getblock(self, idx, index=False):
        """
        Returns the n-th block of the array.
        :param idx: the block number to return. Can be specified as flat index or as tuple.
        :param index: if set to true, additionally the centre coordinates of the block are returned.
        """
        if isinstance(idx, int):
            if idx > self.nblocks:
                raise ValueError("block index too large!")
        elif len(idx) == len(self.shape):
            if np.any(idx > self.shape):
                raise ValueError("block index too large!")
            idx = np.ravel_multi_index(idx, self.shape)
        else:
            raise ValueError("shapes do not match!")
        block = self.blocks[idx]
        if index is True:
            current = self.slices[idx]
            for k in range(self.ndim):
                centre = tuple([(sl.stop + sl.start) / 2 for sl in current])
            block = (block, centre)
        return block

    def getiterblocks(self, **kwargs):
        """
        Returns an iterator of all input blocks.
        :param index: if set to true, additionally the centre coordinates of the blocks are returned.
        """
        idx = 0
        while idx < self.nblocks:
            self.iterindex = idx
            yield self.getblock(idx, **kwargs)
            idx += 1

    def setiterblocks(self, block):
        """
        Sets (accumulates) the blocks provided by getiterblocks() and processed after.
        """
        self.setblock(block, self.iterindex)

    def setblock(self, block, idx):
        """
        Sets (accumulates) the n-th block of the result.
        :param block: the processed block
        :param idx: the block number to return. Can be specified as flat index or as tuple.
        """
        if isinstance(idx, int):
            if idx > self.nblocks:
                raise ValueError("block index too large!")
        elif len(idx) == len(self.shape):
            if np.any(idx > self.shape):
                raise ValueError("block index too large!")
            idx = np.ravel_multi_index(idx, self.shape)
        else:
            raise ValueError("shapes do not match!")
        if self.acc_test[idx] is True:
            raise ValueError("this block has already been set!")
        self.acc_test[idx] = True
        self.acc_data[self.slices[idx]] += block * self.blockweight
        self.acc_blkw[self.slices[idx]] += self.blockweight

    def getresult(self):
        """
        Returns the (accumulated) output array.
        """
        with np.errstate(divide='ignore'):
            result = self.acc_data / self.acc_blkw
        result[np.isnan(result)] = 0.0
        return result

    @staticmethod
    def triangle(length, min=0.01, max=1.0):
        section = length // 2
        if length % 2 != 0:
            val = np.linspace(min, max, section + 1)
            return np.append(val, val[-2::-1])
        else:
            val = np.linspace(min, max, section)
            return np.append(val, val[-1::-1])

