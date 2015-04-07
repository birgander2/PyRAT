import pyrat, os
import logging
from osgeo import gdal
import glob
import xml.etree.ElementTree as ET


class TSX(pyrat.ImportWorker):
    """
    Import of TSX/TDX satellite data.

    **author:** Andreas Reigber\n
    **status:** --beta-- Not all relevant metadata are extracted. Currently only tested with
    TSX single-pol stripmap data, status of all other modes unclear...
    """

    gui = {'menu': 'File|Import spaceborne', 'entry': 'TerraSAR-X'}
    para = [{'var': 'dir', 'value': '', 'type': 'opendir', 'text': 'Product directory'}]

    def __init__(self, *args, **kwargs):
        super(TSX, self).__init__(*args, **kwargs)
        self.name = "TSX / TDX IMPORT"

    def getsize(self, *args, **kwargs):
        files = glob.glob(self.dir+"/*SAR*xml")
        if len(files) > 0:
            file = files[0]
            self.ds = gdal.Open(file)
            self.band = []
            for band in range(self.ds.RasterCount):
                self.band.append(self.ds.GetRasterBand(band+1))
            return self.ds.RasterYSize, self.ds.RasterXSize
        else:
            logging.error("ERROR: product directory not recognised!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []
        for band in self.band:
            array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][0], win_ysize=self.blocksize))
        if len(array) == 1:
            return array[0]
        else:
            return array

    def close(self, *args, **kwargs):
        self.ds = None          # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta  = {}
        metain  = self.ds.GetMetadata()
        meta.update(metain)
        meta['CH_pol'] = []
        for band in self.band:
            metain = band.GetMetadata()
            meta['CH_pol'].append(metain['POLARIMETRIC_INTERP'])

        # Read additional meta data from product xml file
        file = glob.glob(self.dir+"/*SAR*xml")[0]
        tree = ET.parse(file)
        root = tree.getroot()

        meta['sensor'] = root.find('generalHeader/mission').text
        meta['PRF'] = float(root.find('productSpecific/complexImageInfo/commonPRF').text)
        meta['rsf'] = float(root.find('productSpecific/complexImageInfo/commonRSF').text)
        meta['rd'] = float(root.find('productInfo/sceneInfo/rangeTime/firstPixel').text)
        meta['c0'] = 299792458.0
        foo = root.find('*/*/lookDirection').text
        meta['antdir'] = -1 if foo == 'RIGHT' else +1

        return meta
