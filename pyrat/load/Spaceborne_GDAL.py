import pyrat, os
import logging
from osgeo import gdal
import glob
import numpy as np


class ENVISAT(pyrat.ImportWorker):
    """
    Import of ENVISAT satellite data.

    **author:** Andreas Reigber\n
    **status:** --beta-- No metadata are extracted. Mostly untested!
    """

    gui = {'menu': 'File|Import spaceborne', 'entry': 'ENVISAT'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'Product file (*.N1)'}]

    def __init__(self, *args, **kwargs):
        super(ENVISAT, self).__init__(*args, **kwargs)
        self.name = "ENVISAT IMPORT"
        if len(args) == 1:
            self.file = args[0]

    def getsize(self, *args, **kwargs):
        self.ds = gdal.Open(self.file)
        if self.ds is not None:
            self.band = []
            for band in range(self.ds.RasterCount):
                self.band.append(self.ds.GetRasterBand(band + 1))
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
        self.ds = None  # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        meta['sensor'] = "ENVISAT"
        metain = self.ds.GetMetadata()
        meta.update(metain)
        for band in self.band:
            metain = band.GetMetadata()
            meta.update(metain)
        return meta


@pyrat.docstringfrom(ENVISAT)
def envisat(*args, **kwargs):
    return ENVISAT(*args, **kwargs).run(*args, **kwargs)


class PALSAR(pyrat.ImportWorker):
    """
    Import of PALSAR satellite data. Only level 1.1. and 1.5 are supported.

    **author:** Andreas Reigber\n
    **status:** --beta-- No metadata are extracted. Mostly untested!
    """

    gui = {'menu': 'File|Import spaceborne', 'entry': 'PALSAR'}
    para = [{'var': 'dir', 'value': '', 'type': 'opendir', 'text': 'Product directory'}]

    def __init__(self, *args, **kwargs):
        super(PALSAR, self).__init__(*args, **kwargs)
        self.name = "PALSAR IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def getsize(self, *args, **kwargs):
        volfile = glob.glob(self.dir + "/VOL*")
        if len(volfile) > 0:
            self.ds = gdal.Open(volfile[0])
            if self.ds is not None:
                self.band = []
                for band in range(self.ds.RasterCount):
                    self.band.append(self.ds.GetRasterBand(band + 1))
                return len(self.band), self.ds.RasterYSize, self.ds.RasterXSize
            else:
                logging.error("ERROR: product directory not recognised!")
                return False, False
        else:
            logging.error("ERROR: volume file not found!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []
        for band in self.band:
            array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][0], win_ysize=self.blocksize))
        out = np.empty((len(array),) + array[0].shape, dtype=array[0].dtype)
        for k in range(len(array)):
            out[k, ...] = array[k]
        out[~np.isfinite(out)] = 0
        return out.squeeze()

    def close(self, *args, **kwargs):
        self.ds = None  # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        meta['sensor'] = "PALSAR"
        metain = self.ds.GetMetadata()
        meta.update(metain)
        for band in self.band:
            metain = band.GetMetadata()
            meta.update(metain)
        return meta


@pyrat.docstringfrom(PALSAR)
def palsar(*args, **kwargs):
    return PALSAR(*args, **kwargs).run(*args, **kwargs)


class Radarsat2(pyrat.ImportWorker):
    """
    Import of Radarsat-2 satellite data.

    **author:** Andreas Reigber\n
    **status:** --beta-- No metadata are extracted. Mostly untested!
    """

    gui = {'menu': 'File|Import spaceborne', 'entry': 'Radarsat-2'}
    para = [{'var': 'dir', 'value': '', 'type': 'opendir', 'text': 'Product directory'}]

    def __init__(self, *args, **kwargs):
        super(Radarsat2, self).__init__(*args, **kwargs)
        self.name = "RADARSAT-2 IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def getsize(self, *args, **kwargs):
        volfile = glob.glob(self.dir + "/product.xml")
        if len(volfile) > 0:
            self.ds = gdal.Open(volfile[0])
            if self.ds is not None:
                self.band = []
                for band in range(self.ds.RasterCount):
                    self.band.append(self.ds.GetRasterBand(band + 1))
                return len(self.band), self.ds.RasterYSize, self.ds.RasterXSize
            else:
                logging.error("ERROR: product directory not recognised!")
                return False, False
        else:
            logging.error("ERROR: product.xml file not found!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []
        for band in self.band:
            array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][0], win_ysize=self.blocksize))
        out = np.empty((len(array),) + array[0].shape, dtype=array[0].dtype)
        for k in range(len(array)):
            out[k, ...] = array[k]
        out[~np.isfinite(out)] = 0
        return out.squeeze()

    def close(self, *args, **kwargs):
        self.ds = None  # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        meta['sensor'] = "Radarsat-2"
        metain = self.ds.GetMetadata()
        meta.update(metain)
        meta['CH_pol'] = []
        for band in self.band:
            metain = band.GetMetadata()
            meta['CH_pol'].append(metain['POLARIMETRIC_INTERP'])
            meta.update(metain)
        return meta


@pyrat.docstringfrom(Radarsat2)
def radarsat2(*args, **kwargs):
    return Radarsat2(*args, **kwargs).run(*args, **kwargs)


class Sentinel1(pyrat.ImportWorker):
    """
    Very basic import of Sentinel-1 satellite data. The current driver uses GDAL and therefore does not
    perform debursting and combination of subswaths. This routine needs to be improved in future.

    **author:** Andreas Reigber\n
    **status:** --beta-- Mostly untested!
    """

    gui = {'menu': 'File|Import spaceborne', 'entry': 'Sentinel-1 (primitive)'}
    para = [
        {'var': 'dir', 'value': '', 'type': 'opendir', 'text': 'Product directory'},
        {'var': 'swath', 'value': 1, 'type': 'int', 'range': [1, 3], 'text': 'Swath to Load'}
    ]

    def __init__(self, *args, **kwargs):
        super(Sentinel1, self).__init__(*args, **kwargs)
        self.name = "SENTINEL-1 IMPORT"


    def reader(self, *args, **kwargs):
        volfile = glob.glob(self.dir + "/manifest.safe")
        if len(volfile) > 0:
            ds = gdal.Open(volfile[0])
            if ds is not None:
                self.band = []
                if ds.GetMetadata()['PRODUCT_TYPE'] == 'SLC':
                    # Open sub-swath dataset and every band in it
                    sds_list = ds.GetSubDatasets()
                    self.ds = gdal.Open(sds_list[((self.swath-1)*3)+2][0])
                else:
                    # Keep existing functionality for GRD data
                    self.ds = ds
                # Load bands as usual
                for band in range(self.ds.RasterCount):
                    self.band.append(self.ds.GetRasterBand(band + 1))
                # FIXME: These variables are unused
                nswath = len(self.band)
                YSize = [band.YSize for band in self.band]
                XSize = [band.XSize for band in self.band]
            else:
                logging.error("ERROR: product directory not recognised!")
                return False, False
        else:
            logging.error("ERROR: manifest.save file not found!")
            return False, False
        array = []
        for band in self.band:
            array.append(band.ReadAsArray())

        meta = {}
        meta['sensor'] = "Sentinel-1"
        # Attach subdataset metadata with actual swath info
        metain = self.ds.GetMetadata()
        meta.update(metain)
        # TODO: Attach per sub-swath GCP info
        # TODO: Attach polarimetric interpretation to band

        return array, meta

    def close(self, *args, **kwargs):
        self.ds = None  # correct according to GDAL manual!!??


def sentinel1(*args, **kwargs):
    return Sentinel1(*args, **kwargs).run(*args, **kwargs)
