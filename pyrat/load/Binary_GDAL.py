import pyrat
import logging
from osgeo import gdal
import numpy as np


def gdal_import_formats():
    formats = ""
    for i in range(gdal.GetDriverCount()):
        drv = gdal.GetDriver(i)
        if drv.GetMetadataItem(gdal.DCAP_RASTER) and drv.GetMetadataItem(gdal.DCAP_OPEN):
            typ = drv.GetMetadataItem(gdal.DMD_LONGNAME).split("(")[0].strip()
            exts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
            extensions = ""
            if exts is None:
                extensions = "*"
            else:
                for ex in exts.split(" "):
                    if ex == "":
                        ex = "*"
                    extensions += "*." + ex + " "
            formats += typ + " (" + extensions[:-1] + ");;"
    if formats[-1] == ";":
        formats = formats[:-2]
    formats = ";;".join(sorted(formats.split(";;")))
    return formats


class GDALRASTER(pyrat.ImportWorker):
    """
    Import of data in (any?) GDAL raster format.
    Comment: Import of metadata somehow unclear. Probably needs more work!

    **author:** Andreas Reigber\n
    **status:** --alpha-- Mostly untested!
    """
    formats = gdal_import_formats()
    gui = {'menu': 'File|Import raster', 'entry': 'GDAL raster source'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'GDAL data file', 'extensions': formats},
            {'var': 'stack', 'value': True, 'type': 'bool', 'text': 'Stack channels'}]

    def __init__(self, *args, **kwargs):
        super(GDALRASTER, self).__init__(*args, **kwargs)
        self.name = "GDAL IMPORT"
        if len(args) == 1:
            self.file = args[0]

    def getsize(self, *args, **kwargs):
        self.ds = gdal.Open(self.file)
        if self.ds is not None:
            logging.info("GDAL IMPORT: file type identified as " + self.ds.GetDriver().LongName)
            self.band = []
            for band in range(self.ds.RasterCount):
                self.band.append(self.ds.GetRasterBand(band + 1))
            if self.stack is True:
                return self.ds.RasterCount, self.ds.RasterYSize, self.ds.RasterXSize
            else:
                return self.ds.RasterYSize, self.ds.RasterXSize
        else:
            logging.error("ERROR: GDAL format not recognised!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []
        for band in self.band:
            array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][0], win_ysize=self.blocksize))
        if len(array) == 1:
            return array[0]
        else:
            if self.stack is True:
                array = np.stack(array, axis=0)
            return array

    def close(self, *args, **kwargs):
        self.ds = None                                                         # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        metain = self.ds.GetMetadata()
        meta.update(metain)
        for band in self.band:
            metain = band.GetMetadata()
            meta.update(metain)
        return meta


@pyrat.docstringfrom(GDALRASTER)
def gdalraster(*args, **kwargs):
    return GDALRASTER(*args, **kwargs).run(*args, **kwargs)


class GeoTIFF(GDALRASTER):
    """
    GeoTiff format reader (using GDAL)
    """
    gui = {'menu': 'File|Import raster', 'entry': 'GeoTIFF'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'GeoTIFF file', 'extensions': 'GeoTIFF (*.tif *.tiff)'},
            {'var': 'stack', 'value': False, 'type': 'bool', 'text': 'Stack channels'}]

    def __init__(self, *args, **kwargs):
        super(GeoTIFF, self).__init__(*args, **kwargs)
        self.name = "GEOTIFF IMPORT"


@pyrat.docstringfrom(GeoTIFF)
def geotiff(*args, **kwargs):
    return GeoTIFF(*args, **kwargs).run(*args, **kwargs)


class ENVI(GDALRASTER):
    """
    GeoTiff format reader (using GDAL)
    """
    gui = {'menu': 'File|Import raster', 'entry': 'ENVI'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'ENVI file', 'extensions': 'ENVI .hdr Labelled (*.*)'},
            {'var': 'stack', 'value': False, 'type': 'bool', 'text': 'Stack channels'}]

    def __init__(self, *args, **kwargs):
        super(ENVI, self).__init__(*args, **kwargs)
        self.name = "ENVI IMPORT"


@pyrat.docstringfrom(ENVI)
def envi(*args, **kwargs):
    return ENVI(*args, **kwargs).run(*args, **kwargs)

