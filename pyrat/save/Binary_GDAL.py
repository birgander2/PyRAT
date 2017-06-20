import pyrat
import logging
from osgeo import gdal
import numpy as np


def gdal_export_formats():
    formats = ""
    format = ""
    driver = []
    dnames = []
    for i in range(gdal.GetDriverCount()):
        drv = gdal.GetDriver(i)
        if drv.GetMetadataItem(gdal.DCAP_RASTER) and drv.GetMetadataItem(gdal.DCAP_CREATE):
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
                format = typ + " (" + extensions[:-1] + ")"
                formats += format + ";;"
            driver.append(drv.GetDescription())
            dnames.append(format)
    if formats[-1] == ";":
        formats = formats[:-2]
    formats = ";;".join(sorted(formats.split(";;")))
    return formats, driver, dnames


def dict2gdalmeta(meta, nchannel=1):
    """
    Converts meta data (in dict) to gdal format (in strings)
    """
    meta_main = ""
    meta_channel = [""] * nchannel
    for key, val in meta.items():
        if key[0:3] == "CH_":
            for k, vsub in enumerate(val):
                meta_channel[k] += str(key[3:]) + "=" + str(vsub) + " "
        else:
            meta_main += str(key) + "=" + str(val) + " "
    for k, meta in enumerate(meta_channel):
        meta_channel[k] = meta.strip()

    return meta_main.strip(), meta_channel


NP2GDAL_CONVERSION = {
    "uint8": 1,
    "int8": 1,
    "uint16": 2,
    "int16": 3,
    "uint32": 4,
    "int32": 5,
    "float32": 6,
    "float64": 7,
    "complex64": 10,
    "complex128": 11,
}


class GDALRASTER(pyrat.ExportWorker):
    """
    Export to various GDAL raster formats

    When using in CLI mode, the format has to be selected together with the file name
    as a tuple. Example: file=("myfile.dat", "ENVI"). The format string corresponds to
    the GDAL driver name to use - actually the "short name" which you get with
    driver.GetDescription().

    When using in GUI mode, the format has to be selected in the fileselector.

    Comment: Export of metadata should work, but does it really work? Untested!

    **author:** Andreas Reigber\n
    **status:** --alpha--  Mostly untested!
    """
    # todo: Update routine to use block writing?

    formats, driver, dnames = gdal_export_formats()
    gui = {'menu': 'File|Export to raster', 'entry': 'GDAL raster formats'}
    para = [{'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save as GDAL', 'extensions': formats}]

    def __init__(self, *args, **kwargs):
        super(GDALRASTER, self).__init__(*args, **kwargs)
        self.name = "GDAL EXPORT"
        self.nthreads = 1
        if len(args) == 1:
            self.file = args[0]

    def writer(self, array, *args, **kwargs):
        if hasattr(pyrat, "app"):  # select correct GDAL driver
            drv_idx = self.dnames.index(self.file[1])  # GUI mode: long name
            driver = gdal.GetDriverByName(self.driver[drv_idx])
        else:
            driver = gdal.GetDriverByName(self.file[1])  # CLI mode: short name

        if driver is None:
            logging.error("Unknown GDAL driver: " + self.file[1])
            return

        meta = kwargs["meta"]

        shp = array.shape  # reshape to (nchannel, ny, nx)
        array = array.reshape((int(np.prod(shp[0:-2])),) + shp[-2:])
        shp = array.shape
        gdaltype = NP2GDAL_CONVERSION[array.dtype.name]

        ds = driver.Create(self.file[0], shp[2], shp[1], shp[0], gdaltype)
        for k, arr in enumerate(array):
            ds.GetRasterBand(k + 1).WriteArray(arr)

        gmeta, cmeta = dict2gdalmeta(meta, shp[0])
        ds.SetMetadata(gmeta, "")
        for k in range(ds.RasterCount):
            ds.GetRasterBand(k + 1).SetMetadata(cmeta[k], "")
        ds = None  # correct according to GDAL manual!!??


@pyrat.docstringfrom(GDALRASTER)
def gdalraster(*args, **kwargs):
    return GDALRASTER(*args, **kwargs).run(*args, **kwargs)
