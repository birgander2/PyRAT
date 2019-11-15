import pyrat
from pyrat.viewer.tools import sarscale, phascale, cohscale
import logging
import numpy as np
from PIL import Image
from pyrat.tools import colortables


class Pixmap(pyrat.ExportWorker):
    """
    Export to various pixmap formats.

    For possible parameters see source code ;-)
    """
    para = [
        {'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save to :'},
        {'var': 'chscl', 'value': True, 'type': 'bool', 'text': 'Scale channels indiviually'},
        {'var': 'method', 'value': 'amplitude', 'type': 'list',
         'range': ['amplitude', 'intensity', 'phase', 'coherence','minmax'], 'text': 'Method'},
        {'var': 'scaling', 'value': 2.5, 'type': 'float', 'range': [0.1, 20.0], 'text': 'SAR scaling factor'},
        {'var': 'palette', 'value': 'bw linear', 'type': 'list', 'range': colortables()[0], 'text': 'Color table'},
        {'var': 'order', 'type': 'list', 'value': '0', 'range': ['0', '2', '1'], 'text': 'Channel selection'}
    ]
    key = None

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        if self.key:
            self.name = "EXPORT TO " + self.key
        else:
            self.name = "EXPORT TO PIXMAP"
        if len(args) == 1:
            self.file = args[0]
        if self.order == '0':
            self.order = [0, 2, 1]
        elif self.order == '1':
            self.order = [0, 1, 2]
        else:
            self.order = [1, 2, 0]

    def writer(self, array, *args, **kwargs):
        if isinstance(self.file, tuple):  # remove file type if present
            self.file = self.file[0]

        if (self.method == 'amplitude' or self.method == 'intensity') and np.iscomplexobj(array):
            array = np.abs(array)

        array = np.squeeze(array)
        if array.ndim == 4:
            nchannels = np.prod(array.shape[0:array.ndim-2])
            array = np.reshape(array, (nchannels, )+array.shape[array.ndim-2:])

        if array.ndim == 3 and self.chscl == True:
            for k in range(array.shape[0]):
                array[k, ...] = array[k, ...] / np.mean(array[k, ...])

        if self.method == 'intensity' or self.method == 'amplitude':
            if self.method == 'amplitude':
                array **= 0.7
            if self.method == 'intensity':
                array **= 0.35
            out = sarscale(array, factor=self.scaling)
        elif self.method == 'phase':
            out = phascale(array)
        elif self.method == 'coherence':
            out = cohscale(array)
        elif self.method == 'minmax':
            start = array.min()
            end = array.max()
            out = np.uint8(np.clip((array - start) / (end - start) * 255, 0, 255))
        else:
            logging.error("Scaling method unknown")
            return False

        if array.ndim == 3:
            if out.shape[0] < 3:
                oshp = out.shape
                oshp[0] = 3
                out = np.resize(out, oshp)
            if out.shape[0] > 3:
                out = out[0:3, ...]
            out = np.rollaxis(np.rollaxis(out[self.order, ...], 2), 2)
        else:
            out = colortables(self.palette)[1][out]

        try:
            pilimg = Image.fromarray(out)
            pilimg.save(self.file, format=self.key)
            logging.info("FINISHED SAVING IMAGE")
            return True
        except IOError as err:
            logging.error("ERROR:" + str(err))
            return False
        else:
            logging.error("UNKNOWN ERROR")
            return False

    @classmethod
    def guirun(cls, viewer, title=None):
        if not hasattr(pyrat, "app"):
            return

        paras = [
            {'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save to :'},
            {'var': 'zoom', 'value': True, 'type': 'bool', 'text': 'Save only current window content'}
        ]

        if cls.key == 'JPEG':
            paras[0]['extensions'] = 'JPEG file (*.jpg *.jpeg)'
        elif cls.key == 'PNG':
            paras[0]['extensions'] = 'PNG file (*.png)'
        elif cls.key == 'TIFF':
            paras[0]['extensions'] = 'TIFF file (*.tif *.tiff)'
        elif cls.key == 'PDF':
            paras[0]['extensions'] = 'PDF file (*.pdf)'
        elif cls.key == 'EPS':
            paras[0]['extensions'] = 'EPS file (*.eps)'
        else:
            cls.key = None

        wid = pyrat.viewer.Dialogs.FlexInputDialog(paras, parent=viewer, title="Save to " + cls.key, doc=cls.__doc__)
        res = wid.exec_()

        if res == 1:
            para = {}
            for p in paras:
                para[p['var']] = p['value']

            if para['zoom'] is True:  # save only window content
                img = np.squeeze(viewer.img)
                try:
                    pilimg = Image.fromarray(img)
                    filename = para['file'][0] if isinstance(para['file'], tuple) else para['file']
                    pilimg.save(filename, format=cls.key)
                    return True
                except IOError as err:
                    logging.error("ERROR:" + str(err))
                    return False
            else:  # save entire image
                channels = [int(foo.split('D')[1]) for foo in viewer.config['rgblayer']]
                plugin = cls(file=para['file'], method=viewer.config['scaling'], scaling=viewer.sarscale,
                             palette=viewer.display[viewer.current]['palette'], order=channels)
                plugin.run()


@pyrat.docstringfrom(Pixmap)
def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(*args, **kwargs)


class JPG(Pixmap):
    """
    JPG format writer
    """
    gui = {'menu': 'File|Export to pixmap', 'entry': 'JPEG'}
    key = "JPEG"


@pyrat.docstringfrom(JPG)
def jpg(*args, **kwargs):
    JPG(*args, **kwargs).run(*args, **kwargs)


class PNG(Pixmap):
    """
    PNG format writer
    """
    gui = {'menu': 'File|Export to pixmap', 'entry': 'PNG'}
    key = "PNG"


@pyrat.docstringfrom(PNG)
def png(*args, **kwargs):
    PNG(*args, **kwargs).run(*args, **kwargs)


class TIFF(Pixmap):
    """
    TIFF format writer
    """
    gui = {'menu': 'File|Export to pixmap', 'entry': 'TIFF'}
    key = "TIFF"


@pyrat.docstringfrom(TIFF)
def tiff(*args, **kwargs):
    TIFF(*args, **kwargs).run(*args, **kwargs)


class PDF(Pixmap):
    """
    PDF format writer
    """
    gui = {'menu': 'File|Export to pixmap', 'entry': 'PDF'}
    key = "PDF"


@pyrat.docstringfrom(PDF)
def pdf(*args, **kwargs):
    PDF(*args, **kwargs).run(*args, **kwargs)


class EPS(Pixmap):
    """
    EPS format writer
    """
    gui = {'menu': 'File|Export to pixmap', 'entry': 'EPS'}
    key = "EPS"


@pyrat.docstringfrom(EPS)
def eps(*args, **kwargs):
    EPS(*args, **kwargs).run(*args, **kwargs)
