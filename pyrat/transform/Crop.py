import pyrat
from PyQt5 import QtCore, QtWidgets, QtGui
import copy
import numpy as np

class Crop(pyrat.FilterWorker):
    """
    Crop region (optionally using Qrubberband in viewer).
    Crop coordinates: [azimuth start, azimuth size, range start, range size]

    :author: Andreas Reigber
    """
    # todo: With only reading selected block...)
    # todo: Update rubberband when changing coordinates in GUI (add apply button...)

    gui = {'menu': 'Tools', 'entry': 'Crop region'}
    para = [{'var': 'crop', 'type': 'int', 'value': [0, 0, 0, 0], 'text': 'Crop region [pixel]',
             'subtext': ['y start', 'y size', 'x start', 'x size']}]

    def __init__(self, *args, **kwargs):
        super(Crop, self).__init__(*args, **kwargs)
        self.name = "CROP REGION"
        self.blockprocess = False

    def pre(self, *args, **kwargs):
        self.orig_height = pyrat.data.dshape[0]

    def filter(self, array, *args, **kwargs):
        block = list(self.crop)
        if block[0] < 0 or block[0] >= array.shape[-2]:
            block[0] = 0
        if block[2] < 0 or block[2] >= array.shape[-1]:
            block[2] = 0

        if block[1] == 0 or block[0] + block[1] > array.shape[-2]:
            block[1] = array.shape[-2] - block[0]
        if block[3] == 0 or block[2] + block[3] > array.shape[-1]:
            block[3] = array.shape[-1] - block[2]
        return array[..., block[0]:block[0] + block[1], block[2]:block[2] + block[3]]

    def post(self, *args, **kwargs):
        meta = pyrat.data.getAnnotation(layer=self.layer)
        if("geo_min_east" in meta and "orig_max_north" in meta and
            "geo_ps_east" in meta and "geo_ps_north" in meta):
            pyrat.data.setAnnotation({'geo_min_east':
                                        meta["geo_min_east"] +
                                        self.crop[2] * meta["geo_ps_east"],

                                      'geo_min_north':
                                        meta["geo_min_north"] +
                                        (self.orig_height - self.crop[1] - self.crop[0] + 2) * meta["geo_ps_north"]
                                    })

    @classmethod
    def guirun(cls, viewer):
        rubberband = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, viewer)

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        viewer.show_rubberband = True
        while viewer.show_rubberband is True:
            QtCore.QCoreApplication.processEvents()
        crop = viewer.rubberband.geometry()
        QtWidgets.QApplication.restoreOverrideCursor()

        wx = viewer.imageLabel.width()
        wy = viewer.imageLabel.height()
        win_ratio = wx / wy

        ix = viewer.box[1] - viewer.box[0]
        iy = viewer.box[3] - viewer.box[2]
        im_ratio = ix / iy

        if im_ratio >= win_ratio:  # width matches
            scale = ix / wx
        else:
            scale = iy / wy

        xs = int(ix / scale)
        ys = int(iy / scale)
        xo = (wx - xs) // 2
        yo = (wy - ys) // 2

        x1 = viewer.box[0] + int(scale * (crop.x() - xo))
        x2 = x1 + int(scale * crop.width())
        y1 = viewer.box[2] + int(scale * (crop.y() - yo))
        y2 = y1 + int(scale * crop.height())

        para_backup = copy.deepcopy(cls.para)  # keep a deep copy of the default parameters
        cls.para[0]['value'] = [y1, y2 - y1, x1, x2 - x1]

        wid = pyrat.viewer.Dialogs.FlexInputDialog(cls.para, parent=viewer, doc=cls.__doc__)
        res = wid.exec_()
        if res == 1:
            plugin = cls()  # instance with new parameters
            setattr(cls, 'para', para_backup)  # copy back the defaults
            viewer.statusBar.setMessage(message=' ' + plugin.name + ' ', colour='R')
            layers = plugin.run()
            del plugin
            viewer.statusBar.setMessage(message=' Ready ', colour='G')
            viewer.updateViewer()


@pyrat.docstringfrom(Crop)
def crop(*args, **kwargs):
    return Crop(*args, **kwargs).run(*args, **kwargs)


class Regionstat(pyrat.FilterWorker):
    """

    Performs distributed target analysis after cropping a set region:
    Crop coordinates: [azimuth start, azimuth size, range start, range size]

    :author: Joel Amao
    """


    gui = {'menu': 'SAR', 'entry': 'Region statistics'}
    para = [{'var': 'crop', 'type': 'int', 'value': [0, 0, 0, 0], 'text': 'Crop region [pixel]',
             'subtext': ['y start', 'y size', 'x start', 'x size']},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude','intensity'],
             'text': 'SAR data type'}]

    def __init__(self, *args, **kwargs):
        super(Regionstat, self).__init__(*args, **kwargs)
        self.name = "REGION STATISTICS"
        self.blockprocess = False

    def pre(self, *args, **kwargs):
        self.orig_height = pyrat.data.dshape[0]

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']

        block = list(self.crop)
        if block[0] < 0 or block[0] >= array.shape[-2]:
            block[0] = 0
        if block[2] < 0 or block[2] >= array.shape[-1]:
            block[2] = 0

        if block[1] == 0 or block[0] + block[1] > array.shape[-2]:
            block[1] = array.shape[-2] - block[0]
        if block[3] == 0 or block[2] + block[3] > array.shape[-1]:
            block[3] = array.shape[-1] - block[2]

        if array.ndim >= 3:
            nch = np.shape(array)[0]
        else:
            array = [np.newaxis, ...]
            nch = 1


        crop = array[..., block[0]:block[0] + block[1], block[2]:block[2] + block[3]]
        enl = 0
        for i in range(nch):
            if self.type == "amplitude":
                enl1 = np.abs(crop[i, ...])
                beta = (np.sqrt(np.mean((enl1 - (np.mean(enl1))) ** 2))) / np.mean(enl1)
                enl +=  (0.5227 / beta)** 2
            else:
                enl1 = np.abs(crop[i, i, ...])
                beta = (np.sqrt(np.mean((enl1 - (np.mean(enl1))) ** 2))) / np.mean(enl1)
                enl += 1 / (beta ** 2)
        enl /= nch
        print('\n  Data type: ' + str(self.type))
        print('  Number of pixels total: ' + str(np.prod(array.shape)))
        print('  Number of pixels considered: ' + str(np.prod(crop.shape)))
        print('  Ratio of considered pixels: ' + str(100*(np.prod(crop.shape)/
                                                      np.prod(array.shape))) + '%')
        print('\n  Minimum: ' + str(np.min(np.abs(crop))))
        print('  Maximum: ' + str(np.max(np.abs(crop))))
        print('\n  Mean: ' +str(np.mean(np.abs(crop))))
        print('  Median: ' +str(np.median(np.abs(crop))))
        print('  Std-Dev: ' +str(np.std(np.abs(crop))))
        print('  Equivalent Number of Looks: ' +str(enl))
        return crop

    def post(self, *args, **kwargs):
        meta = pyrat.data.getAnnotation(layer=self.layer)
        if("geo_min_east" in meta and "orig_max_north" in meta and
            "geo_ps_east" in meta and "geo_ps_north" in meta):
            pyrat.data.setAnnotation({'geo_min_east':
                                        meta["geo_min_east"] +
                                        self.crop[2] * meta["geo_ps_east"],

                                      'geo_min_north':
                                        meta["geo_min_north"] +
                                        (self.orig_height - self.crop[1] - self.crop[0] + 2) * meta["geo_ps_north"]
                                    })

    @classmethod
    def guirun(cls, viewer):
        rubberband = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, viewer)

        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        viewer.show_rubberband = True
        while viewer.show_rubberband is True:
            QtCore.QCoreApplication.processEvents()
        crop = viewer.rubberband.geometry()
        QtWidgets.QApplication.restoreOverrideCursor()

        wx = viewer.imageLabel.width()
        wy = viewer.imageLabel.height()
        win_ratio = wx / wy

        ix = viewer.box[1] - viewer.box[0]
        iy = viewer.box[3] - viewer.box[2]
        im_ratio = ix / iy

        if im_ratio >= win_ratio:  # width matches
            scale = ix / wx
        else:
            scale = iy / wy

        xs = int(ix / scale)
        ys = int(iy / scale)
        xo = (wx - xs) // 2
        yo = (wy - ys) // 2

        x1 = viewer.box[0] + int(scale * (crop.x() - xo))
        x2 = x1 + int(scale * crop.width())
        y1 = viewer.box[2] + int(scale * (crop.y() - yo))
        y2 = y1 + int(scale * crop.height())

        para_backup = copy.deepcopy(cls.para)  # keep a deep copy of the default parameters
        cls.para[0]['value'] = [y1, y2 - y1, x1, x2 - x1]

        wid = pyrat.viewer.Dialogs.FlexInputDialog(cls.para, parent=viewer, doc=cls.__doc__)
        res = wid.exec_()
        if res == 1:
            plugin = cls()  # instance with new parameters
            setattr(cls, 'para', para_backup)  # copy back the defaults
            viewer.statusBar.setMessage(message=' ' + plugin.name + ' ', colour='R')
            layers = plugin.run()
            del plugin
            viewer.statusBar.setMessage(message=' Ready ', colour='G')
            viewer.updateViewer()


@pyrat.docstringfrom(Regionstat)
def regionstat(*args, **kwargs):
    return Regionstat(*args, **kwargs).run(*args, **kwargs)
