import pyrat
from PyQt5 import QtCore, QtWidgets, QtGui
import copy


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
