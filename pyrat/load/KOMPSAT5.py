__author__ = 'kemp_mc'

import pyrat
import h5py
import logging
import copy
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError


class KOMPSAT5(pyrat.ImportWorker):
    """
    Import of Kompsat-5 satellite data.

    """
    gui = {'menu': 'File|Import spaceborne', 'entry': 'KOMPSAT-5'}

    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'something with *.h5 please'},
            {'var': 'crop', 'value': [0, 0, 0, 0]}]

    def __init__(self, *args, **kwargs):
        super(KOMPSAT5, self).__init__(*args, **kwargs)
        self.name = 'KOMPSAT5 IMPORT'
        if len(args) == 1:
            self.file = args[0]

    def reader(self, *args, **kwargs):
        # read data
        try:
            fileID = h5py.File(self.file, 'r')

        except:
            logging.error('Error: can not read file')
            return None, None

        B001 = fileID['S01']['B001']
        QLK = fileID['S01']['QLK']
        SBI = fileID['S01']['SBI']

        if len(SBI.shape) is 2:
            data = SBI[:, :]
        elif len(SBI.shape) is 3:
            re = SBI[:, :, 0]
            im = SBI[:, :, 1]
            data = re + 1j * im
        else:
            data = SBI

        # crop image
        if not tuple(self.crop) == (0, 0, 0, 0):
            self.crop[1] = data.shape[0] if self.crop[1] is 0 else self.crop[1]
            self.crop[3] = data.shape[1] if self.crop[3] is 0 else self.crop[3]
            data = data[self.crop[2]:self.crop[3], self.crop[0]:self.crop[1]]

        # read meta
        meta = {}
        metapath = self.file.replace('.h5', '_Aux.xml')
        try:
            tree = ET.parse(metapath)
            assert tree is not None

        except (AssertionError, ParseError, Exception):
            logging.error('XML file not found')

        else:
            meta['prf'] = [element.text for element in tree.iter('PRF')][0]
            meta['rs'] = [element.text for element in tree.iter('SamplingRate')][0]

        return data, meta

    @classmethod
    def guirun(cls, viewer):
        para_backup = copy.deepcopy(cls.para)  # keep a deep copy of the default parameters
        res = 1
        if len(cls.para) > 0:
            wid = pyrat.viewer.Dialogs.FlexInputDialog([cls.para[0]], parent=viewer, doc=cls.__doc__)
            res = wid.exec_()
        if res == 1:
            plugin = cls()  # instance with new parameters
            setattr(cls, 'para', para_backup)  # copy back the defaults
            viewer.statusBar.setMessage(message=' ' + plugin.name + ' ', colour='R')
            layers = plugin.run()
            del plugin
            viewer.statusBar.setMessage(message=' Ready ', colour='G')
            viewer.updateViewer()


@pyrat.docstringfrom(KOMPSAT5)
def kompsat5(*args, **kwargs):
    return KOMPSAT5(*args, **kwargs).run(*args, **kwargs)
