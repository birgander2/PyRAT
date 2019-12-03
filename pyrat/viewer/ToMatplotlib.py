import pyrat
import matplotlib.pyplot as plt


class ToMatPlotLib(pyrat.ExportWorker):
    """
    Exports the active layer to matplotlib
    """

    gui = {'menu': 'File',
           'entry': 'Export to matplotlib',
           'before': 'Export to raster'}

    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.name = "Export to matplotlib"

    def writer(self, array, *args, **kwargs):
        if type(array) is list:
            for layer in array:
                self.writer(layer)
        else:
            plt.figure()
            plt.imshow(array)
            plt.colorbar()
            plt.show()


@pyrat.docstringfrom(ToMatPlotLib)
def tomatplotlib(*args, **kwargs):
    return ToMatPlotLib(*args, **kwargs).run(*args, **kwargs)
