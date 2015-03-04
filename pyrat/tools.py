import sys
import logging
import pyrat

class ProgressBar():
    """
    Simple progress bar for the command line.

    :author: Andreas Reigber
    """

    def __init__(self, message, max, width=60):
        terminal_width = width
        self.message = message.ljust(20)
        if terminal_width < 50:
            self.message = self.message[:terminal_width / 2]
        self.width = terminal_width - len(self.message) - 10
        self.max = max
        if hasattr(pyrat, "app"):
            pyrat.app.statusBar.progressbar.reset()

    def __del__(self):
        if hasattr(pyrat, "app"):
            pyrat.app.statusBar.progressbar.setValue(0)
        print()

    def update(self, val):
        """
        Updates the progress bar to the given value of progress
        """
        percent = min(float(val) / self.max, 1.0)
        if hasattr(pyrat, "app"):
            pyrat.app.statusBar.progressbar.setValue(int(percent * 100))

        hashes = '#' * int(round(percent * self.width))
        spaces = ' ' * (self.width - len(hashes))
        retline = "\r" if sys.stdout.isatty() else ""
        if sys.stdout.isatty() or val == 0:
            sys.stdout.write(retline + self.message + ": [{0}] {1}%".format(hashes + spaces, int(round(percent * 100))))
            sys.stdout.flush()


def deshape(shape):
    """
    Extracts layer shape and data shape from a numpy ndarray shape
    :returns: lshape, dshape
    """
    lshape = (1,)
    dshape = tuple(shape)
    if len(shape) == 2:                                               # normal data
        pass
    elif len(shape) == 3:                                             # vector data
        lshape = (shape[0],)                                          # layer shape
        dshape = tuple(shape[1:])                                     # data shape
    elif len(shape) == 4:                                             # matrix data
        lshape = (shape[0], shape[1])                                 # layer shape
        dshape = tuple(shape[2:])                                     # data shape
    else:
        logging.error('Something wrong with array dimensions!')
        lshape = False
        dshape = False
    return lshape, dshape


