from pyrat.load.tools import RatFile


def srat(filename, array, **kwargs):
    """
    Writes a numpy ndarray into a RAT file
    """
    File = RatFile(filename)
    File.write(array, **kwargs)
    del File
