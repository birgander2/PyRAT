import numpy

def shift(arr,*args):
    """
    Imitates IDL's shift function. 
    WARNING 1: Indexing in python is of opposite order (C-type) than in IDL (Fortran-type) and so it is here!
    WARNING 2: Unluckily much slower than in IDL. numpy.roll seems not to be formula-1 quality.

    :author: Andreas Reigber
    """
    oarr = arr.copy()
    for index, offset in enumerate(args):
        oarr = numpy.roll(oarr, offset, axis = index)
    return oarr

    