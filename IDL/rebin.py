import numpy
def rebin(arr, *shape, **kwargs):
    """
    Imitates IDL's rebin function. Allows also phase rebining.

    :author: Andreas Reigber
    """

    phase = False                                      # combination of *args and fixed keywords
    if 'phase' in kwargs: phase = kwargs['phase']      # works only in in python 3!
    
    if len(shape) == 1:                                # allows to pass shape as list/array or normal arguments
        if type(shape[0]) == int: 
            shape = [shape[0]]         
        else:
            shape = list(shape[0])
            
    oarr  = arr.copy()
    oshap = oarr.shape
    for d in range(arr.ndim):
        n = shape[d]
        N = oshap[d]
        if n < N:
            s = list(oarr.shape)
            s.insert(d+1,N//n)
            s[d] = n
            if phase == True:
                oarr = numpy.angle(numpy.exp(1j*oarr.reshape(s)).mean(d+1))
            else:
                oarr = oarr.reshape(s).mean(d+1)
        else:
            oarr = oarr.repeat(n//N, axis=d)
    return oarr
