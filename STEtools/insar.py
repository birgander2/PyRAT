import numpy
import IDL

def coreg(arr1, arr2, sub=False):
    """
    Returns the coregistration offset between two complex arrays

    :author: Andreas Reigber
    :param arr1: The master image as (2D numpy.ndarray)
    :type arr1: complex
    :param arr2: The slave image as (2D numpy.ndarray)
    :type arr2: complex
    :param sub=True: Calculate subpixel offset.
    :type arr2: boolean
    :returns: tuple containing slave image offset in y and x
    """
    s = arr1.shape
    out1 = numpy.fft.rfft2(abs(arr1)).astype('complex64')
    out2 = numpy.fft.rfft2(abs(arr2)).astype('complex64')
    out1 *= numpy.conj(out2)
    out2 = numpy.fft.irfft2(out1).astype('float32')
    pos  = out2.argmax()
    xoff = pos % s[1]
    yoff = pos / s[1]
    if xoff >= s[1] // 2: xoff = xoff - s[1]
    if yoff >= s[0] // 2: yoff = yoff - s[0]
    if sub == True:           # explicit polynominal expansion of the peak
        ysub = 0.5*(out2[yoff-1,xoff]-out2[yoff+1,xoff])/(out2[yoff-1,xoff]+out2[yoff+1,xoff]-2*out2[yoff,xoff])
        xsub = 0.5*(out2[yoff,xoff-1]-out2[yoff,xoff+1])/(out2[yoff,xoff-1]+out2[yoff,xoff+1]-2*out2[yoff,xoff])
        yoff += numpy.round(ysub,2)
        xoff += numpy.round(xsub,2)
    return yoff.astype('float32'),xoff.astype('float32')

def coherence(arr1, arr2, box=(7,7)):
    """
    Returns the interferometric coherence between two slc images
    :author: Andreas Reigber
    :param arr1: The master image as (2D numpy.ndarray)
    :type arr1: complex
    :param arr2: The slave image as (2D numpy.ndarray)
    :type arr2: complex
    :param box=(7,7): Size of the local box used for coherence calculation
    :type arr2: integer
    :returns: interferometric coherence (float)
    """
    coh = numpy.abs(IDL.smooth(arr1*numpy.conj(arr2),box)/numpy.sqrt(IDL.smooth(arr1*numpy.conj(arr1),box)*IDL.smooth(arr2*numpy.conj(arr2),box)))
    return numpy.clip(numpy.nan_to_num(coh),0.0,1.0)  # get rid of numerical inaccuracies!
    
def wrap(array):
    """
    Rewraps an interferometric phase back to the -pi/pi interval
    :author: Andreas Reigber
    :param array: The unwrapped phase (2D numpy.ndarray)
    :type arr1: float
    :returns: the wrapped phase (float)
    """
    return numpy.angle(numpy.exp(1j*array))
    
def unwrap1d(arr):
    """
    Very simple unwrapping of 1D phase functions.
    :author: Andreas Reigber
    :param arr: The wrapped phase (1D numpy.ndarray)
    :type arr: float
    :returns: the (hopefully) unwrapped array
    """
    oarr = arr.copy()
    for idx in range(1,oarr.size):
        oarr[idx] = oarr[idx-1] + numpy.angle(numpy.exp(1j*(arr[idx]-arr[idx-1])))
    return oarr
    
def flattenpha(arr):
    """
    Removes linear phase components from a 2D InSAR phase
    :author: Andreas Reigber
    :param array: The wrapped phase (2D numpy.ndarray)
    :type arr: float
    :returns: the flattened phase (float)
    """
    s = arr.shape
    oarr = numpy.fft.fft2(numpy.exp(1j*arr))
    pos  = abs(oarr).argmax()
    oarr = numpy.roll(oarr, -pos / s[1], axis = 0)
    oarr = numpy.roll(oarr, -pos % s[1], axis = 1)
    return numpy.angle(numpy.fft.ifft2(oarr))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    