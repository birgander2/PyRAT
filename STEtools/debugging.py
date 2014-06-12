import numpy
import STEtools as STE

def mm(array):
    """
    Prints some simple statistics about an array. Helpful when hunting bugs...
    :author: Andreas Reigber    
    """
    print
    print 'Minimum     : ',numpy.min(array)
    print 'Maximum     : ',numpy.max(array)
    print 'Mean        : ',numpy.mean(array)
    print 'Stddev      : ',numpy.std(array)
   
def fftspeed(n):
    """
    Calculates the pseudo-relative speed of a fast(!) fourier transform
    
    :author: Andreas Reigber    
    :param n: length of the signal
    :type arg1: int
    :returns: a value which indicates the execution time of a fast fourier transform. More means longer...
   
    """
    fac = STE.primefactors(n)
    speed = sum(fac)*n
    return(speed)
    
def pleh(array):
    pass
    