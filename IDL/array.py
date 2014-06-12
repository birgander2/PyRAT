import numpy

def findgen(n, step=1):
    """
    Imitates IDL's findgen function. 
    """
    return numpy.arange(0, n, step, dtype='float32')
    
def dindgen(n, step=1):
    """
    Imitates IDL's dindgen function. 
    """
    return numpy.arange(0, n, step, dtype='float64')

def indgen(n, step=1):
    """
    Imitates IDL's indgen function. 
    """
    return numpy.arange(0, n, step, dtype='int16')
    
def lindgen(n, step=1):
    """
    Imitates IDL's lindgen function. 
    """
    return numpy.arange(0, n, step, dtype='int32')
    
def cindgen(n, step=1):
    """
    Imitates IDL's cindgen function. 
    """
    return numpy.arange(0, n, step, dtype='complex64')
    
def fltarr(*args):
    """
    Imitates IDL's fltarr function. 
    """
    return numpy.zeros(args,dtype='float32')
    
def complexarr(*args):
    """
    Imitates IDL's complexarr function. 
    """
    return numpy.zeros(args,dtype='complex64')
    
def intarr(*args):
    """
    Imitates IDL's intarr function. 
    """
    return numpy.zeros(args,dtype='int16')
    
def lonarr(*args):
    """
    Imitates IDL's lonarr function. 
    """
    return numpy.zeros(args,dtype='int32')
    
    
    
    
    
    
    
    
    
    