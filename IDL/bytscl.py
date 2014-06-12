import numpy
def bytscl(img,start=False,end=False):
    """
    Imitates IDL's bytscl function. 
    
    :author: Andreas Reigber
    """
    if start == False: start = numpy.min(img)
    if   end == False:   end = numpy.max(img)
    return numpy.uint8(numpy.clip((img-start)/(end-start)*255,0,255))
   
