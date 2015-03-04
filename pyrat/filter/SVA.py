import pyrat
import numpy as np

class SVA(pyrat.FilterWorker):
    """
    Spatial Variant Apodisation (SVA) filter: Removes SINC function sidelobes
    Author: J. Fischer, DLR/HR
    Date:   24/10/2014
    Docu:   Fischer,Pupeza,Scheiber: Sidelobe Suppression Using the SVA Method for SAR Images and Sounding Radars, EUSAR Proceedings, 2006
    """

    # Comment: Needs reworking before getting accepted for the GUI
    # gui = {'menu': 'SAR|Sidelobe control', 'entry': 'Spatial Variant Apodization (SVA)'}
    para = [
        {'var': 'ov', 'value': [2, 2], 'type': 'int', 'range': [1, 10], 'text': 'Desired Oversampling', 'subtext': ['azimuth', 'range']}
        ]

    def __init__(self, *args, **kwargs):
        super(SVA, self).__init__(*args, **kwargs)
        self.name = "SVA FILTER"

        if 'ov' not in self.__dict__:
            self.ov = [2, 2]
        elif isinstance(self.ov, int):
            self.ov = [self.ov] * 2

        self.blockprocess = True
        self.blockoverlap = self.ov[0]  # = 1 pixel margin x azimuth oversampling factor
     
    def pre(self, *args, **kwargs):
        #pyrat.filter.unweight(ov=self.ov)
        pass
     
    def post(self, *args, **kwargs):
        pass

    def filter(self, array, *args, **kwargs):

        # 1-D Real Arrays
        if array.ndim == 1 and np.isrealobj(array):
            return self.svafilter(array, self.ov)
        
        # 1-D Complex Arrays
        if array.ndim == 1 and np.iscomplexobj(array):
            return self.sva1D(array, self.ov)
        
        # 2-D Complex Arrays
        if array.ndim == 2 and np.iscomplexobj(array):
            return self.sva2D(array, self.ov)
        
        # 3-D Complex Arrays
        if array.ndim == 3 and np.iscomplexobj(array):
            p = array.shape            
            for k in range(0,p[0]):
                array[k,:,:] = self.sva2D(array[k,:,:], self.ov)
            return array
        
        else:
            print("  ERROR: Bad input.")
            return None


    def sva2D(self, image0, ov):

        # SVA for complex valued 2D arrays
        
        print("  Applying SVA with ov=["+str(ov[0])+","+str(ov[1])+"] ... ")
        
        s = image0.shape
        image1 = np.empty((s[0],s[1]), dtype='complex64')
        image2 = np.empty((s[0],s[1]), dtype='complex64')
        
        for k in range(0,s[0]):
            image1[k,:] = self.sva1D( image0[k,:], ov[1] )  # SVA in range
        for k in range(0,s[1]):
            image2[:,k] = self.sva1D( image1[:,k], ov[0] )  # SVA in azimuth

        print("  Applying SVA with ov=["+str(ov[0])+","+str(ov[1])+"] done. ")

        return image2


    def sva1D(self, image, ov):

        # SVA for complex valued 1D arrays

        #k,d = get_kd(image)  # Note that n = k*d

        k = image.size /  ov
        d = ov

        # n = total signal length
        # k = signal bandwidth in [pixel], k = supp(spectrum), k is the 'actual signal length'
        # d = integer oversampling factor = number of pixels between sinc zero crossings, d is the redundancy factor
        
        return self.svafilter(image.real,k,d) + 1j * self.svafilter(image.imag,k,d)


    def svafilter(self, image,k,d):

        # SVA for real valued 1D arrays

        image_sva = np.empty(image.size,float)
        ind = np.arange(0,k) * d

        for off in range(0,d): # run through all "offsets" (between zero crossings)
            image_sva[ind+off] = self.svafilterkernel( image[ind+off], k)

        return image_sva


    def svafilterkernel(self, sig,k):

        # SVA filter kernel for real valued 1D arrays

        np.seterr(divide='ignore',invalid='ignore')
        wu = np.zeros(k)
        for m in range(1,k-1):
            wu[m]=-sig[m]/(sig[m-1]+sig[m+1])
        np.seterr(divide='warn',invalid='warn')
        
        a = np.zeros(k)
        for m in range(1,k-1):
            if wu[m] <= 0: 
                a[m] = sig[m] # mainlobe (preserved)
            elif (wu[m] >= 0.5):  #sidelobe (0 < wu[m] < 0.5) set to zero
                a[m]=sig[m]+0.5*(sig[m-1]+sig[m+1]) # mixed signal (nothing done)

        return a

        
def sva(*args, **kwargs):
    return SVA(*args, **kwargs).run(**kwargs)
        
        
