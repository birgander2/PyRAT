import matplotlib
matplotlib.use('Qt5Agg')
import pylab as plt

import pyrat
import numpy as np
import pylab as py


class Unweight(pyrat.FilterWorker):
    """
    Removes hamming (or any other) window function in frequency domain
    Author: J. Fischer, DLR/HR
    Date:   24/10/2014
    Docu:   Fischer,Pupeza,Scheiber: Sidelobe Suppression Using the SVA Method for SAR Images and Sounding Radars, EUSAR Proceedings, 2006
    """

    # Comment: Needs reworking before getting accepted for the GUI
    # gui = {'menu': 'SAR|Sidelobe control', 'entry': 'Remove spectral weights (UNWEIGHT)'}
    para = [
        {'var': 'ov', 'value': [2, 2], 'type': 'int', 'range': [1, 10], 'text': 'Desired Oversampling', 'subtext': ['azimuth', 'range']},
        {'var': 'dim', 'value': 'both', 'type': 'list', 'range': ['range','azimuth','both'], 'text': 'Dimension'}
        ]

    def __init__(self, *args, **kwargs):
        super(Unweight, self).__init__(*args, **kwargs)
        self.name = "UNWEIGHT FILTER"

        #if 'ov' not in self.__dict__:
            #self.ov = [2, 2]
        #elif isinstance(self.ov, int):
            #self.ov = [self.ov] * 2
     
    def filter(self, array, *args, **kwargs):

        # Processed bandwith in percentages
        #------------------------------------------------------------------------
        sys = kwargs["meta"]        
        
        bw_proc_az = 1.33 *  sys['v0']     / sys['res_az']
        bw_proc_rg = 1.33 * (sys['c0']/2.) / sys['res_rg']
        
        percentage_az = 100.* bw_proc_az / (sys['prf']/sys['pre_az'])
        percentage_rg = 100.* bw_proc_rg /  sys['rsf']
        
        bw  = [percentage_az, percentage_rg]
        print("  bw=["+str(bw[0])+","+str(bw[1])+"] ... ")
        #------------------------------------------------------------------------        

        if array.ndim == 2 and np.iscomplexobj(array):
            return self.unweight2d(array, self.ov, bw)
        
        if array.ndim == 3 and np.iscomplexobj(array):
            p = array.shape
            for k in range(0,p[0]):
                array_temp = self.unweight2d(array[k,:,:], self.ov, bw)
                if k == 0:
                    s = array_temp.shape
                    array_new = np.empty((p[0],s[0],s[1]), dtype='complex64')
                array_new[k,:,:] = array_temp
            return array_new
        
        else:
            print("  ERROR: Bad input.")
            return None

        
    def unweight2d(self, data, ov, bw):
        
        print("  Applying UNWEIGHT with ov=["+str(ov[0])+","+str(ov[1])+"] ... ")
        
        n0 = data.shape[0]
        n1 = data.shape[1]

        # Percentage -> Pixels, secure having integers
        #------------------------------------------------------------
        bw0 = int(np.floor( (bw[0]*n0/100.) /2) * 2)  # even integer
        bw1 = int(np.floor( (bw[1]*n1/100.) /2) * 2)  # even integer
        
        ov0 = int(np.floor(ov[0]))  
        ov1 = int(np.floor(ov[1]))  
        #------------------------------------------------------------
        
        spec0 = py.fft2(data)
        spec0 = np.roll(spec0,data.shape[0]/2, axis=0)
        spec0 = np.roll(spec0,data.shape[1]/2, axis=1)
        

        # Hamming at processed bandwidth
        #-----------------------------------------------
        if 1:
            t0 = np.arange(bw0)/float(bw0)
            t1 = np.arange(bw1)/float(bw1)
            hamming0 = 0.54-0.46*np.cos(2*np.pi*t0)        
            hamming1 = 0.54-0.46*np.cos(2*np.pi*t1)
        
            unham0 = np.zeros(n0)
            unham1 = np.zeros(n1)
            
            unham0[n0/2-bw0/2:n0/2+bw0/2] = hamming0
            unham1[n1/2-bw1/2:n1/2+bw1/2] = hamming1
        #-----------------------------------------------


        spec0_profile0 = np.abs(spec0).mean(axis=1);  maxv0 = 0.95 * np.max(np.abs(spec0_profile0))
        spec0_profile1 = np.abs(spec0).mean(axis=0);  maxv1 = 0.95 * np.max(np.abs(spec0_profile1))


        # Remove doppler shift and range spectrum shift
        #------------------------------------------------------------------------------------------------------
        corr0 = np.abs( py.ifft( py.fft(np.abs(spec0_profile0)) * np.conj(py.fft(np.abs(unham0))) ))
        corr1 = np.abs( py.ifft( py.fft(np.abs(spec0_profile1)) * np.conj(py.fft(np.abs(unham1))) ))
        
        peak0 = np.where(abs(corr0) == np.max(abs(corr0)));  off0 = n0 - peak0[0]  
        peak1 = np.where(abs(corr1) == np.max(abs(corr1)));  off1 = n1 - peak1[0]  
        
        spec0 = np.roll(spec0, off0, axis=0)
        spec0 = np.roll(spec0, off1, axis=1)
        
        spec0_profile0 = np.abs(spec0).mean(axis=1);  maxv0 = 0.95 * np.max(np.abs(spec0_profile0))
        spec0_profile1 = np.abs(spec0).mean(axis=0);  maxv1 = 0.95 * np.max(np.abs(spec0_profile1))        
        #------------------------------------------------------------------------------------------------------
            
        
        # Replace Unhamming filter by profile filter
        #------------------------------------------------------------------------------
        if 1:
            unham0 = self.smooth(spec0_profile0 / maxv0, window_len=11)
            unham1 = self.smooth(spec0_profile1 / maxv1, window_len=11)
        #------------------------------------------------------------------------------


        # Show profiles
        #------------------------------------------------    
        show_plots = False
        if show_plots:
            plt.plot(spec0_profile0,'k-', lw=1, color='blue')
            plt.show()
            
            plt.plot(spec0_profile1,'k-', lw=1, color='red')
            plt.show()        
        #------------------------------------------------    
        
        
        # Compare profiles to hamming filter
        #----------------------------------------------------------------------------------------------------------------
        if show_plots:
            plt.plot(spec0_profile0,'k-', lw=1, color='blue')
            plt.plot(self.smooth(spec0_profile0,window_len=21),'k-', lw=1, color='green')
            plt.plot(maxv0 * unham0,'k--', lw=1, color='red')
            plt.show()

            plt.plot(spec0_profile1,'k-', lw=1, color='blue')
            plt.plot(self.smooth(spec0_profile1,window_len=21),'k-', lw=1, color='green')
            plt.plot(maxv1 * unham1,'k--', lw=1, color='red')
            plt.show()
            
            plt.plot(spec0_profile0[n0/2-bw0/2:n0/2+bw0/2],'k-', lw=1, color='blue')
            plt.plot(maxv0 * unham0[n0/2-bw0/2:n0/2+bw0/2],'k-', lw=1, color='red')
            plt.show()
            
            plt.plot(spec0_profile1[n1/2-bw1/2:n1/2+bw1/2], 'k-', lw=1, color='blue')
            plt.plot(maxv1 * unham1[n1/2-bw1/2:n1/2+bw1/2], 'k-', lw=1, color='red')
            plt.show()

            plt.plot(spec0_profile0[n0/2-bw0/2:n0/2+bw0/2] / (maxv0 * unham0[n0/2-bw0/2:n0/2+bw0/2]),'k-', lw=1, color='blue')
            plt.show()
            
            plt.plot(spec0_profile1[n1/2-bw1/2:n1/2+bw1/2] / (maxv1 * unham1[n1/2-bw1/2:n1/2+bw1/2]),'k-', lw=1, color='blue')
            plt.show()        
        #----------------------------------------------------------------------------------------------------------------
        

        # Unhamming
        #------------------------------------------------------------------
        #print "  mean ..."+str(np.mean(abs(spec0)))
        #print "    Unhamming ..."
        for k in range(0,n1):
            spec0[n0/2-bw0/2:n0/2+bw0/2,k] /= unham0[n0/2-bw0/2:n0/2+bw0/2]   # range (y)
        for k in range(0,n0):                                       
            spec0[k,n1/2-bw1/2:n1/2+bw1/2] /= unham1[n1/2-bw1/2:n1/2+bw1/2]   # azimuth (x)
        #print "    Unhamming done."    
        #print "  mean ..."+str(np.mean(abs(spec0)))
        #------------------------------------------------------------------
        
        
        # Show profiles
        #------------------------------------------------
        if show_plots:            
            abs_spec0 = np.abs(spec0)
            spec0_profile0 = abs_spec0.sum(axis=1)
            spec0_profile1 = abs_spec0.sum(axis=0)
            
            plt.plot(spec0_profile0,'k-', lw=1, color='blue')
            plt.show()
            
            plt.plot(spec0_profile1,'k-', lw=1, color='red')
            plt.show()        
        #------------------------------------------------


        spec0 = np.roll(-spec0,data.shape[0]/2, axis=0)
        spec0 = np.roll(-spec0,data.shape[1]/2, axis=1)
        
        
        # Zero padding
        #--------------------------------------------------------------------------------------
        n0 = spec0.shape[0]
        n1 = spec0.shape[1]
        zeros0 =  np.zeros((bw0 * (ov0-1),n1),float)  + 1j * np.zeros((bw0 * (ov0-1),n1),float)
        
        spec1 = np.concatenate( (spec0[0:bw0/2,:], zeros0, spec0[-bw0/2:,:]), axis=0)  *  ov0
        
        n0 = spec1.shape[0]
        n1 = spec1.shape[1]
        zeros1 =  np.zeros((n0,bw1 * (ov1-1)),float)  + 1j * np.zeros((n0,bw1 * (ov1-1)),float)
        
        spec2 = np.concatenate( (spec1[:,0:bw1/2], zeros1, spec1[:,-bw1/2:]), axis=1)  *  ov1
        #--------------------------------------------------------------------------------------
        
        
        # Show zeros padding results
        #--------------------------------------------------------------------------------
        '''
        plt.imshow(np.abs(spec0), origin='lower', interpolation='none', cmap=plt.cm.BuGn)
        plt.show()
        
        plt.imshow(np.abs(spec1), origin='lower', interpolation='none', cmap=plt.cm.BuGn)
        plt.show()
        
        plt.imshow(np.abs(spec2), origin='lower', interpolation='none', cmap=plt.cm.BuGn)
        plt.show()
        '''
        #--------------------------------------------------------------------------------
        
        
        data = py.ifft2(spec2)

        print("  Applying UNWEIGHT with ov=["+str(ov[0])+","+str(ov[1])+"] done. ")
        
        return data


    def smooth(self, x,window_len=11,window='hanning'):
        
        if x.ndim != 1:
                raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
                raise ValueError("Input vector needs to be bigger than window size.")
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
            
        s = np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:  
            w=eval('np.'+window+'(window_len)')
        
        y = np.convolve(w/w.sum(),s,mode='same')
        
        return y[window_len:-window_len+1]


@pyrat.docstringfrom(Unweight)
def unweight(*args, **kwargs):
    return Unweight(*args, **kwargs).run(**kwargs)
        
        
