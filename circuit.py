# -*- coding: utf-8 -*-
"""
Class "Circuit" defines simple one port object/circuit having
frequency (single point or array), impedance, and possibly a list
of components that constitute the circuit.
"""


import numpy as np
import scipy
import matplotlib.pyplot as plt
import skrf

class Circuit:
    
    def __init__(self, f, Z, Z0=50, draw_smith=True):
        self.Z = Z
        self.f = f
        self.Z0 = Z0
        self.draw_smith=draw_smith
        self.components = [('start', Z)]
        

        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            skrf.plotting.plot_smith(self.refl(self.Z), marker='o', color='k', linestyle=None,  x_label='',y_label='', title='Smith chart, matching network', label='start')
        #self.components is an array of tuples containing the history of the circuit,
        #i.e., all components attached to the circuit

    def __str__(self):
        if not isinstance(self.f, np.ndarray):
            freq, multi = find_si_prefix(self.f)
            ret = "{:.2f}+j{:2f} ohm at {:.2f} {:s}Hz:\n".format(self.Z.real, self.Z.imag, freq, multi)
        else:
            return "not implemented for freq range"
        
        for component in self.components:
            prefix = find_si_prefix(component[1], False)
            ct = '' #component type
            cc = '' #component connection
            unit = ''
            n = component[0]
            if n.find('cap') >= 0:
                ct = 'capacitor'
                unit = 'F'
            elif n.find('ind') >= 0:
                ct = 'inductor'
                unit = 'H'
            elif n.find('line') >= 0:
                ct = 'line'
                unit = '\pi'
            
            if n.find('ser') >= 0:
                cc = 'Series'
            elif n.find('par') >= 0:
                cc = 'Parallel' 
            elif n.find('start') >= 0:
                cc = 'Start'            
            
            ret += '{:s} {:s}: {:.2f} {:s}{:s}\n'.format(cc, ct, prefix[0], prefix[1], unit)
        
        
        
        return ret
            
    def cap(self, C):
        '''
        impedance of a capacitor at frequency of this circuit
        '''
        return cap(self.f, C)
        
    def ind(self, L):
        return ind(self.f, L)
        
        
    def refl(self, Z=None):
        '''
        returns the reflection coefficient of the present circuit
        '''
        if Z is None:
            return refl(self.Z,self.Z0)
        else:
            return refl(Z, self.Z0)
        
    def ser(self, Z):
        self.components.append(('ser',Z))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = np.logspace(np.log10(Z/1000), np.log10(Z), 100)
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance

            scaled_value = find_si_prefix(Z)
            label = 'ser {:2g} {:s}$\Omega$'.format(scaled_value[0], scaled_value[1])
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z += Z #after plotting, add series capacitance to the current impedance

    def par(self, Z):
        self.components.append(('par',Z))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = np.logspace(np.log10(Z*1000), np.log10(Z), 100)
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = par(self.Z, tmp)
            #the array transformed into impedance

            scaled_value = find_si_prefix(Z)
            label = 'par {:2g} {:s}$\Omega$'.format(scaled_value[0], scaled_value[1])
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = par(Z, self.Z) #after plotting, add series capacitance to the current impedance    
        
   
    def sercap(self, C):
        '''
        Set a capacitor in series with the present circuit
        '''
        self.components.append(('sercap',C))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.cap(np.logspace(np.log10(C*1000), np.log10(C), 100)) 
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(C)
            label = 'ser {:2g} {:s}F'.format(scaled_value[0], scaled_value[1])
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = ser(self.Z, self.cap(C)) #after plotting, add series capacitance to the current impedance
        
    def serind(self, L):
        '''
        Set an inductor in series with the present circuit
        '''
        self.components.append(('serind',L))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.ind(np.logspace(np.log10(L/1000), np.log10(L), 100)) 
            #impedance of a capacitance array, starting from very small inductance
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(L)
            label = 'ser {:2g} {:s}H'.format(scaled_value[0], scaled_value[1])
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = ser(self.Z, self.ind(L)) #after plotting, add series inductance to the current impedance
        
    def parcap(self, C):
        '''
        Set a capacitor parallel to the present circuit
        '''
        self.components.append(('parcap',C))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.cap(np.logspace(np.log10(C/1000), np.log10(C), 100)) 
            #impedance of a capacitance array, starting from very small cap
            tmp = par(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(C)
            label = 'ser {:2g} {:s}F'.format(scaled_value[0], scaled_value[1])
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = par(self.Z, self.cap(C)) #after plotting, add series capacitance to the current impedance
        
    def parind(self, L):
        '''
        Set an inductor parallel to the present circuit
        '''
        self.components.append(('parind',L))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.ind(np.logspace(np.log10(L*1000), np.log10(L), 100)) 
            #impedance of a capacitance array, starting from very large inductance
            #(large series capacitance == small change in impedance)
            tmp = par(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(L)
            label = 'par {:2g} {:s}H'.format(scaled_value[0], scaled_value[1])
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = par(self.Z, self.ind(L)) #after plotting, add series inductance to the current impedance

    def serline(self, kl, Z0=None):
        '''
        perform impedance as seen through a transmission line
        '''
        if Z0==None:
            Z0=self.Z0

        self.components.append(('serline',kl, Z0))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = line(self.Z, np.linspace(0, kl, 100), Z0) 
            label = 'line {:2g}'.format(kl)
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = line(self.Z, kl, Z0) 

    def parline(self, kl, Zl=np.inf, Z0=None):
        '''
        add a transmission line parallel to the current impedance (defaults to open self.Z0 ohm line)
        '''
        if Z0==None:
            Z0=self.Z0

        self.components.append(('parline',kl, Zl, Z0))
        if np.size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = par(self.Z, line(Zl, np.linspace(0, kl, 100), Z0))
            label = 'par line {:2g}'.format(kl)
            skrf.plotting.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)

        self.Z = par(self.Z0, line(Zl, kl, Z0)) 
        
    def db(self, noextrastuff=False):
        '''
        For a point frequency return reflection coefficient, for a frequency range
        plot the response curve. noextrastuff means the labels and grids
        '''
        if np.size(self.f)==1:
            return 20*np.log10(abs(self.refl()))
        elif noextrastuff:
            plt.plot(self.f, 20*np.log10(abs(self.refl())))
        else:
            plt.figure()
            plt.plot(self.f, 20*np.log10(np.abs(self.refl())))
            plt.ylabel('Reflection coefficient, dB')
            plt.xlabel('Frequency, Hz')
            plt.grid()
            
    def smith(self, annotations=True,smith_r=1, chart_type='z', x_label='',
            y_label='', title='Smith chart, frequency', show_legend=True,
            axis='equal', ax=None, force_chart = False, *args, **kwargs):
            '''
            plots the current mathcing network on the smith chart as the function of frequency
            '''
            if np.size(self.f)>1:
                plt.figure()
                skrf.plotting.plot_smith(self.refl(), smith_r=smith_r, chart_type=chart_type, x_label=x_label,
                y_label=y_label, title=title, show_legend=show_legend,
                axis=axis, ax=ax, force_chart = force_chart, *args, **kwargs)
                
                if annotations:
                    xy=(np.real(self.refl()[0]),np.imag(self.refl()[0]))
                    plt.annotate('{:.2e}'.format(self.f[0]) , xy=xy,xycoords='data', xytext=(xy[0]/ np.abs(xy[0]), xy[1]/ np.abs(xy[1])), textcoords='data', arrowprops=dict(arrowstyle="->")).draggable()
                    xy=(np.real(self.refl()[-1]),np.imag(self.refl()[-1]))
                    plt.annotate('{:.2e}'.format(self.f[-1]) , xy=xy,xycoords='data', xytext=(xy[0]/ np.abs(xy[0]), xy[1]/ np.abs(xy[1])), textcoords='data', arrowprops=dict(arrowstyle="->")).draggable()
                    
                    ind = np.argmin(np.abs(self.refl()))
                    xy=(np.real(self.refl()[ind]),np.imag(self.refl()[ind]))
                    plt.annotate('{:.2e}\n{:.1f} dB'.format(self.f[ind], 20*np.log10(np.abs(self.refl()[ind]))) , xy=xy,xycoords='data', xytext=(xy[0]/ np.abs(xy[0])+0.2, xy[1]/ np.abs(xy[1])-0.2), textcoords='data', arrowprops=dict(arrowstyle="->")).draggable()

'''
After this point the file defines basic equations used in Circuit class
'''

def cap(f, C):
    '''
    impedance of a capacitor
    '''
    if isinstance(C, np.ndarray) or isinstance(f, np.ndarray):
        return 1/(1j*2*np.pi*f*C)
    elif np.iscomplex(C) or f<0 or C<0:
        raise ValueError #bullshit numbers given
    else:
        return 1/(1j*2*np.pi*f*C)


def ind(f, L):
    '''
    impedance of a inductor
    '''
    
    #TODO: check that input is valid. 
    return 1j*2*np.pi*f*L

def ser(Z1, Z2, Z3=0, Z4=0, Z5=0, Z6=0, Z7=0, Z8=0, Z9=0, Z10=0):
    '''
    impedance of a series of compoents
    '''
    return Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10

def par(Z1, Z2, Z3=np.inf, Z4=np.inf, Z5=np.inf, Z6=np.inf, Z7=np.inf, Z8=np.inf, Z9=np.inf, Z10=np.inf):
    '''
    impedance of parallel compoents
    '''
    return 1/(1/Z1+1/Z2+1/Z3+1/Z4+1/Z5+1/Z6+1/Z7+1/Z8+1/Z9+1/Z10)

def refl(Z, Z0=50):
    if isinstance(Z, np.ndarray):
        return (Z-Z0)/(Z+Z0)
    elif Z==np.inf and Z0<np.inf:
        return 1
    else:
        return (Z-Z0)/(Z+Z0)
    
def find_si_prefix(x, use_latex = True):
    '''
    return (numeric value , si-prefix)
    
    >>> find_si_prefix(5000)
    (5.0, 'k')
    '''
    if use_latex:
        prefixes = ['f', 'p', 'n', '$\mu$', 'm', '', 'k','M','G']
    else:
        prefixes = ['f', 'p', 'n', 'u', 'm', '', 'k','M','G']
    multipliers = np.asarray([1e-15,1e-12, 1e-9, 1e-6,1e-3,1,1e3,1e6,1e9])
    for n, k in enumerate(x/multipliers):
        if 0.1<k<=(100 + 1e-10):
            break #loop the list until we find a prefix that sets the k to be in nice range
    return (k, prefixes[n])
  
def line(Z, kl, Z0=50):
    #from scipy.special import cotdg
    if Z != np.inf:
        Z = Z0 * (Z+1j*Z0*np.tan(kl))/(Z0 + 1j* Z*np.tan(kl))
    else:
        np.seterr(divide='ignore') #lazy solution to situation kl=0 == divide by error
        Z = -1j*Z0*np.true_divide(1,np.tan(kl))# if np.tan(kl)==0 else 1j*np.inf
        #Z = -1j*Z0*np.cos(kl)/np.sin(kl)
        np.seterr(divide='warn')
    return Z     
    

def db(Z, Z0=50):
    if Z==Z0: #avoid log(0)
        return -np.inf
    else:
        return 20*np.log10(abs(refl(Z, Z0)))
        
        
if __name__ == '__main__':
    plt.ion()
    c = Circuit(64e6, 30)
    c.sercap(30e-12)
    c.parind(100e-9)


