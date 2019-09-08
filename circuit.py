# -*- coding: utf-8 -*-
"""
Class "Circuit" defines simple one port object/circuit having
frequency (single point or array), impedance, and possibly a list
of components that constitute the circuit.
"""

from __future__ import division
import numpy as np
from scipy import *
import matplotlib.pyplot as plt
import skrf

class Circuit:
    
    def __init__(self, f, Z, Z0=50, draw_smith=True):
        self.Z = Z
        self.f = f
        self.Z0 = Z0
        self.draw_smith=draw_smith
        self.components = [('ser', Z)]
        #self.components is an array of tuples containing the history of the circuit,
        #i.e., all components attached to the circuit

        
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
        if Z==None:
            return refl(self.Z,self.Z0)
        else:
            return refl(Z, self.Z0)
        
    def ser(self, Z):
        self.components.append(('ser',Z))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = logspace(log10(Z/1000), log10(Z), 100)
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance

            scaled_value = find_si_prefix(Z)
            label = 'ser {:2g} {:s}$\Omega$'.format(scaled_value[0], scaled_value[1])
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z += Z #after plotting, add series capacitance to the current impedance

    def par(self, Z):
        self.components.append(('par',Z))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = logspace(log10(Z*1000), log10(Z), 100)
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = par(self.Z, tmp)
            #the array transformed into impedance

            scaled_value = find_si_prefix(Z)
            label = 'par {:2g} {:s}$\Omega$'.format(scaled_value[0], scaled_value[1])
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = par(Z, self.Z) #after plotting, add series capacitance to the current impedance    
        
   
    def sercap(self, C):
        '''
        Set a capacitor in series with the present circuit
        '''
        self.components.append(('sercap',C))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.cap(logspace(log10(C*1000), log10(C), 100)) 
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(C)
            label = 'ser {:2g} {:s}F'.format(scaled_value[0], scaled_value[1])
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = ser(self.Z, self.cap(C)) #after plotting, add series capacitance to the current impedance
        
    def serind(self, L):
        '''
        Set an inductor in series with the present circuit
        '''
        self.components.append(('serind',L))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.ind(logspace(log10(L/1000), log10(L), 100)) 
            #impedance of a capacitance array, starting from very small inductance
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(L)
            label = 'ser {:2g} {:s}H'.format(scaled_value[0], scaled_value[1])
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = ser(self.Z, self.ind(L)) #after plotting, add series inductance to the current impedance
        
    def parcap(self, C):
        '''
        Set a capacitor parallel to the present circuit
        '''
        self.components.append(('parcap',C))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.cap(logspace(log10(C/1000), log10(C), 100)) 
            #impedance of a capacitance array, starting from very small cap
            tmp = par(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(C)
            label = 'ser {:2g} {:s}F'.format(scaled_value[0], scaled_value[1])
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = par(self.Z, self.cap(C)) #after plotting, add series capacitance to the current impedance
        
    def parind(self, L):
        '''
        Set an inductor parallel to the present circuit
        '''
        self.components.append(('parind',L))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.ind(logspace(log10(L*1000), log10(L), 100)) 
            #impedance of a capacitance array, starting from very large inductance
            #(large series capacitance == small change in impedance)
            tmp = par(self.Z, tmp)
            #the array transformed into impedance
            scaled_value = find_si_prefix(L)
            label = 'par {:2g} {:s}H'.format(scaled_value[0], scaled_value[1])
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = par(self.Z, self.ind(L)) #after plotting, add series inductance to the current impedance

    def serline(self, kl, Z0=None):
        '''
        perform impedance as seen through a transmission line
        '''
        if Z0==None:
            Z0=self.Z0

        self.components.append(('serline',kl, Z0))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = line(self.Z, linspace(0, kl, 100), Z0) 
            label = 'line {:2g}'.format(kl)
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
            
        self.Z = line(self.Z, kl, Z0) 

    def parline(self, kl, Zl=inf, Z0=None):
        '''
        add a transmission line parallel to the current impedance (defaults to open self.Z0 ohm line)
        '''
        if Z0==None:
            Z0=self.Z0

        self.components.append(('parline',kl, Zl, Z0))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = par(self.Z, line(Zl, linspace(0, kl, 100), Z0))
            label = 'par line {:2g}'.format(kl)
            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network', label=label)
        print line(Zl, kl, Z0) 
        self.Z = par(self.Z0, line(Zl, kl, Z0)) 
        
'''
After this point the file defines basic equations used in Circuit class
'''

def cap(f, C):
    '''
    impedance of a capacitor
    '''
    return 1/(2j*pi*f*C)

def ind(f, L):
    
    return 2j*pi*f*L


def ser(Z1, Z2, Z3=0, Z4=0, Z5=0, Z6=0, Z7=0, Z8=0, Z9=0, Z10=0):
    '''
    impedance of a series of compoents
    '''
    return Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10

def par(Z1, Z2, Z3=inf, Z4=inf, Z5=inf, Z6=inf, Z7=inf, Z8=inf, Z9=inf, Z10=inf):
    '''
    impedance of parallel compoents
    '''
    return 1/(1/Z1+1/Z2+1/Z3+1/Z4+1/Z5+1/Z6+1/Z7+1/Z8+1/Z9+1/Z10)

def refl(Z1, Z0):
    return (Z1-Z0)/(Z1+Z0)
    
def find_si_prefix(x, use_latex = True):
    if use_latex:
        prefixes = ['f', 'p', 'n', '$\mu$', 'm', '', 'k','M','G']
    else:
        prefixes = ['f', 'p', 'n', 'u', 'm', '', 'k','M','G']
    multipliers = asarray([1e-15,1e-12, 1e-9, 1e-6,1e-3,1,1e3,1e6,1e9])
    for n, k in enumerate(x/multipliers):
        if 0.1<k<=(100 + 1e-10):
            break #loop the list until we find a prefix that sets the k to be in nice range
    return (k, prefixes[n])
  
def line(Z, kl, Z0=50):
    if Z!=Inf:
        Z = Z0 * (Z+1j*Z0*tan(kl))/(Z0 + 1j* Z*tan(kl))
    else:
        Z = -1j*Z0*1/tan(kl)
    return Z      

if __name__ == '__main__':
    plt.ion()
    
    c = Circuit(64e6,35)
    #c.sercap(100e-12)
    #c.sercap(110e-12)
    #c.ser(30)
    #c.serind(100e-9)
    c.parcap(100e-12)
    c.par(100)
    c.serind(30e-9)
    #c.parind(30e-9)
    c.serline(0.3)
    c.parline(0.8)
    print c.components
    print c.Z
