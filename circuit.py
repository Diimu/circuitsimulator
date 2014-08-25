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
        
    def ser(self, Z):
        self.components.append(('ser',Z))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = logspace(log10(Z/1000), log10(Z), 100)
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance

            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network')
            
        self.Z += Z #after plotting, add series capacitance to the current impedance

    def par(self, Z):
        self.components.append(('par',Z))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = logspace(log10(Z*1000), log10(Z), 100)
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = par(self.Z, tmp)
            #the array transformed into impedance

            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network')
            
        self.Z = par(Z, self.Z) #after plotting, add series capacitance to the current impedance    
        
        
    def refl(self, Z=None):
        '''
        returns the reflection coefficient of the present circuit
        '''
        if Z==None:
            return refl(self.Z,self.Z0)
        else:
            return refl(Z, self.Z0)
    
   
    def sercap(self, C):
        '''
        Set a capacitor parallel to the present circuit
        '''
        self.components.append(('sercap',C))
        if size(self.f)==1 and self.draw_smith: #single frequency -> plot transition on smith chart
            tmp = self.cap(logspace(log10(C*1000), log10(C), 100)) 
            #impedance of a capacitance array, starting from very large cap
            #(large series capacitance == small change in impedance)
            tmp = ser(self.Z, tmp)
            #the array transformed into impedance

            skrf.plot_smith(self.refl(tmp), x_label='',y_label='', title='Smith chart, matching network')
            
        self.Z += self.cap(C) #after plotting, add series capacitance to the current impedance
        



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
    

if __name__ == '__main__':
    plt.ion()
    
    c = Circuit(64e6,35)
    c.sercap(100e-12)
    c.sercap(100e-12)
    c.ser(30)
    print c.components

