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
    
    def __init__(self, f, Z, Z0=50):
        self.Z = Z
        self.f = f
        self.Z0 = Z0
        
    def cap(self, C):
        '''
        impedance of a capacitor at frequency of this circuit
        '''
        return cap(self.f, C)
        
        
def cap(f, C):
    '''
    impedance of a capacitor
    '''
    return 1/(1j*2*pi*f*C)
