import numpy as np
import iminuit
from iminuit import Minuit

def chi2(a, b):
    return a**2 + (b-3)**2

m = Minuit(chi2, a=1, b=1)
m.migrad()
