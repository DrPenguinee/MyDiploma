from Playground import chi2
import numpy as np
from iminuit import Minuit


m = Minuit(chi2, tuple(np.zeros(3)))
m.errordef = Minuit.LEAST_SQUARES
m.migrad()
print(m.values)