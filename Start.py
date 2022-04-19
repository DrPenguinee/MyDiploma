import numpy as np
from iminuit import Minuit
from Diploma import chi2

m = Minuit(chi2, tuple(np.zeros(20)))
m.errordef = Minuit.LEAST_SQUARES

m.migrad()
print(m.values)