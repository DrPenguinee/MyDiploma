import numpy as np
from iminuit import Minuit
from Diploma import chi2

values = np.zeros(20)

values[0] = 0.374e1
values[1] = 0
values[2] = 0
values[3] = 0.293
values[4] = 0
values[5] = 0.57e03
values[6] = 0
values[7] = 0
values[8] = 0
values[9] = 0.25
values[10] = 18535.
values[11] = 18635.
values[12] = 0.6e-5
values[13] = 1
values[14] = 0
values[15] = 0
values[16] = 0
values[17] = 10000.
values[18] = 0
values[19] = 0

m = Minuit(chi2, tuple(values))
m.errordef = Minuit.LEAST_SQUARES
m.limits = [(0, None) for i in range(20)]
m.fixed = True
m.fixed[1] = False

m.migrad()
print(m.values)