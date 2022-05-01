from scipy import interpolate
import numpy as np

x = np.array([1., 2., 3., 4.])
y = np.array([1., 4., 9., 16.])

x = x[::-1]
y = y[::-1]

sp = interpolate.interp1d(x, y, kind='cubic')
print(sp(0.5))

