import matplotlib.pyplot as plt
import numpy as np

files = ['mnu2_mainz.dat', 'mnu2_katrin.dat', 'mnu2_troitsk.dat']

mnu2 = dict.fromkeys(files)
for name in files:
    mnu2[name] = []

file = files[2]
with open(file, 'r') as f:
    reader = f.readlines()
    for row in reader:
        data = list(map(float, row.split()))
        mnu2[file].append(data[1])

plt.hist(mnu2[file])
plt.grid()
plt.show()
        