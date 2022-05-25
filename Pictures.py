from cProfile import label
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

files = ['mnu2_mainz.dat', 'mnu2_katrin.dat', 'mnu2_troitsk.dat']
names = ['Майнц', 'КАТРИН', 'Троицк']

index = 0
file = files[index]
title = names[index]

mnu2 = dict.fromkeys(files)
for name in files:
    mnu2[name] = []

for file in files:
    with open(file, 'r') as f:
        reader = f.readlines()
        for row in reader:
            data = list(map(float, row.split()))
            mnu2[file].append(data[1])


plot1 = plt.subplot2grid((3, 1), (0, 0), colspan=1)
plot2 = plt.subplot2grid((3, 1), (1, 0), colspan=1)
plot3 = plt.subplot2grid((3, 1), (2, 0), colspan=1)

i = 0
data = mnu2[file]
for plot in (plot1, plot2, plot3):
    data = mnu2[files[i]]
    mu, sigma = np.mean(data), np.std(data)
    gauss = np.random.normal(mu, sigma, 1_000_000)
    plot.hist(gauss, bins=25, density=True, label="Normal", color="orange")
    plot.hist(data, bins=25, density=True, label=names[i])
    plot.grid()
    plot.xaxis.set_major_locator(MultipleLocator(base=0.25))
    plot.yaxis.set_major_locator(MultipleLocator(base=0.5))
    plot.set_xlabel("$m^2_\\nu$")
    plot.set_ylabel("N")
    plot.legend()
    # plot.set_title(names[i])
    i += 1


plt.tight_layout()
plt.show()
        