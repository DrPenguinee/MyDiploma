from scipy.stats import kstest
from scipy.stats import norm, skew, kurtosis, skewtest, kurtosistest
import numpy as np
import matplotlib.pyplot as plt

files = ['mnu2_mainz.dat', 'mnu2_katrin.dat', 'mnu2_troitsk.dat']

mnu2 = dict.fromkeys(files)
for name in files:
    mnu2[name] = []

file = files[2]
with open(file, 'r') as f:
    reader = f.readlines()
    i = 0
    for row in reader:
        i += 1
        data = list(map(float, row.split()))
        mnu2[file].append(data[1])
        if i == 10000:
            break

mu, sigma = 0, 0.19
# data = np.random.normal(mu, sigma, size=100_000)
data = mnu2[file]
print(kstest(data, "norm", (mu,sigma) ))
print(skewtest(data))
print(kurtosistest(data))

count, bins_count = np.histogram(data, bins=50)
  
# finding the PDF of the histogram using count values
pdf = count / sum(count)
  
# using numpy np.cumsum to calculate the CDF
# We can also find using the PDF values by looping and adding
cdf = np.cumsum(pdf)
  
# plotting PDF and CDF
# plt.plot(bins_count[1:], pdf, color="red", label="PDF")
# plt.plot(bins_count[1:], cdf, "g-o", label="CDF")
plt.step(bins_count[1:], cdf, "g", where='mid', label="CDF")

x = np.linspace(-0.8, 0.8, 100)
plt.plot(x, norm.cdf(x, mu, sigma), label='normal')
plt.grid()
plt.title(file)

info = [f"skew = {skew(data):.4f}\n kurt = {kurtosis(data):.4f}"]
plt.legend(title='\n'.join(info))
plt.show()