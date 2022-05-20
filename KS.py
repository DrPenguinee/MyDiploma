from scipy.stats import kstest
from scipy.stats import norm, skew, kurtosis, skewtest, kurtosistest
from scipy import stats as st
import numpy as np
import matplotlib.pyplot as plt

N = 100_000

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
        if i == N:
            break

mu, sigma = 0, 0.18
# data = np.random.normal(mu, sigma, size=100_000)
data = mnu2[file]
dist = getattr(st, "norm")        
param = dist.fit(data)
mu, sigma = param
print(param)

count, bins_count = np.histogram(data, bins=100)
pdf = count / sum(count)  
cdf = np.cumsum(pdf)
  

KStest, p = kstest(data, "norm", param )
print(skewtest(data))
print(kurtosistest(data))
# plotting PDF and CDF
# plt.plot(bins_count[1:], pdf, color="red", label="PDF")
# plt.plot(bins_count[1:], cdf, "g-o", label="CDF")
plt.step(bins_count[1:], cdf, "g", where='mid', label="CDF of data")

x = np.linspace(-0.8, 0.8, 100)
plt.plot(x, norm.cdf(x, mu, sigma), label=f'normal, \n$\sigma$ = {sigma} \n$\mu$ = {mu}')
plt.grid()
plt.title(file + f"  N = {N}")

info = [f"skew = {skew(data):.4f}\n kurt = {kurtosis(data, fisher=False):.4f}\n KStest = {KStest:.4f}\n pvalue = {p:.4f}"]
plt.legend(title='\n'.join(info))
plt.show()