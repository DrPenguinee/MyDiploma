from unittest import result
from scipy.stats import kstest
from scipy.stats import norm, skew, kurtosis, skewtest, kurtosistest
from scipy import stats as st
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import pandas as pd
import texttable
import latextable

files = ['mainz', 'katrin', 'troitsk']

files_k = [0.1, 0.2, 0.25, 0.4, 0.5, 1, 2, 4]
files_k = list(map(str, files_k))

files_m = [0.1, 0.25, 0.5, 1, 2, 4, 6, 7, 8, 9, 10, 11]
files_m = list(map(str, files_m))

files_t = [0.1, 0.25, 0.5, 0.8, 1, 2, 4, 6, 8, 10]
files_t = list(map(str, files_t))

file = files[1]
mnu2 = dict.fromkeys(files + files_k)
mnu2_m = dict.fromkeys(files_m)
mnu2_t = dict.fromkeys(files_t)

for key in files:
    mnu2[key] = []
    with open("mnu2_" + key + ".dat", 'r') as f:
        reader = f.readlines()
        i = 0
        for row in reader:
            i += 1
            data = list(map(float, row.split()))
            mnu2[key].append(data[1])

for key in files_k:
    if key == '1':
        mnu2[key] = mnu2['katrin'].copy()
        continue
    mnu2[key] = []
    with open("katrin mnu2 s = " + key + ".dat", 'r') as f:
        reader = f.readlines()
        i = 0
        for row in reader:
            i += 1
            data = list(map(float, row.split()))
            mnu2[key].append(data[1])

            
for key in files_m:
    if key == '1':
        mnu2_m[key] = mnu2['mainz'].copy()
        continue
    mnu2_m[key] = []
    with open("mainz mnu2 s = " + key + ".dat", 'r') as f:
        reader = f.readlines()
        i = 0
        for row in reader:
            i += 1
            data = list(map(float, row.split()))
            mnu2_m[key].append(data[1])


for key in files_t:
    if key == '1':
        mnu2_t[key] = mnu2['troitsk'].copy()
        continue
    mnu2_t[key] = []
    with open("troitsk mnu2 s = " + key + ".dat", 'r') as f:
        reader = f.readlines()
        i = 0
        for row in reader:
            i += 1
            data = list(map(float, row.split()))
            mnu2_t[key].append(data[1])


results = texttable.Texttable()
results.set_cols_align(["c"] * 6)
results.set_cols_dtype(["t"] * 6)
results.add_row(["Formula", "$m^2_\\nu$", "$\\sigma$", "$D_n$", "$\\alpha$", "skew"])

results_k = texttable.Texttable()
results_k.set_cols_align(["r"] * 6)
results_k.set_cols_dtype(["t"] * 6)
results_k.add_row(["$\\sigma^2_{err}$", "$m^2_\\nu$", "$\\sigma$", "$D_n$", "$\\alpha$", "skew"])

results_m = texttable.Texttable()
results_m.set_cols_align(["r"] * 6)
results_m.set_cols_dtype(["t"] * 6)
results_m.add_row(["$\\sigma^2_{err}$", "$m^2_\\nu$", "$\\sigma$", "$D_n$", "$\\alpha$", "skew"])

results_t = texttable.Texttable()
results_t.set_cols_align(["r"] * 6)
results_t.set_cols_dtype(["t"] * 6)
results_t.add_row(["$\\sigma^2_{err}$", "$m^2_\\nu$", "$\\sigma$", "$D_n$", "$\\alpha$", "skew"])

for name in files:
    for N in [100_000]:
        data = mnu2[name][:N]

        # x = np.sort(data)
        # y = 1. * np.arange(len(data))/(len(data)-1)
        # vals = list(curve_fit(lambda z, mu, sigma: norm.cdf(z, loc=mu, scale=sigma), x, y))
        # vals[1] = np.sqrt(np.diag(vals[1]))

        # mu, sigma = vals[0]
        mu, sigma = np.mean(data), np.std(data)
        skewness = skew(data)
        KS, pvalue = kstest(data, norm.cdf, (mu, sigma))
        if pvalue > 0.001:
            results.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.3f}", f"{skewness:.4f}"])
        else:
            results.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.1e}", f"{skewness:.4f}"])

for name in files_k:
    data = mnu2[name]

    mu, sigma = np.mean(data), np.std(data)
    skewness = skew(data)
    KS, pvalue = kstest(data, norm.cdf, (mu, sigma))
    if pvalue > 0.001:
        results_k.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.3f}", f"{skewness:.4f}"])
    else:
        results_k.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.1e}", f"{skewness:.4f}"])
        
for name in files_m:
    data = mnu2_m[name]

    mu, sigma = np.mean(data), np.std(data)
    skewness = skew(data)
    KS, pvalue = kstest(data, norm.cdf, (mu, sigma))
    if pvalue > 0.001:
        results_m.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.3f}", f"{skewness:.4f}"])
    else:
        results_m.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.1e}", f"{skewness:.4f}"])

for name in files_t:
    data = mnu2_t[name]

    mu, sigma = np.mean(data), np.std(data)
    skewness = skew(data)
    KS, pvalue = kstest(data, norm.cdf, (mu, sigma))
    if pvalue > 0.001:
        results_t.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.3f}", f"{skewness:.4f}"])
    else:
        results_t.add_row([name, f"{mu:.2e}", f"{sigma:.3f}", f"{KS:.5f}", f"{pvalue:.1e}", f"{skewness:.4f}"])

with open("results of fitting.dat", "w") as f:
    print(results.draw(), file=f)
    print('\n katrin', file=f)
    print(results_k.draw(), file=f)
    print('\n mainz', file=f)
    print(results_m.draw(), file=f)
    print('\n troitsk', file=f)
    print(results_t.draw(), file=f)

with open("latex results.dat", "w") as f:
    print(latextable.draw_latex(results), file=f)
    print('\n katrin', file=f)
    print(latextable.draw_latex(results_k), file=f)
    print('\n mainz', file=f)
    print(latextable.draw_latex(results_m), file=f)
    print('\n troitsk', file=f)
    print(latextable.draw_latex(results_t), file=f)
"""
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
"""