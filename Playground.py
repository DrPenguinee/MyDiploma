import csv
import numpy as np

epoint, Edata, error = np.zeros(301), np.zeros(301), np.zeros(301)

# with open('knm1.dat', 'r') as f, open('res.dat', 'w') as out:
#             reader = csv.reader(f, delimiter=' ')
#             m = 0
#             for row in reader:
#                 epoint[m] = 1000 * float(row[0])
#                 Edata[m] = float(row[1])
#                 error[m] = float(row[2])
#                 print("%3d %10.1f %11.6f %13.8f" % (m, epoint[m], Edata[m], error[m]), file=out)
#                 m+=1

with open('res.dat', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        print(row)