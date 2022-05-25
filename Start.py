import numpy as np
from iminuit import Minuit
import MonteCarlo as MC

values = np.zeros(20)

values[0] = 0.374e1        # e0
values[1] = 0.964e2        # s
values[2] = 0              # snm
values[3] = 0.293          # back
values[4] = 0              # backpar
values[5] = 0.57e03        # step
values[6] = 0              # eend
values[7] = 0              # shift
values[8] = 0              # shift2
values[9] = 0.25           # thick
values[10] = 18535.        # emin
values[11] = 18635.        # emax
values[12] = 0.6e-5        # dtime
values[13] = 1             # ntype
values[14] = 0             # snm2
values[15] = 0             # hnupr
values[16] = 0             # unvis
values[17] = 10000.        # cut
values[18] = 0             # pos_ml
values[19] = 0             # prob_ml

m = Minuit(MC.chi2, tuple(values))
m.errordef = Minuit.LEAST_SQUARES
m.fixed = True
m.fixed[2] = False
m.limits[2] = (-3.95, 20)

sigmas_t = []

for s in sigmas_t:
    MC.file_index = 3
    with open(f"troitsk mnu2 s = {s}.dat", 'w') as g:
        MC.sigma = s
        for i in range(10**5):
            MC.mc_flag = 1
            m.migrad()
            print('{:6d} {:f}'.format(i+1, m.values[2]), file=g)

sigmas_k = [0.2, 0.4, 0.8]

# for s in sigmas_k:
#     MC.file_index = 2
#     with open(f"katrin mnu2 s = {s}.dat", 'w') as g:
#         MC.sigma = s
#         for i in range(10**5):
#             MC.mc_flag = 1
#             m.migrad()
#             print('{:6d} {:f}'.format(i+1, m.values[2]), file=g)


print(m.values)
