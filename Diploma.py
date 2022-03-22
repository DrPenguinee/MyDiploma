from tkinter import E
import numpy as np
import csv

em = 510999.06
alpha = 1./137.0359895
thick = .0
shift = shift2 = .0
epoint = npoint = .0
IX, IY, IZ = 1, 255, 25555

def chi2(npar, grad, fval, xval, iflag):
    epoint = Edata = error = espnm = np.zeros(301)
    espec = backval = endef = eml = np.zeros(301)
    trspec = hnu = espnm2 = fakesp = np.zeros(301)
    rancor = tspec = np.zeros(301)
    ifl = 0

    e0 = xval[0]+18570.0                # e0 - spectrum endpoint
    s = xval[1]                         # s - spectrum square
    snm = xval[2]                       # snm - squared neutrino mass
    back = xval[3]                      # back - background value
    backpar = xval[4]                   # backpar - background parameter
    step = xval[5]                      # step - value of endpoint effect
    eend = xval[6] + 18000.0            # eend - endpoint of step effect
    shift = xval[7]                     # shift - probability shift of excitation-ionization losess
    shift2 = xval[8]                    # shift2 - mainstate probability shift
    thick = xval[9]                     # thick - source thickness factor
    emin = xval[10]                     # emin - lower limit of analysis interval
    emax = xval[11]                     # emax - upper limit of analysis interval
    dtime = xval[12]                    # dtime - counting system deadtime
    dtime0 = 6.e-6
    ntype = int(xval[13] + 0.1)         # ntype = 1 - calculation by splines
                                        # ntype = 2 - direct calculation
    snm2 = xval[14]**2                  # snm2 - squared heavy neutrino mass
    hnupr = xval[15]                    # hnupr - heavy neutrino probability
    unvis = xval[16]                    # unvis - unvisible correction shift
    cut = xval[17]                      # cut - point "glitches" cut
    pos_ml = e0 - xval[18]              # pos_ml - missing level energy
    prob_ml = xval[19]                  # prob_ml - missing level probability
    acc = .0002

    if ifl == 0:
        ifl += 1
        
        Num_MC = 0                           # MC simulations counter

        # data file read --------------------------------------------------------
        npoint = 27

        with open('knm1.dat', 'r') as f:
            reader = csv.reader(f, delimiter=' ')
            m = 0
            for row in reader:
                epoint[m] = 1000 * float(row[0])
                Edata[m] = float(row[1])
                error[m] = float(row[2])
                m+=1
                print("%3d %10.1f %11.6f %13.8f\n" % (m, epoint[m], Edata[m], error[m]))

        
        # Arrays calculation
        for i in range(1, 101):
            entr = float(i)
            tran = transmission(entr)
            # fstail(e0 - entr - 165., e0, tail)
            specint(e0 - entr, vesp, e0, snm, vespnm, thick)
            print(entr, tran, tail)
        # Arrays calculated

        # Faked spectrum formation--------------------------------------
        w = 96.4
        with open('run_kat_MC2.dat', 'w') as f:
            print("Simulation with Elow=", Emin, file=f)
            print("%s %s %s %s" % ("HV", "Freq", "Err", "Time"))
            for i in range(1, npoint+1):
                e = epoint[i]
                specint(e, espec[i], e0, snm, espnm[i], thick)
                specint(e, hnu[i], e0, snm2, espnm2[i], thick)
                specmlint(e, eml[i], pos_ml, e0)
                trapbackground(e, trspec[i], e0)
                endeffect(e, eend, step, endef[i])
                background (e, e0, back, backpar, backval[i])
                backval[i] = backval[i]
                fakesp[i] = w * (espec[i] + hnu[i]*hnupr + eml[i]*prob_ml) + backval[i] + endef[i]
                Edata[i] = fakesp[i]
                print(e, Edata[i], error[i])
        # Faked spectrum formed--------------------------------------------

    if iflag == 5:
        Num_MC += 1
        for i in range(1, npoint+1):
            Edata[i] = fakesp[i] + error[i]*np.random.normal()
            rancor[i] = (Edata[i] - fakesp[i])/error[i]
        for i in range(1, npoint+1):
            print(rancor[i])
        print(Num_MC)

    # string 188
    if iflag == 3:
        pass

    fval = 0.
    sespec = 0.
    sml = 0.
    strap = 0.
    sbackval = 0.
    sendef = 0.
    for i in range(1, npoint+1):
        e = epoint[i]
        if iflag != 3:
            if e < emin or e > emax:
                continue
        if ntype == 1:
            spectint(e, espec[i], e0, snm, espnm[i], thick)
            spectint(e, hnu[i], e0, snm, espnm[i], thick)
        if ntype == 2:
            expspectrum(e, espec[i], e0, snm, ACC)
            expspectrum(e, hnu[i], e0, snm2, ACC)
        specmlint(e, eml[i], pos_ml, e0)
        background(e, e0, back, backpar, backval[i])
        endeffect(e, eend, step, endef[i])
        # if
        sespec = sespec + espec[i] + hnu[i]*hnupr
        sml = sml +eml[i]*prob_ml
        sbackval = sbackval + backval[i]
        sendef = sendef + endef[i]
    w = (s - sbackval)/sespec
    for i in range(1, npoint+1):
        e = epoint[i]
        theor = w * (espec[i] + hnu[i]*hnupr + trspec[i]*backpar + eml[i]*prob_ml) \
            + backval[i] + endef[i]
        exper = Edata[i] * (1 - Edata[i]*dtime0)/(1 - Edata[i]*dtime)
        if iflag == 3:
            with open('graph_kat.dat', 'w') as f:
                f.write(
                    i, e-18000.0, 
                    (exper-w*(espec[i]-espnm[i]+hnupr*(hnu[i]-espnm2[i]))-backval[i])/error[i],
                    1.0, endef[i]/error[i], w*espnm[i]/error[i], w * espnm2[i]*hnupr/error[i]                
                )
            with open('graph_katd.dat', 'w') as f:
                f.write(
                    i, e-18000.0, exper, w * (espec[i]-espnm[i]+hnupr*(hnu[i]-espnm2[i])) - backval[i],
                    endef[i], w * espnm[i], w * espnm2[i] *hnupr, w * eml[i] * prob_ml
                )

        if e < (emin - 0.01) or e > emax:
            continue
        dif = (exper - theor)/error[i]
        if abs(dif) > cut:
            continue
        fval += dif ** 2
    return

    
def specmlint(e, eml, e0, emax):
    pass

def expmlspectrum(E, Specml, E0, ac):
    pass

def specint(E, ESP, E0, SNM, espnm, thick):
    pass



def transmission(entr):
    return
    #
    #



def endeffect(e, eend, step, endef):
    energy = eend - e
    transmission(energy, tran)
    endef = tran * step
    return

def background(e, e0, back, backpar, backval):
    backval = back + backpar * (18575.0 - e)/40.0
    return

